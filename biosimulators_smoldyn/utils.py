''' Methods for parsing configurations of Smoldyn simulations for use with SED-ML
'''

__author__ = 'Jonathan Karr'
__email__ = 'karr@mssm.edu'

__all__ = [
    'read_simulation',
]


import os
import re
import types
from typing import Optional, List
from biosimulators_smoldyn.data_model import Simulation as smoldynSim, SimulationInstruction
from biosimulators_utils.model_lang.smoldyn.validation import validate_model
from biosimulators_utils.config import Config
from biosimulators_utils.sedml.data_model import (
    SedDocument,
    ModelAttributeChange,
    Variable,
    Symbol,
    Simulation,
    UniformTimeCourseSimulation,
    Algorithm,
    Task
)
from biosimulators_utils.utils.core import flatten_nested_list_of_strings


def get_parameters_variables_outputs_for_simulation(
        model_filename: str,
        model_language: str,
        simulation_type: Simulation,
        algorithm_kisao_id=None,
        change_level=SedDocument,
        native_ids=False,
        native_data_types=False,
        config=None
        ):
    """ Get the possible observables for a simulation of a model

    Args:
        model_filename (:obj:`str`): path to model file
        model_language (:obj:`str`): model language (e.g., ``urn:sedml:language:sbml``)
        simulation_type (:obj:`types.Type`): subclass of :obj:`Simulation`
        algorithm_kisao_id (:obj:`str`): KiSAO id of the algorithm for simulating the model (e.g., ``KISAO_0000019``
            for CVODE)
        change_level (:obj:`types.Type`, optional): level at which model changes will be made (:obj:`SedDocument` or :obj:`Task`)
        native_ids (:obj:`bool`, optional): whether to return the raw id and name of each model component rather than the suggested name
            for the variable of an associated SED-ML data generator
        native_data_types (:obj:`bool`, optional): whether to return new_values in their native data types
        config (:obj:`Config`, optional): whether to fail on missing includes

    Returns:
        :obj:`list` of :obj:`ModelAttributeChange`: possible attributes of a model that can be changed and their default values
        :obj:`list` of :obj:`Simulation`: simulation of the model
        :obj:`list` of :obj:`Variable`: possible observables for a simulation of the model
    """
    # check model file exists and is valid
    if not isinstance(model_filename, str):
        raise ValueError('`{}` is not a path to a model file.'.format(model_filename))

    if not os.path.isfile(model_filename):
        raise FileNotFoundError('Model file `{}` does not exist.'.format(model_filename))

    errors, _, (smoldyn_model, model_config) = validate_model(model_filename, config=config)
    if errors:
        raise ValueError('Model file `{}` is not a valid Smoldyn file.\n  {}'.format(
            model_filename, flatten_nested_list_of_strings(errors).replace('\n', '\n  ')))

    if simulation_type not in [UniformTimeCourseSimulation]:
        raise NotImplementedError('`simulation_type` must be `OneStepSimulation` or `UniformTimeCourseSimulation`')

    sim = UniformTimeCourseSimulation(
        id='simulation',
        initial_time=smoldyn_model.start,
        output_start_time=smoldyn_model.start,
        output_end_time=smoldyn_model.stop,
        number_of_steps=int((smoldyn_model.stop - smoldyn_model.start) / smoldyn_model.dt),
        algorithm=Algorithm(
            kisao_id=algorithm_kisao_id or 'KISAO_0000057',
        ),
    )

    # get parameters and observables
    model = read_simulation(model_filename)

    params = []
    for instruction in model.instructions:
        if not (
            instruction.macro.startswith('define ')
            # or instruction.macro.startswith('define_global ')
            or instruction.macro.startswith('difc ')
            or instruction.macro.startswith('difc_rule ')
            or instruction.macro.startswith('difm ')
            or instruction.macro.startswith('difm_rule ')
            or instruction.macro.startswith('drift ')
            or instruction.macro.startswith('drift_rule ')
            or instruction.macro.startswith('surface_drift ')
            or instruction.macro.startswith('surface_drift_rule ')
        ):
            continue

        if native_data_types:
            id = instruction.macro.partition(' ')[2]
        else:
            id = re.sub(r'[^a-zA-Z0-9_]', '_', instruction.macro)

        if native_data_types:
            if (
                instruction.macro.startswith('define ')
                # or instruction.macro.startswith('define_global ')
                or instruction.macro.startswith('difc ')
                or instruction.macro.startswith('difc_rule ')
            ):
                new_value = float(instruction.arguments)
            elif (
                instruction.macro.startswith('difm ')
                or instruction.macro.startswith('difm_rule ')
                or instruction.macro.startswith('drift ')
                or instruction.macro.startswith('drift_rule ')
                or instruction.macro.startswith('surface_drift ')
                or instruction.macro.startswith('surface_drift_rule ')
            ):
                new_value = [float(val) for val in instruction.arguments.split(' ')]
        else:
            new_value = instruction.arguments

        if instruction.macro.partition(' ')[2] != 'all':
            params.append(ModelAttributeChange(
                id=id,
                name=None if native_ids else instruction.description,
                target=instruction.macro,
                new_value=new_value,
            ))

    smoldyn_model.addOutputData('counts')
    smoldyn_model.addCommand(cmd='molcount counts', cmd_type='E')
    for compartment in model.compartments:
        data_id = 'counts_cmpt_' + compartment
        smoldyn_model.addOutputData(data_id)
        smoldyn_model.addCommand(cmd='molcountincmpt ' + compartment + ' ' + data_id, cmd_type='E')
    for surface in model.surfaces:
        data_id = 'counts_surf_' + surface
        smoldyn_model.addOutputData(data_id)
        smoldyn_model.addCommand(cmd='molcountonsurf ' + surface + ' ' + data_id, cmd_type='E')

    smoldyn_model.run(stop=1e-12, dt=1., overwrite=True, display=False, quit_at_end=False)

    data_id = 'counts'
    species_counts = smoldyn_model.getOutputData(data_id, True)[0][1:]
    for species, count in zip(model.species, species_counts):
        params.append(ModelAttributeChange(
            id=species if native_ids else 'initial_count_species_{}'.format(re.sub('[^a-zA-Z0-9_]', '_', species)),
            name=None if native_ids else 'Initial count of species "{}"'.format(species),
            target="fixmolcount {}".format(species),
            new_value=count if native_data_types else str(count),
        ))

    for compartment in model.compartments:
        data_id = 'counts_cmpt_' + compartment
        species_counts = smoldyn_model.getOutputData(data_id, True)[0][1:]
        for species, count in zip(model.species, species_counts):
            params.append(ModelAttributeChange(
                id="{}.{}".format(species, compartment) if native_ids else 'initial_count_species_{}_compartment_{}'.format(
                    re.sub('[^a-zA-Z0-9_]', '_', species), compartment),
                name=None if native_ids else 'Initial count of species "{}" in compartment "{}"'.format(species, compartment),
                target="fixmolcountincmpt {} {}".format(species, compartment),
                new_value=new_value,
            ))

    for surface in model.surfaces:
        data_id = 'counts_surf_' + surface
        species_counts = smoldyn_model.getOutputData(data_id, True)[0][1:]
        for species, count in zip(model.species, species_counts):
            params.append(ModelAttributeChange(
                id="{}.{}".format(species, surface) if native_ids else 'initial_count_species_{}_surface_{}'.format(
                    re.sub('[^a-zA-Z0-9_]', '_', species), surface),
                name=None if native_ids else 'Initial count of species "{}" in surface "{}"'.format(species, surface),
                target="fixmolcountonsurf {} {}".format(species, surface),
                new_value=new_value,
            ))

    vars = []

    # Add global commands (not specific to a given species)
    vars.append(Variable(
        id='execution_time' if native_ids else 'execution_time',
        name='Execution time',
        target='executiontime',
    ))

    vars.append(Variable(
        id='molecule_list',
        name='Molecule list',
        target='listmols',
    ))

    '''vars.append(Variable(
        id=None if native_ids else 'time',
        name=None if native_ids else 'Time',
        symbol=Symbol.time.value,
    ))'''
    for species in model.species:
        vars.append(Variable(
            id=species if native_ids else 'count_species_{}'.format(re.sub('[^a-zA-Z0-9_]', '_', species)),
            name=None if native_ids else 'Count of species "{}"'.format(species),
            target="molcount {}".format(species),
        ))
        for compartment in model.compartments:
            vars.append(Variable(
                id="{}.{}".format(species, compartment) if native_ids else 'count_species_{}_compartment_{}'.format(
                    re.sub('[^a-zA-Z0-9_]', '_', species), compartment),
                name=None if native_ids else 'Count of species "{}" in compartment "{}"'.format(species, compartment),
                target="molcountincmpt {} {}".format(species, compartment),
            ))
        for surface in model.surfaces:
            vars.append(Variable(
                id="{}.{}".format(species, surface) if native_ids else 'count_species_{}_surface_{}'.format(
                    re.sub('[^a-zA-Z0-9_]', '_', species), surface),
                name=None if native_ids else 'Count of species "{}" in surface "{}"'.format(species, surface),
                target="molcountonsurf {} {}".format(species, surface),
            ))

    return (params, [sim], vars, [])


def read_simulation(filename):
    """ Read the configuration for a Smoldyn simulation

    Args:
        filename (:obj:`str`): path to a configuration for a Smoldyn simulation

    Returns:
        :obj:`Simulation`: data structure which represents the configuration of the Smoldyn simulation
    """
    sim = Simulation()
    param_group_counts = {}
    with open(filename, 'r') as file:
        for line in file:
            # separate comments
            line, _, comment = line.partition('#')
            line = line.strip()
            comment = comment.strip()

            # remove consecutive spaces
            line = re.sub(' +', ' ', line)

            _read_simulation_line(line, param_group_counts, sim)

            if re.match(r'^end_file( |$)', line):
                break

    return sim


def _read_simulation_line(line, macro_counts, sim):
    """ Parse a line of a configuration of a Smoldyn simulation

    Args:
        line (:obj:`str`): line of a configuration of a Smoldyn simulation
        macro_counts (:obj:`dict`): dictionary used to count instances of the same macro
        sim (:obj:`Simulation`): data structure which represents the configuration of the Smoldyn simulation
    """
    if line.startswith('species '):
        sim.species.extend(line.split(' ')[1:])

    elif line.startswith('start_compartment '):
        sim.compartments.append(line.partition(' ')[2])

    elif line.startswith('start_surface '):
        sim.surfaces.append(line.partition(' ')[2])

    for pattern in CONFIG_DECLARATION_PATTERNS:
        match = re.match(pattern['regex'], line)
        if match:
            if pattern.get('macro', None):
                group = pattern['macro']['group'](match)
                if group not in macro_counts:
                    macro_counts[group] = -1
                macro_counts[group] += 1
                i_group = macro_counts[group]

                sim.instructions.append(SimulationInstruction(
                    id=pattern['macro']['id'](match, i_group),
                    description=pattern['macro']['description'](match, i_group),
                    macro=pattern['macro']['macro'](match),
                    arguments=pattern['macro']['arguments'](match),
                ))


CONFIG_DECLARATION_PATTERNS = [
    {
        'regex': r'^(dim) (.*?)$',
        'macro': {
            'group': lambda match: 'number_dimensions',
            'id': lambda match, i_group: 'number_dimensions',
            'description': lambda match, i_group: 'Number of dimensions',
            'macro': lambda match: match.group(1),
            'arguments': lambda match: match.group(2),
        }
    },
    {
        'regex': r'^low_wall ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'low_wall_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'low_{}_wall_{}'.format(match.group(1), i_group + 1),
            'description': lambda match, i_group: 'Low {} wall {}'.format(match.group(1), i_group + 1),
            'macro': lambda match: 'low_wall {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^high_wall ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'high_wall_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'high_{}_wall_{}'.format(match.group(1), i_group + 1),
            'description': lambda match, i_group: 'High {} wall {}'.format(match.group(1), i_group + 1),
            'macro': lambda match: 'high_wall {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^boundaries ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'boundaries_{}'.format(match.group(1)),
            'id': lambda match, i_group: '{}_boundary'.format(match.group(1)),
            'description': lambda match, i_group: '{} boundary'.format(match.group(1).upper()),
            'macro': lambda match: 'boundaries {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^define ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'define_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'value_parameter_{}'.format(match.group(1)),
            'description': lambda match, i_group: 'Value of parameter "{}"'.format(match.group(1)),
            'macro': lambda match: 'define {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^difc ([^ \(\)]+) (.*?)$',
        'macro': {
            'group': lambda match: 'difc_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'diffusion_coefficient_species_{}'.format(match.group(1)),
            'description': lambda match, i_group: 'Diffusion coefficient of species "{}"'.format(match.group(1)),
            'macro': lambda match: 'difc {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^difc ([^ \(\)]+)\(([^ \(\)]+)\) (.*?)$',
        'macro': {
            'group': lambda match: 'difc_{}_{}'.format(match.group(1), match.group(2)),
            'id': lambda match, i_group: 'diffusion_coefficient_species_{}_state_{}'.format(match.group(1), match.group(2)),
            'description': lambda match, i_group: 'Diffusion coefficient of species "{}" in state "{}"'.format(match.group(1), match.group(2)),
            'macro': lambda match: 'difc {}({})'.format(match.group(1), match.group(2)),
            'arguments': lambda match: match.group(3),
        },
    },
    {
        'regex': r'^difc_rule ([^ \(\)]+) (.*?)$',
        'macro': {
            'group': lambda match: 'difc_rule_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'diffusion_coefficient_rule_species_{}'.format(re.sub('[^a-zA-Z0-9_]', '_', match.group(1))),
            'description': lambda match, i_group: 'Diffusion coefficient rule for species "{}"'.format(match.group(1)),
            'macro': lambda match: 'difc_rule {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^difc_rule ([^ \(\)]+)\(([^ \(\)]+)\) (.*?)$',
        'macro': {
            'group': lambda match: 'difc_rule_{}_{}'.format(
                match.group(1), match.group(2)),
            'id': lambda match, i_group: 'diffusion_coefficient_rule_species_{}_state_{}'.format(
                re.sub('[^a-zA-Z0-9_]', '_', match.group(1)), match.group(2)),
            'description': lambda match, i_group: 'Diffusion coefficient rule for species "{}" in state "{}"'.format(
                match.group(1), match.group(2)),
            'macro': lambda match: 'difc_rule {}({})'.format(
                match.group(1), match.group(2)),
            'arguments': lambda match: match.group(3),
        },
    },
    {
        'regex': r'^difm ([^ \(\)]+) (.*?)$',
        'macro': {
            'group': lambda match: 'difm_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'membrane_diffusion_coefficient_species_{}'.format(match.group(1)),
            'description': lambda match, i_group: 'Membrane diffusion coefficient of species "{}"'.format(match.group(1)),
            'macro': lambda match: 'difm {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^difm ([^ \(\)]+)\(([^ \(\)]+)\) (.*?)$',
        'macro': {
            'group': lambda match: 'difm_{}_{}'.format(match.group(1), match.group(2)),
            'id': lambda match, i_group: 'membrane_diffusion_coefficient_species_{}_state_{}'.format(
                match.group(1), match.group(2)),
            'description': lambda match, i_group: 'Membrane diffusion coefficient of species "{}" in state "{}"'.format(
                match.group(1), match.group(2)),
            'macro': lambda match: 'difm {}({})'.format(match.group(1), match.group(2)),
            'arguments': lambda match: match.group(3),
        },
    },
    {
        'regex': r'^difm_rule ([^ \(\)]+) (.*?)$',
        'macro': {
            'group': lambda match: 'difm_rule_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'membrane_diffusion_coefficient_rule_species_{}'.format(
                re.sub('[^a-zA-Z0-9_]', '_', match.group(1))),
            'description': lambda match, i_group: 'Membrane diffusion coefficient rule for species "{}"'.format(match.group(1)),
            'macro': lambda match: 'difm_rule {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^difm_rule ([^ \(\)]+)\(([^ \(\)]+)\) (.*?)$',
        'macro': {
            'group': lambda match: 'difm_rule_{}_{}'.format(match.group(1), match.group(2)),
            'id': lambda match, i_group: 'membrane_diffusion_coefficient_rule_species_{}_state_{}'.format(
                re.sub('[^a-zA-Z0-9_]', '_', match.group(1)), match.group(2)),
            'description': lambda match, i_group: 'Membrane diffusion coefficient rule for species "{}" in state "{}"'.format(
                match.group(1), match.group(2)),
            'macro': lambda match: 'difm_rule {}({})'.format(match.group(1), match.group(2)),
            'arguments': lambda match: match.group(3),
        },
    },
    {
        'regex': r'^drift ([^ \(\)]+) (.*?)$',
        'macro': {
            'group': lambda match: 'drift_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'drift_species_{}'.format(match.group(1)),
            'description': lambda match, i_group: 'Drift of species "{}"'.format(match.group(1)),
            'macro': lambda match: 'drift {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^drift ([^ \(\)]+)\(([^ \(\)]+)\) (.*?)$',
        'macro': {
            'group': lambda match: 'drift_{}_{}'.format(match.group(1), match.group(2)),
            'id': lambda match, i_group: 'drift_species_{}_state_{}'.format(
                match.group(1), match.group(2)),
            'description': lambda match, i_group: 'Drift of species "{}" in state "{}"'.format(
                match.group(1), match.group(2)),
            'macro': lambda match: 'drift {}({})'.format(match.group(1), match.group(2)),
            'arguments': lambda match: match.group(3),
        },
    },
    {
        'regex': r'^drift_rule ([^ \(\)]+) (.*?)$',
        'macro': {
            'group': lambda match: 'drift_rule_{}'.format(match.group(1)),
            'id': lambda match, i_group: 'drift_rule_species_{}'.format(re.sub('[^a-zA-Z0-9_]', '_', match.group(1))),
            'description': lambda match, i_group: 'Drift rule for species "{}"'.format(match.group(1)),
            'macro': lambda match: 'drift_rule {}'.format(match.group(1)),
            'arguments': lambda match: match.group(2),
        },
    },
    {
        'regex': r'^drift_rule ([^ \(\)]+)\(([^ \(\)]+)\) (.*?)$',
        'macro': {
            'group': lambda match: 'drift_rule_{}_{}'.format(match.group(1), match.group(2)),
            'id': lambda match, i_group: 'drift_rule_species_{}_state_{}'.format(
                re.sub('[^a-zA-Z0-9_]', '_', match.group(1)), match.group(2)),
            'description': lambda match, i_group: 'Drift rule for species "{}" in state "{}"'.format(
                match.group(1), match.group(2)),
            'macro': lambda match: 'drift_rule {}({})'.format(match.group(1), match.group(2)),
            'arguments': lambda match: match.group(3),
        },
    },
    {
        'regex': r'^surface_drift ([^ \(\)]+) ([^ ]+) ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'surface_drift_{}_{}_{}'.format(
                match.group(1), match.group(2), match.group(3)),
            'id': lambda match, i_group: 'surface_drift_species_{}_surface_{}_shape_{}'.format(
                match.group(1), match.group(2), match.group(3)),
            'description': lambda match, i_group: 'Surface drift of species "{}" on surface "{}" with panel shape "{}"'.format(
                match.group(1), match.group(2), match.group(3)),
            'macro': lambda match: 'surface_drift {} {} {}'.format(
                match.group(1), match.group(2), match.group(3)),
            'arguments': lambda match: match.group(4),
        },
    },
    {
        'regex': r'^surface_drift ([^ \(\)]+)\(([^ \(\)]+)\) ([^ ]+) ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'surface_drift_{}_{}_{}_{}'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'id': lambda match, i_group: 'surface_drift_species_{}_state_{}_surface_{}_shape_{}'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'description': lambda match, i_group: 'Surface drift of species "{}" in state "{}" on surface "{}" with panel shape "{}"'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'macro': lambda match: 'surface_drift {}({}) {} {}'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'arguments': lambda match: match.group(5),
        },
    },
    {
        'regex': r'^surface_drift_rule ([^ \(\)]+) ([^ ]+) ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'surface_drift_rule_{}_{}_{}'.format(
                match.group(1), match.group(2), match.group(3)),
            'id': lambda match, i_group: 'surface_drift_rule_species_{}_surface_{}_panel_{}'.format(
                re.sub('[^a-zA-Z0-9_]', '_', match.group(1)), match.group(2), match.group(3)),
            'description': lambda match, i_group: 'Surface drift rule for species "{}" on surface "{}" of panel shape "{}"'.format(
                match.group(1), match.group(2), match.group(3)),
            'macro': lambda match: 'surface_drift_rule {} {} {}'.format(
                match.group(1), match.group(2), match.group(3)),
            'arguments': lambda match: match.group(4),
        },
    },
    {
        'regex': r'^surface_drift_rule ([^ \(\)]+)\(([^ \(\)]+)\) ([^ ]+) ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'surface_drift_rule_{}_{}_{}_{}'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'id': lambda match, i_group: 'surface_drift_rule_species_{}_state_{}_surface_{}_panel_{}'.format(
                re.sub('[^a-zA-Z0-9_]', '_', match.group(1)), match.group(2), match.group(3), match.group(4)),
            'description': lambda match, i_group: 'Surface drift rule for species "{}" in state "{}" on surface "{}" of panel shape "{}"'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'macro': lambda match: 'surface_drift_rule {}({}) {} {}'.format(
                match.group(1), match.group(2), match.group(3), match.group(4)),
            'arguments': lambda match: match.group(5),
        },
    },
    {
        'regex': r'^mol ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'mol_{}'.format(match.group(2)),
            'id': lambda match, i_group: 'initial_count_species_{}'.format(re.sub(r'[^a-zA-Z0-9_]', '_', match.group(2))),
            'description': lambda match, i_group: 'Initial count of species "{}"'.format(match.group(2)),
            'macro': lambda match: 'mol {}'.format(match.group(2)),
            'arguments': lambda match: match.group(1),
        },
    },
    {
        'regex': r'^compartment_mol ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'compartment_mol_{}'.format(match.group(2)),
            'id': lambda match, i_group: 'initial_count_species_{}'.format(re.sub(r'[^a-zA-Z0-9_]', '_', match.group(2))),
            'description': lambda match, i_group: 'Initial count of species "{}"'.format(match.group(2)),
            'macro': lambda match: 'compartment_mol {}'.format(match.group(2)),
            'arguments': lambda match: match.group(1),
        },
    },
    {
        'regex': r'^surface_mol ([^ ]+) (.*?)$',
        'macro': {
            'group': lambda match: 'surface_mol_{}'.format(match.group(2)),
            'id': lambda match, i_group: 'initial_count_species_{}'.format(re.sub(r'[^a-zA-Z0-9_]', '_', match.group(2))),
            'description': lambda match, i_group: 'Initial count of species "{}"'.format(match.group(2)),
            'macro': lambda match: 'surface_mol {}'.format(match.group(2)),
            'arguments': lambda match: match.group(1),
        },
    },
]
