from biosimulators_smoldyn.utils import get_parameters_variables_outputs_for_simulation
from biosimulators_smoldyn.combine import validate_variables


model_fp = 'biosimulators_smoldyn/Min1/Andrews-Min1.txt'

model_lang = 'urn:sedml:language:smoldyn'


parameters = get_parameters_variables_outputs_for_simulation(model_fp, model_lang)

variables = parameters[2]

validation = validate_variables(variables)

print(validation)

