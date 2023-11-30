from biosimulators_smoldyn.utils import get_parameters_variables_outputs_for_simulation


model_fp = 'biosimulators_smoldyn/Andrews-Min1.txt'

model_lang = 'urn:sedml:language:smoldyn'


parameters = get_parameters_variables_outputs_for_simulation(model_fp, model_lang)

variables = parameters[2]

