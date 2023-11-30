from biosimulators_smoldyn.utils import get_parameters_variables_outputs_for_simulation
from biosimulators_smoldyn.combine import validate_variables
from biosimulators_smoldyn.combine import exec_sed_doc
from biosimulators_utils.report.data_model import ReportFormat
from biosimulators_utils.config import Config
import os


model_fp = 'biosimulators_smoldyn/Min1/model.txt'

sed_fp = 'biosimulators_smoldyn/Min1/simulation.sedml'

model_lang = 'urn:sedml:language:smoldyn'


parameters = get_parameters_variables_outputs_for_simulation(model_fp, model_lang)

variables = parameters[2]

validation = validate_variables(variables)

config = Config(REPORT_FORMATS=[ReportFormat.csv])

results, log = exec_sed_doc(doc=sed_fp, working_dir='biosimulators_smoldyn/Min1', base_out_path=os.getcwd(), config=config)

print(results)
