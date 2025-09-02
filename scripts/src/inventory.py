from .common.packages import ValueEnum


class ExperimentName(ValueEnum):
    fruitfly_compart_model_20250731 = 'fruitfly_compart_model_20250731'
    fruitfly_compart_model_20250815 = 'fruitfly_compart_model_20250815'
    fruitfly_compart_model_20250815_gut = 'fruitfly_compart_model_20250815_gut'
    fruitfly_compart_infusion_model_20250815 = 'fruitfly_compart_infusion_model_20250815'
    fruitfly_compart_infusion_model_20250815_gut = 'fruitfly_compart_infusion_model_20250815_gut'

# default_data_model_name = ExperimentName.fruitfly_compart_model_20250731
# default_data_model_name = ExperimentName.fruitfly_compart_model_20250815
# default_data_model_name = ExperimentName.fruitfly_compart_model_20250815_gut
# default_data_model_name = ExperimentName.fruitfly_compart_infusion_model_20250815
default_data_model_name = ExperimentName.fruitfly_compart_infusion_model_20250815_gut


class MFARunningMode(ValueEnum):
    flux_analysis = 'flux_analysis'
    raw_experimental_data_plotting = 'raw_experimental_data_plotting'
    result_process = 'result_process'
    solver_output = 'solver_output'


data_model_comment = {
    ExperimentName.fruitfly_compart_model_20250731:
        'Data from fruit flies',
    ExperimentName.fruitfly_compart_model_20250815:
        'Data from fruit flies',
    ExperimentName.fruitfly_compart_model_20250815_gut:
        'Data from fruit flies',
    ExperimentName.fruitfly_compart_infusion_model_20250815:
        'Data from fruit flies',
    ExperimentName.fruitfly_compart_infusion_model_20250815_gut:
        'Data from fruit flies',
}
