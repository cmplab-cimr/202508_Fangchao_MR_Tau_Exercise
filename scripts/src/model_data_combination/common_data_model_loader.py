from ..inventory import ExperimentName


def common_data_model_function_loader(model_name):
    # if model_name == ExperimentName.fruitfly_compart_model_20250731:
    #     from . import fruitfly_compart_model_20250731 as data_model_object
    if model_name == ExperimentName.fruitfly_compart_model_20250815:
        from . import fruitfly_compart_model_20250815 as data_model_object
    elif model_name == ExperimentName.fruitfly_compart_model_20250815_gut:
        from . import fruitfly_compart_model_20250815_gut as data_model_object
    elif model_name == ExperimentName.fruitfly_compart_infusion_model_20250815:
        from . import fruitfly_compart_infusion_model_20250815 as data_model_object
    elif model_name == ExperimentName.fruitfly_compart_infusion_model_20250815_gut:
        from . import fruitfly_compart_infusion_model_20250815_gut as data_model_object
    elif model_name == 'test_model_data':
        from . import test_model_data as data_model_object
    else:
        raise ValueError()
    return data_model_object
