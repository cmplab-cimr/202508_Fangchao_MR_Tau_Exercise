from .fruitfly_compart_model_20250815 import *
from .fruitfly_compart_model_20250815 import (
    return_data_content as old_return_data_content,
    tissue_mid_name_list_dict_constructor as old_tissue_mid_name_list_dict_constructor,
    tissue_flux_name_list_dict_constructor as old_tissue_flux_name_list_dict_constructor
)


model_name = ModelList.complete_model_gut
user_defined_model = model_loader(model_name)
mfa_model_obj = common_model_constructor(user_defined_model)

def return_data_content(running_mode, average=False, **kwargs):
    return old_return_data_content(
        running_mode, average=average, no_sucrose_2=False, sucrose_no_exercise=True, **kwargs)


def tissue_mid_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    return old_tissue_mid_name_list_dict_constructor(current_model_name=current_model_name, **kwargs)


def tissue_flux_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    return old_tissue_flux_name_list_dict_constructor(current_model_name=current_model_name, **kwargs)
