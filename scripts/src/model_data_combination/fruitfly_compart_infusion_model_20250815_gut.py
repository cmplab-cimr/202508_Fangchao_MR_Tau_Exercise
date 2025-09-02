from .fruitfly_compart_model_20250815_gut import *
from .fruitfly_compart_infusion_model_20250815 import UserDefinedMFASettings as OldUserDefinedMFASettings


model_name = ModelList.complete_model_gut_infusion
user_defined_model = model_loader(model_name)
mfa_model_obj = common_model_constructor(user_defined_model)


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    optimization_num_for_each_suffix_dict = {
        None: 10000,
    }
user_defined_mfa_settings = UserDefinedMFASettings()
