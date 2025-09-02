from .fruitfly_compart_model_20250815 import *
from .fruitfly_compart_model_20250815 import UserDefinedMFASettings as OldUserDefinedMFASettings


model_name = ModelList.complete_model_whole_body_infusion
user_defined_model = model_loader(model_name)
mfa_model_obj = common_model_constructor(user_defined_model)


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    user_defined_specific_flux_range_dict = {
        **OldUserDefinedMFASettings.user_defined_specific_flux_range_dict,
        'PYR_supplement_net': ((keyword.whole_body, keyword.gut), (-200, 200)),
        'CIT_supplement_net': ((keyword.whole_body, keyword.gut), (-200, 200)),
    }
    optimization_num_for_each_suffix_dict = {
        None: 11000,
    }
user_defined_mfa_settings = UserDefinedMFASettings()
