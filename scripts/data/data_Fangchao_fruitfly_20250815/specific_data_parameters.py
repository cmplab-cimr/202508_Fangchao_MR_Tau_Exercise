from scripts.src.common.config import Direct, Keywords as CommonKeywords
from scripts.src.core.common.classes import TransformDict
from scripts.src.core.data.data_class import MIDData, MFAData, InputMetaboliteData
from ..common_input_metabolites import input_mid_data_processor
from ..config import np, pd, CoreConstants, ModelKeyword, CompleteDataset, standardize_metabolite_name
from ..data_Fangchao_fruitfly_20250731.specific_data_parameters import (
    SpecificParameters as OldSpecificParameters, Keyword as OldKeyword, glucose_6_labeled_input_metabolite_dict,
    data_metabolite_to_standard_name_dict as old_data_metabolite_to_standard_name_dict)


data_metabolite_to_standard_name_dict = TransformDict({
    **old_data_metabolite_to_standard_name_dict,
    'c5h9no3(5)': 'l-hydroxyproline',
    'c6h14o12p2(7)': 'fructose 1,6-bisphosphate',
    'fructose-1,6-bisphosphate': 'fructose 1,6-bisphosphate',
    'c20h32n6o12s2': 'l-palmitoylcarnitine',
})

class Keyword(OldKeyword):
    sucrose_2 = 'sucrose_2'


class SpecificParameters(OldSpecificParameters):
    data_name = 'data_Fangchao_fruitfly_20250815'

    def __init__(self):
        super(SpecificParameters, self).__init__()
        self.mixed_compartment_list = ('c', 'm')
        self.special_compartment_tissue_dict = {}

        self.current_direct = f'{Direct.data_direct}/{self.data_name}'
        self.pickle_file_path = f'{self.current_direct}/complete_data.pickle'

        self.input_metabolite_data_dict = input_mid_data_processor(glucose_6_labeled_input_metabolite_dict)
        self.combined_data = False
        self.deconvolve_other_atoms = False

        first_sheet_name_list = ['WB - + full', 'gut- full', 'gut (ctrl, MR-Tau) - + full']
        file_name_sheet_name_dict = {
            '13C6-Sucrose WB Gut full.xlsx': first_sheet_name_list,
        }
        self.first_sheet_name_list = first_sheet_name_list
        self.file_name_list = list(file_name_sheet_name_dict.keys())
        self.file_path_list = [f'{self.current_direct}/{file_name}' for file_name in file_name_sheet_name_dict.keys()]
        index_col_name = 'Name'
        gut_mfa_excluded_metabolite_dict = {
            'fructose 1,6-bisphosphate': None,
            'glucose 6-phosphate': None,
            'fructose 6-phosphate': None,
            'a-ketoglutarate': None,
            '2-phosphoglycerate/3-phosphoglycerate': None,
        }
        self.mfa_exclude_metabolites_dict = {
            Keyword.gut: gut_mfa_excluded_metabolite_dict,
        }
        self.condition_standard_name_dict = {
            'Ctrl': Keyword.ctrl,
            'MR': Keyword.mr,
            'Tau': Keyword.tau,
            'MR-Tau': Keyword.mr_tau,
            'MR+Tau': Keyword.mr_tau,
        }
        self.sheet_name_tissue_dict = {
            'WB - + full': Keyword.whole_body,
            'gut- full': Keyword.gut,
            'gut (ctrl, MR-Tau) - + full': Keyword.gut,
        }
        self.excluded_column_name_dict = {
            'ID': None,
            'MZ': None,
        }
        complete_data_parameter_dict_dict = {}
        for (file_name, each_file_sheet_name_list), file_path in zip(
                file_name_sheet_name_dict.items(), self.file_path_list):
            for sheet_index, each_sheet_name in enumerate(each_file_sheet_name_list):
                if sheet_index == 0:
                    skip_row_num = 0
                else:
                    skip_row_num = 1
                complete_data_parameter_dict_dict[f'{file_name}_{each_sheet_name}'] = {
                    'xlsx_file_path': file_path,
                    'xlsx_sheet_name': each_sheet_name,
                    'skip_row_num': skip_row_num,
                    'index_col_name': index_col_name,
                    'mixed_compartment_list': self.mixed_compartment_list,
                    'to_standard_name_dict': data_metabolite_to_standard_name_dict,
                    'specific_dimension_metabolite_dict': self.specific_dimension_metabolite_dict,
                    'compound_replace_dict': self.compound_replace_dict,
                    'excluded_column_name_dict': self.excluded_column_name_dict,
                    'decipher_column_name': self._decipher_column_name,
                    'mfa_exclude_metabolites_dict': self.mfa_exclude_metabolites_dict,
                    'special_compartment_tissue_dict': self.special_compartment_tissue_dict,
                }
        self._complete_data_parameter_dict_dict = complete_data_parameter_dict_dict

    def add_data_sheet(self, sheet_name, current_data_dict):
        if sheet_name.endswith(self.first_sheet_name_list[-1]):
            labeling_name = Keyword.sucrose_2
        else:
            labeling_name = Keyword.sucrose
        if labeling_name not in self.complete_dataset:
            self.complete_dataset[labeling_name] = {}
        for exercise, each_exercise_dict in current_data_dict.items():
            if exercise not in self.complete_dataset[labeling_name]:
                self.complete_dataset[labeling_name][exercise] = {}
            for condition_name, each_condition_data in each_exercise_dict.items():
                if condition_name not in self.complete_dataset[labeling_name][exercise]:
                    self.complete_dataset[labeling_name][exercise][condition_name] = {}
                for index_num, each_index_data in each_condition_data.items():
                    if index_num not in self.complete_dataset[labeling_name][exercise][condition_name]:
                        self.complete_dataset[labeling_name][exercise][condition_name][index_num] = each_index_data
                    else:
                        self.complete_dataset[labeling_name][exercise][condition_name][index_num].update(each_index_data)

    def return_multi_tissue_mfa_data_dict(self):
        current_input_metabolite_data_dict = input_mid_data_processor(glucose_6_labeled_input_metabolite_dict)
        final_mfa_data_dict = {}
        final_result_information_dict = {}
        for labeling, each_labeling_data_dict in self.complete_dataset.items():
            for exercise, each_exercise_dict in each_labeling_data_dict.items():
                for condition, each_condition_data_dict in each_exercise_dict.items():
                    for index_num, current_mid_data_dict in each_condition_data_dict.items():
                        project_name = self.project_name_generator(labeling, exercise, condition, index_num)
                        mfa_data_obj = MFAData(
                            project_name, current_mid_data_dict, current_input_metabolite_data_dict,
                            combined_data=self.combined_data, ratio_dict_to_objective_func=None)
                        final_mfa_data_dict[project_name] = mfa_data_obj
                        final_result_information_dict[project_name] = {
                            Keyword.labeling: labeling,
                            Keyword.exercise: exercise,
                            Keyword.condition: condition,
                            Keyword.index: index_num,
                        }
        return final_mfa_data_dict, final_result_information_dict

    def _decipher_column_name(self, raw_column_name):
        if raw_column_name.endswith('+'):
            exercise = Keyword.with_exec
        elif raw_column_name.endswith('-'):
            exercise = Keyword.no_exec
        else:
            raise ValueError()
        condition_str, index_str = raw_column_name[:-1].split(' ')
        condition_name = self.condition_standard_name_dict[condition_str]
        assert index_str.isdigit()
        index_num = int(index_str)
        return exercise, condition_name, index_num

    @property
    def _color_dict(self):
        from scripts.figures.common.config import ColorConfig
        color_dict_without_index = {
            (Keyword.no_exec, Keyword.ctrl): ColorConfig.normal_blue,
            (Keyword.no_exec, Keyword.mr): ColorConfig.orange,
            (Keyword.no_exec, Keyword.tau): ColorConfig.light_green,
            (Keyword.no_exec, Keyword.mr_tau): ColorConfig.light_bright_mixed_mid_red,
            (Keyword.with_exec, Keyword.ctrl): ColorConfig.dark_blue,
            (Keyword.with_exec, Keyword.mr): ColorConfig.dark_orange,
            (Keyword.with_exec, Keyword.tau): ColorConfig.green,
            (Keyword.with_exec, Keyword.mr_tau): ColorConfig.mid_red,
        }
        common_legend_color_dict = {
            f'{condition}{'+' if exercise == Keyword.with_exec else '-'}': color
            for (exercise, condition), color in color_dict_without_index.items()
        }
        legend_color_dict = {
            Keyword.sucrose: common_legend_color_dict,
            Keyword.sucrose_2: common_legend_color_dict,
        }
        return color_dict_without_index, legend_color_dict

    def result_information_dict_analysis(self, complete_result_information_dict):
        keyword = Keyword
        data_name = SpecificParameters.data_name

        color_dict_without_index, _ = self._color_dict
        color_order_dict = {}
        total_index = 0
        each_group_replicate_num = 3
        for (exercise_name, condition_name), color in color_dict_without_index.items():
            for current_index in range(1, 1 + each_group_replicate_num):
                color_order_dict[(exercise_name, condition_name, current_index)] = (total_index + current_index, color)
            total_index += each_group_replicate_num

        return_parameter_list = []
        for result_label, current_result_information_dict in complete_result_information_dict.items():
            major_key = current_result_information_dict[keyword.labeling]
            exercise = current_result_information_dict[keyword.exercise]
            exer_str = '+' if exercise == keyword.with_exec else '-'
            condition = current_result_information_dict[keyword.condition]
            index = current_result_information_dict[keyword.index]
            minor_key_list = []
            minor_key_str = f'{condition}_{index}{exer_str}'
            order_index, current_color = color_order_dict[(exercise, condition, index)]
            return_parameter_list.append(
                (result_label, major_key, minor_key_list, minor_key_str, current_color, order_index))
        return return_parameter_list

    def tissue_mid_name_list_dict_constructor(self, **kwargs):
        mid_name_list_in_one_tissue = [
            [
                'glucose', 'fructose 6-phosphate/glucose 6-phosphate', '2-phosphoglycerate/3-phosphoglycerate',
                'pyruvate', 'lactate', 'glycerol 3-phosphate', 'alanine', ],
            [
                'citrate/isocitrate', 'glutamate', 'glutamine', 'succinate',
                'malate', 'aspartate', ],
        ]
        complete_tissue_mid_name_list_dict = {
            Keyword.whole_body: mid_name_list_in_one_tissue,
            Keyword.gut: mid_name_list_in_one_tissue,
        }
        tissue_mid_name_list_dict = complete_tissue_mid_name_list_dict
        return (
            tissue_mid_name_list_dict, self.specific_dimension_metabolite_dict,
            Keyword.display_mid_name_dict, Keyword.display_tissue_name_dict)

    def experimental_figure_config_dict_generator(self):
        _, legend_color_dict = self._color_dict
        return {
            'figure_size': (8, 18),
            # 'size': [0.8, 1.1],
            # 'scale': 0.7,
            'legend': True,
            'legend_color_dict': legend_color_dict,
        }
