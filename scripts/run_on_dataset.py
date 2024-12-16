import os
import shutil
from cf_exi_insight import run_cf_exi_insight

def merge_af_output_folders(af_folder, output_folder):
    """
    Merge af_output folders from af_iter_0 and af_iter_1 into one folder for each complex.

    Args:
        af_folder (str): Path to the base folder containing complex folders.
        output_folder (str): Path to the output folder where merged af_output will be stored.
    """
    for complex_name in os.listdir(af_folder):
        complex_folder = os.path.join(af_folder, complex_name)
        if os.path.isdir(complex_folder):
            merged_output_folder = os.path.join(output_folder, complex_name)
            af_iter_0_path = os.path.join(complex_folder, 'af_iter_0', 'af_output')
            af_iter_1_path = os.path.join(complex_folder, 'af_iter_1', 'af_output')

            if not os.path.exists(merged_output_folder):
                os.makedirs(merged_output_folder)

            for af_iter_path in [af_iter_0_path, af_iter_1_path]:
                if os.path.exists(af_iter_path):
                    for file_name in os.listdir(af_iter_path):
                        full_file_name = os.path.join(af_iter_path, file_name)
                        if os.path.isfile(full_file_name):
                            shutil.copy(full_file_name, merged_output_folder)

def process_complex_folders(pdb_folder, af_folder, output_folder):
    """
    Process each complex folder to merge af_output folders and run cf_exi_insight.

    Args:
        af_folder (str): Path to the base folder containing complex folders.
        output_folder (str): Path to the output folder where merged af_output will be stored.
    """
    for complex_file in os.listdir(pdb_folder):
        if complex_file.endswith('.pdb'):
            complex_name = os.path.splitext(complex_file)[0]
            af_current_folder = os.path.join(af_folder, complex_name)
            if os.path.isdir(af_current_folder):
                # run cf_exi_insight on pdb file and merged af_output
                pdb_path = os.path.join(pdb_folder, complex_file)
                current_output_path = os.path.join(output_folder, complex_name)
                if not os.path.exists(current_output_path):
                    os.makedirs(current_output_path)
                run_cf_exi_insight(pdb_path, af_current_folder, current_output_path)
if __name__ == '__main__':
    pdb_folder = '/cs/usr/bshor/sci/projects/af_combdock/runs/20221121_my_dataset_full/input_complexes'
    af_folder = '/cs/usr/bshor/sci/projects/af_combdock/runs/20221121_my_dataset_full/output'
    af_output_folder = '/cs/usr/tsori/cf_exinsight/af_ouputs'
    merge_af_output_folders(af_folder, af_output_folder)
    process_complex_folders(pdb_folder, af_output_folder, '/cs/usr/tsori/cf_exinsight/outputs')