import os
import shutil

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
            merged_output_folder = os.path.join(output_folder, complex_name, 'merged_af_output')
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

def process_complex_folders(af_folder, output_folder):
    """
    Process each complex folder to merge af_output folders and run cf_exi_insight.

    Args:
        af_folder (str): Path to the base folder containing complex folders.
        output_folder (str): Path to the output folder where merged af_output will be stored.
    """
    merge_af_output_folders(af_folder, output_folder)
    for complex_name in os.listdir(af_folder):
        complex_folder = os.path.join(af_folder, complex_name)
        if os.path.isdir(complex_folder):
            merged_output_folder = os.path.join(output_folder, complex_name, 'merged_af_output')
            # run cf_exi_insight on pdb file and merged af_output
            pdb_path = os.path.join(pdb_folder, complex_name, '.pdb')
            current_output_path = os.path.join(output_folder, complex_name)
            os.system(f"python {pdb_path} {current_output_path} {merged_output_folder}")

if __name__ == '__main__':
    pdb_folder = 'path/to/pdb_folders'
    af_folder = 'path/to/alpha_fold_folders'
    output_folder = 'path/to/output_location'
    merge_af_output_folders(af_folder, output_folder)