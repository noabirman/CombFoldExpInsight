import os
import sys

from pdb_to_json import set_b_factor_to_100

#a test for b_factor function
if __name__ == '__main__':
    if len(sys.argv) == 2:
        pdb_path = os.path.abspath(sys.argv[1])
        set_b_factor_to_100(pdb_path)
    else:
        print("usage: <script> pdb_path, output_path")