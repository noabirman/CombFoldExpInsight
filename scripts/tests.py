import os
import sys
from cf_exi_insight import update_chain_list

from cf_exi_insight import set_b_factor_to_100, update_b_factors_preserve_metadata

#a test for b_factor function
if __name__ == '__main__':
    if len(sys.argv) == 2:
        combfold_exi_rep_path = os.path.abspath(sys.argv[1])
        chain_list_path = os.path.join(combfold_exi_rep_path, '_unified_representation', 'assembly_output',
                                       'chain.list')
        update_chain_list(chain_list_path)
    else:
        print("usage: <script> pdb_path, output_path")