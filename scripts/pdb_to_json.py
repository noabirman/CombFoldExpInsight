import os
import sys
from typing import Dict
import libs.utils_classes
from run_on_pdbs import get_chain_to_seq
from libs.utils_classes import SubunitInfo, save_subunits_info  # Make sure SubunitInfo is available


def create_subunit_info_from_chain_seq(chain_to_seq: Dict[str, str]) -> Dict[str, SubunitInfo]:
    subunits_info = {}

    for chain_id, sequence in chain_to_seq.items():
        # Find the index of the first non-'X' residue
        try:
            start_res = sequence.index(next(filter(lambda x: x != 'X', sequence))) + 1
        except StopIteration:
            # If the sequence is all 'X', we skip this chain
            print(f"Skipping chain {chain_id} with no valid residues")
            continue

        # Create SubunitInfo object for this chain
        subunit_info = SubunitInfo(
            name=chain_id,  # Use the chain_id as the name
            chain_names=[chain_id],  # Use chain_id for chain_names
            start_res=start_res,  # Start residue index
            sequence=sequence  # Sequence from the original dict
        )

        # Add to the dictionary with chain_id as the key
        subunits_info[chain_id] = subunit_info

    return subunits_info


import re
from typing import List, Tuple


def find_long_x_sequences(text) -> List[Tuple[int, int, int]]:
    # Find all sequences of 'X' characters
    x_sequences = re.finditer(r'X+', text)
    # Filter for sequences with length greater than 100 and store their positions and lengths
    result = [(match.start(), match.end(), match.end() - match.start()) for match in x_sequences if
              (match.end() - match.start()) > 100]
    return result


def merge_exi_and_full_seq_into_subunit_info(experimental: Dict[str, str], full_seq: Dict[str, str]) -> Dict[
    str, SubunitInfo]:
    # combined_seq = {}
    subunits_info = {}
    for chain_id, sequence in experimental.items():
        missing_subunits = find_long_x_sequences(sequence)
        if not missing_subunits:  # in case no holes, keep the original seq and chain_id
            # combined_seq[chain_id + "EXI"] = sequence
            subunit_info = SubunitInfo(name=chain_id, chain_names=[chain_id + "EXI"], start_res=1, sequence=sequence)
        chain_num = 1
        start_index = 0
        for hole in missing_subunits:
            # before hole, if there is sequence:
            if hole[0] != 0:
                # combined_seq[chain_id + str(chain_num) + "EXI"] = sequence[start_index:hole[0]]
                subunit_name = chain_id + str(chain_num) + "EXI"
                subunits_info[subunit_name] = SubunitInfo(name=chain_id, chain_names=[subunit_name],
                                                          start_res=start_index + 1,  # base 1?
                                                          sequence=sequence[start_index:hole[0]])
                chain_num += 1
            # filling hole:
            # combined_seq[chain_id + str(chain_num) + "MIS"] = full_seq[chain_id][hole[0]:hole[2]]
            subunit_name = chain_id + str(chain_num) + "MIS"
            subunits_info[subunit_name] = SubunitInfo(name=chain_id, chain_names=[subunit_name],
                                                      start_res=hole[0] + 1,  # base 1?
                                                      sequence=full_seq[chain_id][hole[0]:hole[2]])
            chain_num += 1
            if hole == missing_subunits[-1]:  # if we are in the last hole
                if hole[2] < len(sequence):
                    # combined_seq[chain_id + str(chain_num) + "EXI"] = sequence[hole[2]:]
                    subunit_name = chain_id + str(chain_num) + "EXI"
                    subunits_info[subunit_name] = SubunitInfo(name=chain_id, chain_names=[subunit_name],
                                                              start_res=hole[2] + 1,  # base 1?
                                                              sequence=sequence[hole[2]:])
    return subunits_info


if __name__ == '__main__':
    if len(sys.argv) == 3:
        # first part: take a pdb file and make subunits.json from the experimental data with XXXX where missing location of Amino Acids
        pdb_path, output_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
        chain_to_seq = get_chain_to_seq(pdb_path, use_seqres=False)  # Corrected use_seqres
        subunits_info = create_subunit_info_from_chain_seq(chain_to_seq)
        save_subunits_info(subunits_info,
                           os.path.join(output_path, 'subunits.json'))  # Corrected sub_units -> subunits_info
        # Second Part: define new subunits where there is a sequence of more than 100 X's and replace the X's with the known sequence
        full_known_seq = get_chain_to_seq(pdb_path, use_seqres=True)  # now use_seqres
        # Remove the '6I3M:' prefix from chain IDs in full_known_seq
        cleaned_full_known_seq = {
            chain_id.split(":")[1]: seq for chain_id, seq in full_known_seq.items()
        }
        #exp_seq_after_fiiling = add_missing_subunits_seq(chain_to_seq, cleaned_full_known_seq)
        merged_subunits_info = merge_exi_and_full_seq_into_subunit_info(chain_to_seq, cleaned_full_known_seq)

            # check if filling was good
        #print("Key-Value pairs after filling:", exp_seq_after_fiiling.items())

        #combined_subunits_info = create_subunit_info_from_chain_seq(exp_seq_after_fiiling)
        save_subunits_info(merged_subunits_info, os.path.join(output_path, 'merged_subunits.json'))
    else:
        print("usage: <script> pdb_path, output_path")
