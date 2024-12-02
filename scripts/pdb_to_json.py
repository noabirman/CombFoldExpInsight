import os
import sys
from typing import Dict
from run_on_pdbs import get_chain_to_seq
from libs.utils_classes import SubunitInfo, save_subunits_info, SubunitsInfo  # Make sure SubunitInfo is available
from Bio.PDB import PDBParser, PDBIO
import re
from typing import List, Tuple


def get_first_non_x_index(sequence: str) -> int:
    """
    Get the index of the first non-'X' character in a sequence.

    Args:
        sequence (str): The sequence to search.

    Returns:
        int: The index of the first non-'X' character, or -1 if none is found.
    """
    for i, aa in enumerate(sequence):
        if aa != 'X':
            return i
    return -1


def create_subunit_info_from_chain_seq(chain_to_seq: Dict[str, str]) -> SubunitsInfo:
    """
    Create a SubunitsInfo object from a dictionary of chain sequences.

    Args:
        chain_to_seq (Dict[str, str]): A dictionary mapping chain IDs to sequences.

    Returns:
        SubunitsInfo: A dictionary mapping chain IDs to SubunitInfo objects.
    """
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


def find_long_x_sequences(sequence) -> List[Tuple[int, int]]:
    """
    Find all sequences of 'X' characters longer than 100 in a given text.

    Args:
        sequence (str): The text to search.

    Returns:
        List[Tuple[int, int]]: A list of tuples containing the start and end indices of each long 'X' sequence.
    """
    # Find all sequences of 'X' characters
    x_sequences = re.finditer(r'X+', sequence)
    # Filter for sequences with length greater than 100 and store their positions and lengths
    result = [(match.start(), match.end()) for match in x_sequences if
              (match.end() - match.start()) > 100]
    return result


def merge_exi_and_full_seq_into_subunit_info(experimental: Dict[str, str], full_seq: Dict[str, str]) -> Tuple[SubunitsInfo, SubunitsInfo]:
    """
    Merge experimental sequences with full sequences to create SubunitInfo objects.

    Args:
        experimental (Dict[str, str]): A dictionary mapping chain IDs to experimental sequences.
        full_seq (Dict[str, str]): A dictionary mapping chain IDs to full sequences.

    Returns:
        Dict[str, SubunitInfo]: A dictionary mapping subunit names to SubunitInfo objects.
    """
    merged_subunits_info = {}
    exi_subunits_info = {}
    for chain_id, sequence in experimental.items():
        missing_subunits = find_long_x_sequences(sequence) # Find long 'X' sequences
        if not missing_subunits:  # in case no holes, keep the original seq and chain_id
            subunit_name = chain_id + "EXI" # EXI for existing
            start_res_index = get_first_non_x_index(sequence) + 1 # base 1
            new_subunit = SubunitInfo(name=subunit_name, chain_names=[chain_id],
                                                      start_res=start_res_index,
                                                      sequence=sequence[start_res_index:])
            exi_subunits_info[subunit_name] = new_subunit
            merged_subunits_info[subunit_name] = new_subunit
        # if there are holes, we need to create new subunits
        subunit_num = 1
        start_index = 0
        for hole in missing_subunits:
            if hole[0] != 0:  # if there is a sequence before the hole
                subunit_name = chain_id + str(subunit_num) + "EXI"
                start_res_index = get_first_non_x_index(sequence[start_index:hole[0]]) + start_index + 1 # base 1
                new_subunit = SubunitInfo(name=subunit_name, chain_names=[chain_id],
                                                          start_res=start_res_index + 1,  # base 1
                                                          sequence=sequence[start_res_index:hole[0]])
                exi_subunits_info[subunit_name] = new_subunit
                merged_subunits_info[subunit_name] = new_subunit
                subunit_num += 1
            # filling hole:
            subunit_name = chain_id + str(subunit_num) + "MIS"
            merged_subunits_info[subunit_name] = SubunitInfo(name=subunit_name, chain_names=[chain_id],
                                                      start_res=hole[0] + 1,  # base 1
                                                      sequence=full_seq[chain_id][hole[0]:hole[1]])
            subunit_num += 1
            if hole == missing_subunits[-1]:  # if we are in the last hole
                if hole[1] < len(sequence):
                    subunit_name = chain_id + str(subunit_num) + "EXI"
                    new_subunit = SubunitInfo(name=subunit_name, chain_names=[chain_id],
                                                              start_res=hole[1] + 1,  # base 1
                                                              sequence=sequence[hole[1]:])
                    exi_subunits_info[subunit_name] = new_subunit
                    merged_subunits_info[subunit_name] = new_subunit
            start_index = hole[1]
    return merged_subunits_info, exi_subunits_info


def set_b_factor_to_100(pdb_path: str):
    """
    Set the B-factor of all atoms in a PDB file to 100.

    Args:
        pdb_path (str): The path to the PDB file.
    """
    # Parse the PDB file
    parser = PDBParser()
    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
    structure = parser.get_structure(pdb_id, pdb_path)

    # Modify the B-factor of all atoms
    for atom in structure.get_atoms():
        atom.bfactor = 100.0  # Set a new B-factor value

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_id + "_b_factor_100.pdb")


def update_b_factors_preserve_metadata(input_pdb: str, output_pdb: str, new_b_factor: float = 100.0) -> None:
    """
    Update the B-factors of atoms in a PDB file while preserving metadata.

    Args:
        input_pdb (str): The path to the input PDB file.
        output_pdb (str): The path to the output PDB file.
        new_b_factor (float, optional): The new B-factor value. Defaults to 100.0.
    """
    with open(input_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Update the B-factor column (columns 61-66 in the PDB format)
                updated_line = line[:60] + f"{new_b_factor:6.2f}" + line[66:]
                outfile.write(updated_line)
            else:
                # Write metadata and other lines unchanged
                outfile.write(line)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        # first part: take a pdb file and make subunits.json from the experimental data with XXXX where missing location of Amino Acids
        pdb_path, output_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
        # chain_to_seq = get_chain_to_seq(pdb_path, use_seqres=False)  # Corrected use_seqres
        # subunits_info = create_subunit_info_from_chain_seq(chain_to_seq)
        # save_subunits_info(subunits_info,
        #                    os.path.join(output_path, 'subunits.json'))  # Corrected sub_units -> subunits_info
        # Second Part: define new subunits where there is a sequence of more than 100 X's and replace the X's with the known sequence
        full_known_seq = get_chain_to_seq(pdb_path, use_seqres=True)  # now use_seqres
        # Remove the '6I3M:' prefix from chain IDs in full_known_seq
        cleaned_full_known_seq = {
            chain_id.split(":")[1]: seq for chain_id, seq in full_known_seq.items()
        }
        # exp_seq_after_fiiling = add_missing_subunits_seq(chain_to_seq, cleaned_full_known_seq)
        merged_subunits_info, exi_subunits_info = merge_exi_and_full_seq_into_subunit_info(chain_to_seq, cleaned_full_known_seq)

        # check if filling was good
        # print("Key-Value pairs after filling:", exp_seq_after_fiiling.items())

        # combined_subunits_info = create_subunit_info_from_chain_seq(exp_seq_after_fiiling)
        save_subunits_info(merged_subunits_info, os.path.join(output_path, 'merged_subunits.json'))
        save_subunits_info(exi_subunits_info, os.path.join(output_path, 'exi_subunits.json'))
    else:
        print("usage: <script> pdb_path, output_path")