from __future__ import annotations
from typing import List
import re
import os
import math
from copy import deepcopy
from collections import defaultdict, OrderedDict

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio import BiopythonWarning
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.mmcifio import MMCIFIO, mmcif_order
import warnings
warnings.filterwarnings('ignore')

def flatten(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result

def read_pdb_simple(filename, targets=None):
    with open(filename, 'r') as f:
        lines = f.readlines()

    coords = []
    for line in lines:
        if line.startswith("ATOM"):
            atom_type = line[12:16]
            if targets is not None:
                if atom_type not in targets:
                    continue

            coord = [float(line[i:i+8]) for i in [30, 38, 46]]
            coords.append(coord)
    return np.asarray(coords).astype(np.float32)



def read_pdb(filename, single=False):
    with open(filename, 'r') as f:
        lines = f.readlines()
    models = []
    chains = []
    chain = []
    for line in lines:
        if line.startswith("MODEL") and len(chains) > 0:
            models.append(chains)
            chains = []

        if line.startswith("TER"):
            if len(chain) > 0:
                chains.append(np.asarray(chain, dtype=np.float32))
                chain = []

        if line.startswith("ATOM"):
            xyz = [float(x) for x in [line[i:i+8] for i in [30, 38, 46]]]
            chain.append(np.asarray(xyz))

    if len(chains) > 0:
        models.append(chains)
        chains = []

    if single:
        models = flatten(models)

    return models



"""Nucleotides"""
nuc_type_1 = ["A", "G", "C", "U", "T", "I", "N", "X"]
nuc_type_3 = [
    "A", "G", "C", "U", "DT", "I", "N", "X", 
]
nuc_type_1_to_3 = {
    "A": "A",
    "G": "G",
    "C": "C",
    "U": "U",
    "T": "DT",
    "I": "I",
    "N": "N",
    "X": "X",
}
nuc_type_3_to_1 = {
    "A": "A",
    "G": "G",
    "C": "C",
    "U": "U",
    "DA": "A",
    "DG": "G",
    "DC": "C",
    "DT": "T",
    "I": "I",
    "N": "N",
    "X": "X",
}
nuc_type_1_to_index = {
    "A": 0,
    "G": 1,
    "C": 2,
    "U": 3,
    "T": 4,
    "I": 5,
    "N": 6,
    "X": 7,
}
nuc_type_index_to_1 = {
    0: "A",
    1: "G",
    2: "C",
    3: "U",
    4: "T",
    5: "I",
    6: "N",
    7: "X",
}
nuc_type_1_to_atoms = OrderedDict({
    "A": ["P"  , "OP1", "OP2",
          "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
          "N9" , "C8" , "N7" , "C5" , "C6" , "N6" , "N1" , "C2" , None , "N3" , "C4" ,
         ],

    "G": ["P"  , "OP1", "OP2",
          "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
          "N9" , "C8" , "N7" , "C5" , "C6" , "O6" , "N1" , "C2" , "N2" , "N3" , "C4" ,
         ],

    "C": ["P"  , "OP1", "OP2",
          "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
          "N1" , "C2" , "O2" , "N3" , "C4" , "N4" , "C5" , "C6" , None , None , None ,
         ],

    "U": ["P"  , "OP1", "OP2",
          "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
          "N1" , "C2" , "O2" , "N3" , "C4" , "O4" , "C5" , "C6" , None , None , None ,
         ],

    "T": ["P"  , "OP1", "OP2",
          "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", None , "C1'",
          "N1" , "C2" , "O2" , "N3" , "C4" , "O4" , "C5" , "C7" , "C6" , None , None ,
         ],

})


chain_names = []
# One-letter chain id
for i in range(26):
    chain_names.append(chr(ord('A') + i))
for i in range(26):
    chain_names.append(chr(ord('a') + i))
for i in range(10):
    chain_names.append(str(i))
# Two-letter chain id
for i in range(26):
    for k in range(26):
        chain_names.append( chr(ord('A')+i) + chr(ord('A')+k))
for i in range(26):
    for k in range(26):
        chain_names.append( chr(ord('a')+i) + chr(ord('a')+k))
# In total we have 26 * 2 + 10 + 2 * 26 * 26 ~= 1500 chain ids



# 0702
class CIFXIO(MMCIFIO):
    def _save_dict(self, out_file):
        label_seq_id = deepcopy(self.dic["_atom_site.auth_seq_id"])
        auth_seq_id = deepcopy(self.dic["_atom_site.auth_seq_id"])
        self.dic["_atom_site.label_seq_id"] = label_seq_id
        self.dic["_atom_site.auth_seq_id"] = auth_seq_id

        # Adding missing "pdbx_formal_charge", "auth_comp_id", "auth_atom_id" to complete a record
        N = len(self.dic["_atom_site.group_PDB"])
        self.dic["_atom_site.pdbx_formal_charge"] = ["?"]*N
        self.dic["_atom_site.auth_comp_id"] = deepcopy(self.dic["_atom_site.label_comp_id"])
        self.dic["_atom_site.auth_asym_id"] = deepcopy(self.dic["_atom_site.label_asym_id"])
        self.dic["_atom_site.auth_atom_id"] = deepcopy(self.dic["_atom_site.label_atom_id"])

        # Handle an extra space at the end of _atom_site.xxx
        _atom_site = mmcif_order["_atom_site"]
        _atom_site = [x.strip() + " " for x in _atom_site]
        mmcif_order["_atom_site"] = _atom_site

        new_dic = defaultdict()
        for k, v in self.dic.items():
            if k[:11] == "_atom_site.":
                new_k = k.strip() + " "
            else:
                new_k = k
            new_dic[new_k] = v
        self.dic = new_dic

        return super()._save_dict(out_file)


# Get atom23
def read_atom23_na(filename, ignore_hetatm=False, check_dna=False):
    # Read file
    if filename.endswith("pdb"):
        parser = PDBParser()
    elif filename.endswith("cif"):
        parser = MMCIFParser()
    else:
        raise Exception("Error only support pdb/cif file")
    structure = parser.get_structure('pdb', filename)

    # Only use model 0
    model = structure[0]

    # Extract residue information
    atom23_pos = []
    atom23_mask = []
    res_types = []
    res_indices = []
    chain_indices = []
    res_is_dna = []

    for i, chain in enumerate(model):
        prev_residue_number = None
        for residue_index, residue in enumerate(chain):
            residue_number = residue.get_id()[1]
            if prev_residue_number is None or residue_number != prev_residue_number:
                resname = residue.get_resname().strip()
                hetfield, resseq, icode = residue.get_id()

                # If hetatom
                if ignore_hetatm and hetfield != " ":
                    continue

                # 0128
                is_dna = False
                if not check_dna:
                    if len(resname) == 2 and resname[0] == "D":
                        resname = resname[1]
                    if resname == "T":
                        resname = "U"
                else:
                    if len(resname) == 2 and resname[0] == "D":
                        resname = resname[1]
                        is_dna = True

                # If resname is not standard residues
                if resname not in nuc_type_1:
                    resname = "U"

                # 0128
                res_is_dna.append(is_dna)

                res_types.append(nuc_type_1_to_index[resname])
                res_indices.append(residue_number)
                prev_residue_number = residue_number

            coords = []
            mask = []
            # Get atom types for current residue
            atom_names = nuc_type_1_to_atoms[resname]
            while len(atom_names) < 23:
                atom_names.append(None)

            for atom_index in atom_names:
                try:
                    atom = residue[atom_index]
                    coords.append(atom.get_coord())
                    mask.append(1)
                except KeyError:
                    coords.append([float("nan") for _ in range(3)])
                    mask.append(0)
            atom23_pos.append(coords)
            atom23_mask.append(mask)
            chain_indices.append(i)

    # Convert to NumPy arrays
    atom23_pos = np.array(atom23_pos).astype(np.float32)
    atom23_mask = np.array(atom23_mask).astype(np.int32)
    res_types = np.array(res_types).astype(np.int32)
    res_indices = np.array(res_indices).astype(np.int32)
    chain_indices = np.array(chain_indices).astype(np.int32)
    # 0129
    res_is_dna = np.asarray(res_is_dna).astype(np.int32)

    # 0128
    if not check_dna:
        return atom23_pos, atom23_mask, res_types, res_indices, chain_indices
    else:
        return atom23_pos, atom23_mask, res_types, res_indices, chain_indices, res_is_dna


def fix_quotes(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        # Handle '' issues
        fixed_lines = [re.sub(r"'([A-Z0-9]+)''", '"\\1\'"', line) for line in lines]
        temp_filename = filename + ".tmp"
        with open(temp_filename, 'w') as f:
            f.writelines(fixed_lines)
        os.replace(temp_filename, filename)
    except IOError as e:
        print(f"Error processing file: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
    

# chain atom23 to pdb
def chain_atom23_to_pdb(filename, atom23_pos, atom23_mask, res_types, chain='A', occupancy=None, bfactors=None, is_dna=False):
    assert len(atom23_pos) == len(atom23_mask)
    assert len(atom23_pos) == len(res_types)

    struct = StructureBuilder()
    struct.init_structure("1")
    struct.init_seg("1")
    struct.init_model("1")
    struct.init_chain(chain)

    if bfactors is None:
        bfactors = np.ones(len(atom23_pos)) * 100.
    if occupancy is None:
        occupancy = np.ones(len(atom23_pos))


    n_atom = 0
    for i in range(len(atom23_pos)):
        res_name_1 = res_types[i]
        res_name_1 = "C" if res_name_1 == "N" else res_name_1
        res_name_3 = nuc_type_1_to_3[res_name_1]
        bfactor = bfactors[i]
        atom_names = nuc_type_1_to_atoms[res_name_1]

        field_name = " " if res_name_1 != "N" else "H"

        # If is DNA
        if is_dna:
            if res_name_3 == 'U':
                res_name_3 = 'T'
            res_name_3 = "D" + res_name_3

        struct.init_residue(res_name_3, field_name, i+1, " ")
        for atom_name, pos, mask in zip(
            atom_names, atom23_pos[i], atom23_mask[i]
        ):
            if atom_name is None or \
               mask < 1 or \
               np.any(np.isnan( pos )):
                continue

            struct.set_line_counter(n_atom+1)
            struct.init_atom(
                name=atom_name,
                coord=pos,
                b_factor=bfactor,
                occupancy=1.0,
                altloc=" ",
                fullname=atom_name,
                element=atom_name[0],
            )
            n_atom += 1

    struct = struct.get_structure()
    if filename.endswith("cif"):
        io = CIFXIO()
        io.set_structure(struct)
        io.save(filename)
        fix_quotes(filename)
    else:
        io = PDBIO()
        io.set_structure(struct)
        io.save(filename, write_end=False)

    #io = PDBIO()
    #io.set_structure(struct)
    #io.save(filename, write_end=False)


# Write multiple chains to a pdb file
def chains_atom23_to_pdb(
        filename : str, 
        chains_atom23_pos : List[np.ndarray], 
        chains_atom23_mask : List[np.ndarray], 
        chains_res_types : list[np.ndarray], 
        chains_occupancy=None, 
        chains_bfactors=None, 
        chains_is_dna=None,
    ):
    # For different chains
    assert len(chains_atom23_pos) == len(chains_atom23_mask)
    assert len(chains_atom23_pos) == len(chains_res_types)
    if chains_is_dna is not None:
        assert len(chains_atom23_pos) == len(chains_is_dna)
    else:
        chains_is_dna = [False] * len(chains_atom23_pos)

    struct = StructureBuilder()
    struct.init_structure("1")
    struct.init_seg("1")
    struct.init_model("1")

    n_atom = 0
    for k in range(len(chains_atom23_pos)):
        # For each nuc
        struct.init_chain(chain_names[k])

        atom23_pos = chains_atom23_pos[k]
        atom23_mask = chains_atom23_mask[k]
        res_types = chains_res_types[k]

        for i in range(len(atom23_pos)):
            # For each atom
            res_name_1 = res_types[i]
            res_name_1 = "C" if res_name_1 == "N" else res_name_1
            res_name_3 = nuc_type_1_to_3[res_name_1]
            bfactor = 100.0
            atom_names = nuc_type_1_to_atoms[res_name_1]
    
            field_name = " " if res_name_1 != "N" else "H"
    
            # If is DNA
            if chains_is_dna[k]:
                if res_name_3 == 'U':
                    res_name_3 = 'T'
                res_name_3 = "D" + res_name_3
    
            struct.init_residue(res_name_3, field_name, i+1, " ")
            for atom_name, pos, mask in zip(
                atom_names, atom23_pos[i], atom23_mask[i]
            ):
                if atom_name is None or \
                   mask < 1 or \
                   np.any(np.isnan( pos )):
                    continue
    
                struct.set_line_counter(n_atom+1)
                struct.init_atom(
                    name=atom_name,
                    coord=pos,
                    b_factor=bfactor,
                    occupancy=1.0,
                    altloc=" ",
                    fullname=atom_name,
                    element=atom_name[0],
                )
                n_atom += 1

    struct = struct.get_structure()
    if filename.endswith("cif"):
        io = CIFXIO()
        io.set_structure(struct)
        io.save(filename)
        fix_quotes(filename)
    else:
        io = PDBIO()
        io.set_structure(struct)
        io.save(filename, write_end=False)


    #io = PDBIO()
    #io.set_structure(struct)
    #io.save(filename, write_end=False)



# Only for eval_final.py
# Write multiple chains to a pdb file
def res_chains_atom23_to_pdb(
        filename : str, 
        chains_atom23_pos : List[np.ndarray], 
        chains_atom23_mask : List[np.ndarray], 
        chains_res_types : list[np.ndarray], 
        chains_occupancy=None, 
        chains_bfactors=None, 
        chains_is_dna=None,
    ):
    # For different chains
    assert len(chains_atom23_pos) == len(chains_atom23_mask)
    assert len(chains_atom23_pos) == len(chains_res_types)
    if chains_is_dna is not None:
        assert len(chains_atom23_pos) == len(chains_is_dna)
    else:
        chains_is_dna = []
        for i in range(len(chains_atom23_pos)):
            chains_is_dna.append([False] * len(chains_atom23_pos[i]))

    struct = StructureBuilder()
    struct.init_structure("1")
    struct.init_seg("1")
    struct.init_model("1")

    n_atom = 0
    for k in range(len(chains_atom23_pos)):
        # For each nuc
        struct.init_chain(chain_names[k])

        atom23_pos = chains_atom23_pos[k]
        atom23_mask = chains_atom23_mask[k]
        res_types = chains_res_types[k]

        for i in range(len(atom23_pos)):
            # For each atom
            res_name_1 = res_types[i]
            res_name_1 = "C" if res_name_1 == "N" else res_name_1
            res_name_3 = nuc_type_1_to_3[res_name_1]
            bfactor = chains_bfactors[k][i]
            atom_names = nuc_type_1_to_atoms[res_name_1]
    
            field_name = " " if res_name_1 != "N" else "H"
    
            # If is DNA
            if chains_is_dna[k][i]:
                if len(res_name_3) == 1:
                    if res_name_3 == 'U':
                        res_name_3 = 'T'
                    res_name_3 = "D" + res_name_3
    
            struct.init_residue(res_name_3, field_name, i+1, " ")
            for atom_name, pos, mask in zip(
                atom_names, atom23_pos[i], atom23_mask[i]
            ):
                if atom_name is None or \
                   mask < 1 or \
                   np.any(np.isnan( pos )):
                    continue
    
                struct.set_line_counter(n_atom+1)
                struct.init_atom(
                    name=atom_name,
                    coord=pos,
                    b_factor=bfactor,
                    occupancy=1.0,
                    altloc=" ",
                    fullname=atom_name,
                    element=atom_name[0],
                )
                n_atom += 1

    struct = struct.get_structure()
    if filename.endswith("cif"):
        io = CIFXIO()
        io.set_structure(struct)
        io.save(filename)
        fix_quotes(filename)
    else:
        io = PDBIO()
        io.set_structure(struct)
        io.save(filename, write_end=False)

    #io = PDBIO()
    #io.set_structure(struct)
    #io.save(filename, write_end=False)
