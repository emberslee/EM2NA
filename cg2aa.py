import os
import tempfile
import subprocess
import numpy as np
from pdbio import read_atom23_na, nuc_type_1_to_atoms
from collections import OrderedDict

# 0702
# Fix bug in N1/N9 assignment
def write_residues_as_pdb(filename, res_pos, res_types, chain='A', ter=False):
    assert len(res_types) >= len(res_pos)
    f = open(filename, 'w')
    for i in range(len(res_pos)):
        if res_types[i] in ['A', 'G', 'DA', 'DG']:
            atom_name = " N9 "
        else:
            atom_name = " N1 "
       
        f.write("ATOM  {:5d} {:4s}   {:1s}{:>2s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(i+1, " P  "   , res_types[i], chain, i+1, res_pos[i, 0][0], res_pos[i, 0][1], res_pos[i, 0][2]))
        f.write("ATOM  {:5d} {:4s}   {:1s}{:>2s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(i+1, " C4'"   , res_types[i], chain, i+1, res_pos[i, 1][0], res_pos[i, 1][1], res_pos[i, 1][2]))
        f.write("ATOM  {:5d} {:4s}   {:1s}{:>2s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(i+1, atom_name, res_types[i], chain, i+1, res_pos[i, 2][0], res_pos[i, 2][1], res_pos[i, 2][2]))
        if ter:
            f.write("TER\n")

    if not ter:
        f.write("TER\n")

# 0702
# Optional, sort the pdb lines
def sort_lines_by_atom_name(fin, fout):
    with open(fin, 'r') as f:
        lines = f.readlines()

    res_dict = OrderedDict()
    for line in lines:
        line = line.strip("\n")
        if line.startswith("ATOM") and len(line) >= 54:
            chain_id = line[20:22].strip()
            res_id = line[22:26].strip()
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            # ignore H
            if atom_name[0] == 'H':
                continue
            label = chain_id + "_" + res_id + "_" + res_name
            if label not in res_dict:
                res_dict[label] = OrderedDict()

            res_dict[label][atom_name] = line

    new_lines = []
    for label, atom_dict in res_dict.items():
        res_name = label.split("_")[-1]
        
        for atom_name in nuc_type_1_to_atoms[res_name]:
            if atom_name is None:
                continue

            line = atom_dict[atom_name]
            new_lines.append(line)
    
    n_atom = 1
    with open(fout, 'w') as f:
        for line in new_lines:
            new_line = line[:6] + "{:>5d}".format(n_atom) + line[11:]
            f.write(new_line + '\n')
            n_atom += 1

def full_atom(atom3_pos, res_types=None, lib_dir="~/", flag="", natype='rna', temp_dir=None):
    if res_types is None:
        res_types = ['U'] * len(atom3_pos)
    assert len(atom3_pos) == len(res_types)

    # Which program to use
    if natype in ['dna', 'DNA']:
        program = 'cg2aa_dna'
    else:
        program = 'cg2aa'

    # Full atom
    with tempfile.TemporaryDirectory() as __temp_dir:
        if temp_dir is None:
            temp_dir = __temp_dir

        # Write atom3 to PDB
        fin = os.path.join(temp_dir, f"{flag}3p.pdb")
        write_residues_as_pdb(fin, atom3_pos, res_types, ter=False)

        # Rebuild all atom
        fout = os.path.join(temp_dir, f"{flag}aa.pdb")
        cmd = f"{lib_dir}/bin/{program} {fin} {fout} 4 "
        print("# Running command {}".format(cmd))
        result = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            print("# Cannot rebuild full-atom")
            print("# STDOUT")
            print(result.stdout)
            print("# STDERR")
            print(result.stderr)

        # Read PDB
        #foutx = os.path.join(temp_dir, f"{flag}aa_sorted.pdb")
        #sort_lines_by_atom_name(fout, foutx)
        #atom23_pos, atom23_mask, _, _, _ = read_atom23_na(foutx, ignore_hetatm=True)

        atom23_pos, atom23_mask, _, _, _ = read_atom23_na(fout, ignore_hetatm=True)

    return atom23_pos, atom23_mask


