#Use to collect, visualise atomic structures
from ase.io import read

import numpy as np #Mathematical calculations

#Check length consistency of files
def check_atom_consistency(item, ini_q, fin_q):

    # Bader vs Bader
    if len(ini_q) != len(fin_q): #i.e., Check that Bader charge lengths are consistent
        raise ValueError( #If not, raise an error
            f"Atom count mismatch in {item['transition']}: "
            f"ini={len(ini_q)}, fin={len(fin_q)}"
        )

    # POSCAR vs Bader
    ini_len = len(read(item["ini_structure"]))
    fin_len = len(read(item["fin_structure"]))

    if (fin_len != len(fin_q)) or (ini_len != len(ini_q)): #i.e., Check that Bader charge and POSCAR/CONTCAR lengths are consistent
        raise ValueError( #If not, raise an error
            f"Mismatch between POSCAR/CONTCAR and Bader atoms in {item['transition']}"
        )

#Check that Bader coordinates and POSCAR/CONTCAR coordinates are equivalent
def check_bader_alignment(atoms, acf_coords, tol=1e-3):

    #N.B. This tolerance is to not be confused with that of Bader charge difference.
    #This difference is to check that atoms in POSCAR/CONTCAR and ACF.dat are the same.
    #Here, difference should not be too large, and hence it is defaulted.

    pos = atoms.get_positions()  # Cartesian from ASE

    #Calculate difference based on absolute distance
    diff = np.linalg.norm(pos - acf_coords, axis=1)

    #Should be within threshold
    if np.max(diff) > tol:
        raise ValueError( #Else raise an error
            f"ACF.dat coordinates do not match POSCAR/CONTCAR "
            f"(max diff = {np.max(diff):.4f} Å)"
        )