#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#   Copyright 2026 Katarzyna Szleper-Shendy <katarzyna.szleper@polsl.pl> 
# 

import os
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial.distance import cdist



EXPOSED_SASA_THRESHOLD = 10
MIN_PATHWAY_LENGTH = 3
PROTON_ACCEPTOR_DISTANCE = 4.0  # Distance in Angstroms (from phenolic O to acceptor)
TYR_CU_DISTANCE = 3.5  # Distance in Angstroms (from tyrosine OH to copper)

# === Residue definitions ===
IMPORTANT_ATOMS = {
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "NE1", "CD1", "CG"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    # "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    # "CYS": ["SG"],
    # "MET": ["SD"]
}

@dataclass
class Residue:
    chain: str = field(hash=True)
    index: int = field(hash=True)
    name: str = field(hash=True)
    sasa: float
    coords: np.ndarray
    atoms: list
    def __repr__(self):
        return "{} {:>3d} {} \t{:>8.3f} {:>8.3f} {:>8.3f}".format(self.chain, self.index,
                self.name, self.coords[0], self.coords[1], self.coords[2])

def compute_relative_sasa(structure):
    """Relative solvent accessibility (%) per residue: the SASA of the whole residue over the
    residue type's maximum (Tien et al. 2013), as pipeline step 5 defines it. The pipeline
    passes step 5's values instead of calling this (hole_hopping_ehpath._run); until 2026-10-01
    the side-chain SASA was divided by these whole-residue maxima, which is no relative
    accessibility of anything."""
    MAX_SASA = {
    "A": 129,  # Alanine
    "R": 274,  # Arginine
    "N": 195,  # Asparagine
    "D": 193,  # Aspartate
    "C": 167,  # Cysteine
    "Q": 225,  # Glutamine
    "E": 223,  # Glutamate
    "G": 104,  # Glycine
    "H": 224,  # Histidine
    "I": 197,  # Isoleucine
    "L": 201,  # Leucine
    "K": 236,  # Lysine
    "M": 224,  # Methionine
    "F": 240,  # Phenylalanine
    "P": 159,  # Proline
    "S": 155,  # Serine
    "T": 172,  # Threonine
    "W": 285,  # Tryptophan
    "Y": 263,  # Tyrosine
    "V": 174   # Valine
    }
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    sr = ShrakeRupley()
    sr.compute(structure, level="A")
    data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                rescode = three_to_one.get(resname)
                if not rescode or rescode not in MAX_SASA:
                    continue
                resid = residue.id[1]
                residue_sasa = sum(atom.sasa for atom in residue)
                max_sasa = MAX_SASA[rescode]
                rel_sasa = (residue_sasa / max_sasa) * 100 if max_sasa else 0
                data.append({
                    "resid": resid,
                    "resname": resname,
                    "SASA": round(residue_sasa, 2),
                    "Relative SASA": round(rel_sasa, 2)
                })
    return pd.DataFrame(data)

def load_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    return structure

METAL_NAMES = {"CU", "ZN", "NI", "MN", "FE", "CO", "MO", "MG"}

def find_copper(model):
    """Return (residue, coord) for any Cu/Zn/Ni/Mg/… ion; priority Cu > Zn > others."""
    found = {}
    for chain in model:
        for residue in chain:
            rname = residue.get_resname().strip().upper().rstrip("0123456789")
            if rname in METAL_NAMES:
                for atom in residue.get_atoms():
                    aname = atom.get_name().strip().upper().rstrip("0123456789")
                    if aname in METAL_NAMES:
                        found.setdefault(rname, (residue, atom.get_coord()))
                        break
    for metal in ("CU", "ZN", "NI", "CO", "MN", "FE", "MO", "MG"):
        if metal in found:
            return found[metal]
    return None, None

def get_centroid(residue, atom_names):
    coords = [residue[atom].get_coord() for atom in atom_names if atom in residue]
    return np.mean(coords, axis=0) if len(coords) == len(atom_names) else None

def find_histidine_brace(model, cu_coord):
    his_residues = []
    for chain in model:
        for residue in chain:
            if residue.get_resname() == "HIS":
                centroid = get_centroid(residue, IMPORTANT_ATOMS["HIS"])
                if centroid is not None:
                    dist = np.linalg.norm(centroid - cu_coord)
                    his_residues.append((dist, (chain.id, residue.id, "HIS"), centroid))
    his_residues.sort()
    return his_residues[:2]

def has_proton_acceptor_near_tyr_oh(tyr_residue, model, his_brace):
    """
    Check if a tyrosine residue has a proton acceptor within 4Å of its phenolic O atom.
    Proton acceptors: carboxylate O from Asp/Glu, imidazole N from His, backbone O from carbonyl group.
    
    Args:
        tyr_residue: BioPython residue object for tyrosine
        model: BioPython model object
        his_brace: List of histidine brace residues to exclude
        
    Returns:
        bool: True if proton acceptor found within 4Å, False otherwise
    """
    # Get the OH atom from tyrosine (phenolic oxygen)
    if "OH" not in tyr_residue:
        return False
    
    oh_atom = tyr_residue["OH"]
    oh_coord = oh_atom.get_coord()
    
    # Get list of brace residue IDs to exclude
    brace_ids = [r[1] for r in his_brace]
    
    # Check all residues in the model for proton acceptors
    for chain in model:
        for residue in chain:
            resname = residue.get_resname()
            full_id = (chain.id, residue.id, resname)
            
            # Skip the tyrosine itself and histidine brace
            if full_id == (tyr_residue.get_parent().id, tyr_residue.id, "TYR"):
                continue
            if full_id in brace_ids:
                continue
            
            # Check Asp carboxylate O atoms (OD1, OD2)
            if resname == "ASP":
                for atom_name in ["OD1", "OD2"]:
                    if atom_name in residue:
                        acceptor_coord = residue[atom_name].get_coord()
                        dist = np.linalg.norm(oh_coord - acceptor_coord)
                        if dist <= PROTON_ACCEPTOR_DISTANCE:
                            return True
            
            # Check Glu carboxylate O atoms (OE1, OE2)
            elif resname == "GLU":
                for atom_name in ["OE1", "OE2"]:
                    if atom_name in residue:
                        acceptor_coord = residue[atom_name].get_coord()
                        dist = np.linalg.norm(oh_coord - acceptor_coord)
                        if dist <= PROTON_ACCEPTOR_DISTANCE:
                            return True
            
            # Check His imidazole N atoms (ND1, NE2)
            elif resname == "HIS":
                for atom_name in ["ND1", "NE2"]:
                    if atom_name in residue:
                        acceptor_coord = residue[atom_name].get_coord()
                        dist = np.linalg.norm(oh_coord - acceptor_coord)
                        if dist <= PROTON_ACCEPTOR_DISTANCE:
                            return True
            
            # Check backbone O from carbonyl group (present in all residues)
            if "O" in residue:
                acceptor_coord = residue["O"].get_coord()
                dist = np.linalg.norm(oh_coord - acceptor_coord)
                if dist <= PROTON_ACCEPTOR_DISTANCE:
                    return True
    
    return False

def is_tyr_oh_near_copper(tyr_residue, cu_coord):
    """
    Check if a tyrosine residue's phenolic OH is within 3.5Å of copper.
    
    Args:
        tyr_residue: BioPython residue object for tyrosine
        cu_coord: numpy array with copper coordinates
        
    Returns:
        bool: True if OH is within 3.5Å of copper, False otherwise
    """
    if cu_coord is None:
        return False
    
    # Get the OH atom from tyrosine (phenolic oxygen)
    if "OH" not in tyr_residue:
        return False
    
    oh_atom = tyr_residue["OH"]
    oh_coord = oh_atom.get_coord()
    
    # Calculate distance to copper
    dist = np.linalg.norm(oh_coord - cu_coord)
    return dist <= TYR_CU_DISTANCE

def classify_residues(model, sasa_df, cu_coord, his_brace, cu_residue,
                      end_residue=None, sasa_threshold=None):
    classified = {
        "bridging": [],
        "donor": [],
        "acceptor": [],
        }
    all_residues = []
    # First pass: gather all valid residues
    for chain in model:
        for residue in chain:
            resname = residue.get_resname()
            resid = residue.id
            # Restrict to catalytic domain if a boundary is given
            if end_residue is not None and resid[1] > end_residue:
                continue
            full_id = (chain.id, resid, resname)
            new_res = Residue(chain=chain.id, index=residue.id[1], name=residue.get_resname(),
                              sasa=0.0, coords=np.zeros(3), atoms=residue.child_list)
            if resname not in IMPORTANT_ATOMS:
                continue
            if full_id in [r[1] for r in his_brace]:
                continue  # skip histidine brace
            centroid = get_centroid(residue, IMPORTANT_ATOMS[resname])
            if centroid is None:
                continue
            new_res.coords = centroid
            sasa_entry = sasa_df[(sasa_df["resid"] == resid[1])]
            sasa_val = sasa_entry["Relative SASA"].values[0] if not sasa_entry.empty else 0
            new_res.sasa = sasa_val
            all_residues.append(new_res)
    # Second pass: classify residues
    threshold = sasa_threshold if sasa_threshold is not None else EXPOSED_SASA_THRESHOLD
    for res in all_residues:
        # Modified order and classification procedure to set only to one kind
        # Also add exposed residues to donor
        if res.sasa >= threshold:
            classified["donor"].append(res)
        # Add CYS, MET, TRP to bridge residues (skip HIS)
        elif res.name in ["CYS", "MET", "TRP"]:
            classified["bridging"].append(res)
        # For TYR: include if exposed (handled above) OR OH near copper OR has proton acceptor nearby
        # Note: exposed TYR residues are already added to donor in the first if statement
        elif res.name == "TYR":
            # Find the actual residue object in the model for checking
            tyr_residue_obj = None
            for chain in model:
                for residue in chain:
                    if (chain.id == res.chain and residue.id[1] == res.index and 
                        residue.get_resname() == "TYR"):
                        tyr_residue_obj = residue
                        break
                if tyr_residue_obj:
                    break
            
            if tyr_residue_obj:
                # First check if OH is near copper (3.5Å) - if so, include in bridge
                if is_tyr_oh_near_copper(tyr_residue_obj, cu_coord):
                    classified["bridging"].append(res)
                # Otherwise, check for nearby proton acceptor (exposed TYR already handled above)
                elif has_proton_acceptor_near_tyr_oh(tyr_residue_obj, model, his_brace):
                    classified["bridging"].append(res)
                # Otherwise, exclude this tyrosine (not near copper, no proton acceptor, and not exposed)
    
    # Add donor residues to bridge residues (except histidines)
    for donor_res in classified["donor"]:
        if donor_res.name != "HIS":
            # Check if not already in bridging to avoid duplicates
            if donor_res not in classified["bridging"]:
                classified["bridging"].append(donor_res)
    
    # Add copper as acceptor and also to bridge residues
    if cu_residue is not None:
        cu_chain = cu_residue.get_parent().id if cu_residue.get_parent() else " "
        cu_index = cu_residue.id[1] if isinstance(cu_residue.id, tuple) else cu_residue.id
        cu_res = Residue(chain=cu_chain, index=cu_index, 
                        name="CU", sasa=0.0, coords=cu_coord, atoms=cu_residue.child_list)
        classified["acceptor"].append(cu_res)
        # classified["bridging"].append(cu_res)  # Also add copper to bridge residues
    return classified

def print_classified(classified):
    print("Donor:")
    for r in classified["donor"]:
        print("\t", r)
    print("Bridge residues:")
    for r in classified["bridging"]:
        print("\t", r)
    print("Acceptor:")
    for r in classified["acceptor"]:
        print("\t", r)

def print_EH_input(classified, filename, out_dir):
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    def _format(at):
        resname = at.get_parent().resname
        element = at.element
        atom_name = at.name
        chain_id = at.get_parent().get_parent().id
        # Capitalize copper residue name, atom name, and element to "CU"
        if (resname.upper() in METAL_NAMES or element.upper() in METAL_NAMES
                or atom_name.upper() in METAL_NAMES):
            resname = "CU"
            atom_name = "CU"
            element = "CU"
        return "ATOM,{},{},{},{},{},{:5.3f},{:5.3f},{:5.3f},1,{},{}\n".format(at.serial_number,
            atom_name, resname, chain_id,
            at.get_parent().id[1], at.get_coord()[0], at.get_coord()[1], at.get_coord()[2],
            at.bfactor, element)
    
    basename = os.path.join(out_dir, filename)
    
    # Bridge residues - write to bridge file
    output = ""
    for res in classified["bridging"]:
        for at in res.atoms:
            output += _format(at)
    outname = basename[:-4] + "_bridge.csv"
    with open(outname, "w") as bridgef:
        bridgef.write(output)
    
    # Donor - save each donor in a separate file
    # print("Donor:")
    for res in classified["donor"]:
        output = ""
        for at in res.atoms:
            output += _format(at)
        # Create filename with resname and residue index: {basename}_donor_{resname}{resid}.csv
        outname = basename[:-4] + "_donor_{}{}.csv".format(res.name, res.index)
        with open(outname, "w") as donorf:
            donorf.write(output)
    
    # Acceptor (copper)
    # print("Acceptor:")
    output = ""
    for res in classified["acceptor"]:
        for at in res.atoms:
            output += _format(at)
    outname = basename[:-4] + "_acceptor.csv"
    with open(outname, "w") as acceptorf:
        acceptorf.write(output)


if __name__ == "__main__":
    src_lpmo = "1.PDB_structures"
    pdb_files = os.listdir(src_lpmo)
    pdb_files = [os.path.join(src_lpmo, pf) for pf in pdb_files]
    pdb_files.sort()
    # pdb_files = ["1.PDB_structures/050_SmAA10A_8RRY.pdb", "1.PDB_structures/087_MtAA9E_9BJQ.pdb", "1.PDB_structures/037_BaAA10A_5IJU.pdb", "1.PDB_structures/018_NcAA9C_4D7U.pdb"]
    for pdbf in pdb_files:
        structure = load_structure(pdbf)
        sasa_df = compute_relative_sasa(structure)
        for model in structure:
            cu_residue, cu_coord = find_copper(model)
            his_brace = find_histidine_brace(model, cu_coord)
        classified = classify_residues(model, sasa_df, cu_coord, his_brace, cu_residue)
        print("-"*30 + os.path.basename(pdbf) + "-"*30)
        print_classified(classified)
        print_EH_input(classified, os.path.basename(pdbf), "EH_files")
        print("-"*80)
