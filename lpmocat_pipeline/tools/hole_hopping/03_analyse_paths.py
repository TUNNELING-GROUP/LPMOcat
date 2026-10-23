#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#   Copyright 2026 Katarzyna Szleper-Shendy <katarzyna.szleper@polsl.pl> 
# 


import os
import re
import numpy as np
from collections import defaultdict
from Bio.PDB import PDBParser


def get_copper_coords(acceptor_csv):
    """
    Read copper coordinates from acceptor CSV file.
    Returns (x, y, z) tuple or None if not found.
    """
    if not os.path.isfile(acceptor_csv):
        return None
    
    with open(acceptor_csv, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 9 and parts[2] == 'CU':  # atomtype is CU
                try:
                    x = float(parts[6])
                    y = float(parts[7])
                    z = float(parts[8])
                    return np.array([x, y, z])
                except (ValueError, IndexError):
                    continue
    return None


def get_residue_type(pdb_file, resid):
    """
    Get the residue type (resname) for a given residue ID.
    Returns resname string or None if not found.
    """
    if not os.path.isfile(pdb_file):
        return None
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        resid_int = int(resid)
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.id[1] if isinstance(residue.id, tuple) else residue.id
                    if res_id == resid_int:
                        return residue.get_resname()
        return None
    except Exception as e:
        return None


def get_residue_object(pdb_file, resid):
    """
    Get the residue object for a given residue ID.
    Returns BioPython residue object or None if not found.
    """
    if not os.path.isfile(pdb_file):
        return None
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        resid_int = int(resid)
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.id[1] if isinstance(residue.id, tuple) else residue.id
                    if res_id == resid_int:
                        return residue
        return None
    except Exception as e:
        return None


def is_tyr_oh_near_copper(tyr_residue, cu_coord, distance_threshold=3.5):
    """
    Check if a tyrosine residue's phenolic OH is within distance_threshold of copper.
    
    Args:
        tyr_residue: BioPython residue object for tyrosine
        cu_coord: numpy array with copper coordinates
        distance_threshold: Distance threshold in Angstroms (default 3.5)
    
    Returns:
        bool: True if OH is within threshold of copper, False otherwise
    """
    if cu_coord is None or tyr_residue is None:
        return False
    
    # Get the OH atom from tyrosine (phenolic oxygen)
    if "OH" not in tyr_residue:
        return False
    
    oh_atom = tyr_residue["OH"]
    oh_coord = oh_atom.get_coord()
    
    # Calculate distance to copper
    dist = np.linalg.norm(oh_coord - cu_coord)
    return dist <= distance_threshold


def get_donor_residue_coords(pdb_file, donor_resid):
    """
    Get the centroid coordinates of important atoms in the donor residue.
    Returns (x, y, z) numpy array or None if not found.
    """
    if not os.path.isfile(pdb_file):
        return None
    
    # Important atoms for different residue types (from 01_KS_classify_resids.py)
    IMPORTANT_ATOMS = {
        "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "NE1", "CD1", "CG"],
        "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
        "CYS": ["SG"],
        "MET": ["SD"]
    }
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        donor_resid_int = int(donor_resid)
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check if this is the donor residue
                    res_id = residue.id[1] if isinstance(residue.id, tuple) else residue.id
                    if res_id == donor_resid_int:
                        resname = residue.get_resname()
                        
                        # Get important atoms for this residue type
                        atom_names = IMPORTANT_ATOMS.get(resname, ["CA"])  # Default to CA if not in list
                        
                        # Get coordinates of important atoms
                        coords = []
                        for atom_name in atom_names:
                            if atom_name in residue:
                                coords.append(residue[atom_name].get_coord())
                        
                        if coords:
                            # Return centroid
                            return np.mean(coords, axis=0)
                        else:
                            # Fallback to CA if no important atoms found
                            if "CA" in residue:
                                return residue["CA"].get_coord()
        return None
    except Exception as e:
        print(f"Warning: Could not read donor residue {donor_resid} from {pdb_file}: {e}")
        return None


def calculate_distance(coord1, coord2):
    """
    Calculate Euclidean distance between two coordinate arrays.
    """
    if coord1 is None or coord2 is None:
        return None
    return np.linalg.norm(coord1 - coord2)


def filter_his_donor_paths(paths, pdb_file, cu_coord):
    """
    Filter paths: 
    1. Exclude paths with only 2 residues (length 2)
    2. If histidine is the donor (first residue), 
       the second residue must be tyrosine, otherwise exclude the path.
    3. If histidine is donor and second residue is tyrosine,
       exclude the path if that tyrosine's OH is within 3.5Å of copper.
    
    Args:
        paths: List of tuples (path, time) where path is [resid1, resid2, ...]
        pdb_file: Path to PDB file to look up residue types
        cu_coord: numpy array with copper coordinates
    
    Returns:
        Filtered list of paths
    """
    filtered_paths = []
    
    for path, time in paths:
        # Exclude paths with only 2 residues
        if len(path) <= 2:
            continue
        
        donor_resid = path[0]
        second_resid = path[1]
        
        donor_type = get_residue_type(pdb_file, donor_resid)
        second_type = get_residue_type(pdb_file, second_resid)
        
        # If donor is histidine, second residue must be tyrosine
        if donor_type == "HIS":
            if second_type == "TYR":
                # Check if this tyrosine is close to copper - if so, exclude the path
                tyr_residue = get_residue_object(pdb_file, second_resid)
                if is_tyr_oh_near_copper(tyr_residue, cu_coord, distance_threshold=3.5):
                    # Exclude this path: HIS donor with TYR that is close to copper
                    continue
                else:
                    filtered_paths.append((path, time))
            # Otherwise, exclude this path (second residue is not TYR)
        else:
            # Donor is not histidine, include the path
            filtered_paths.append((path, time))
    
    return filtered_paths


def parse_log_file(log_path):
    """
    Parse a log file and extract paths with their times.
    Returns a tuple: (approx_paths, exact_paths)
    where each is a list of tuples: [(path, time), ...]
    """
    approx_paths = []
    exact_paths = []
    
    with open(log_path, 'r') as f:
        lines = f.readlines()
    
    current_section = None
    for line in lines:
        line = line.strip()
        
        # Detect section headers
        if "approximate mean residence time" in line.lower():
            current_section = "approx"
            continue
        elif "exact mean residence time" in line.lower():
            current_section = "exact"
            continue
        
        # Parse path lines (format: [12, 3, 93, 94, 160, 179] 1.1937535809875595e-05)
        if line.startswith('[') and current_section:
            # Extract path and time
            match = re.match(r'\[(.*?)\]\s+(.+)', line)
            if match:
                path_str = match.group(1)
                time_str = match.group(2)
                
                # Parse path as list of integers
                # Handle np.int64(129) format from newer numpy versions
                def _parse_int(s):
                    s = s.strip()
                    m = re.match(r'np\.int\d+\((\d+)\)', s)
                    return int(m.group(1)) if m else int(s)
                path = [_parse_int(x) for x in path_str.split(',')]
                # Parse time as float
                time = float(time_str)
                
                if current_section == "approx":
                    approx_paths.append((path, time))
                elif current_section == "exact":
                    exact_paths.append((path, time))
    
    return approx_paths, exact_paths


def analyze_eh_results(eh_results_dir, output_dir=None, eh_files_dir="EH_files", pdb_dir="1.PDB_structures", protein_filter=None):
    """
    Analyze EH result files, group by protein, and write sorted paths to files.
    If protein_filter is given (e.g. '004_AfAA9_7OVA'), only process that protein's logs.
    """
    if output_dir is None:
        output_dir = "EH_analysis"

    os.makedirs(output_dir, exist_ok=True)

    # Dictionary to store paths by protein
    # Structure: {protein_id: {'approx': [(path, time, donor_resid), ...],
    #                           'exact': [(path, time, donor_resid), ...]}}
    protein_data = defaultdict(lambda: {'approx': [], 'exact': []})

    # Process log files (optionally filtered to a single protein)
    log_files = [f for f in os.listdir(eh_results_dir)
                 if f.endswith('.log') and os.path.isfile(os.path.join(eh_results_dir, f))]
    if protein_filter:
        log_files = [f for f in log_files if f.startswith(protein_filter + '_donor_')]
    log_files.sort()
    
    for log_file in log_files:
        log_path = os.path.join(eh_results_dir, log_file)
        
        # Extract protein ID and donor resid from filename
        # Format: 037_BaAA10A_5IJU_donor_HIS77.log or 037_BaAA10A_5IJU_donor_TYR12.log
        base_name = log_file.replace('.log', '')
        
        # Try to match the new format: {protein_id}_donor_{resname}{resid}
        match = re.search(r'^(.+)_donor_([A-Z]+)(\d+)$', base_name)
        if match:
            protein_id = match.group(1)
            donor_resname = match.group(2)
            donor_resid = match.group(3)
        else:
            # Fallback: try old format or simple split
            parts = base_name.rsplit('_', 1)
            if len(parts) == 2:
                protein_id = parts[0]
                donor_resid = parts[1]
            else:
                protein_id = base_name
                donor_resid = "unknown"
        
        # Parse the log file
        approx_paths, exact_paths = parse_log_file(log_path)
        
        # Get PDB file path and copper coordinates for filtering
        pdb_file = os.path.join(pdb_dir, f"{protein_id}.pdb")
        acceptor_csv = os.path.join(eh_files_dir, f"{protein_id}_acceptor.csv")
        cu_coords = get_copper_coords(acceptor_csv)
        
        # Filter paths: if HIS is donor, second residue must be TYR, and TYR must not be close to copper
        approx_paths = filter_his_donor_paths(approx_paths, pdb_file, cu_coords)
        exact_paths = filter_his_donor_paths(exact_paths, pdb_file, cu_coords)
        
        # Add paths with donor resid info
        for path, time in approx_paths:
            protein_data[protein_id]['approx'].append((path, time, donor_resid))
        
        for path, time in exact_paths:
            protein_data[protein_id]['exact'].append((path, time, donor_resid))
    
    # Write output files for each protein
    for protein_id, data in protein_data.items():
        print(f"Processing protein: {protein_id} (approx paths: {len(data['approx'])}, exact paths: {len(data['exact'])})")
        output_file = os.path.join(output_dir, f"{protein_id}_paths.txt")
        
        # Get copper coordinates for this protein
        acceptor_csv = os.path.join(eh_files_dir, f"{protein_id}_acceptor.csv")
        cu_coords = get_copper_coords(acceptor_csv)
        
        # Get PDB file path
        pdb_file = os.path.join(pdb_dir, f"{protein_id}.pdb")
        
        # Cache distances per donor resid
        distance_cache = {}
        resname_cache = {}
        
        with open(output_file, 'w') as f:
            f.write(f"Path Analysis for Protein: {protein_id}\n")
            f.write("=" * 80 + "\n\n")
            
            # Write approximate mean residence time paths
            f.write("Paths ranked according to approximate mean residence time:\n")
            f.write("-" * 80 + "\n")
            
            # Sort by time
            sorted_approx = sorted(data['approx'], key=lambda x: x[1])
            
            for idx, (path, time, donor_resid) in enumerate(sorted_approx, 1):
                # Build path string with residue names
                path_str_parts = []
                for resid in path:
                    if resid not in resname_cache:
                        rn = get_residue_type(pdb_file, resid)
                        if rn is None:
                            rn = "CU" if resid == 999 else "UNK"
                        resname_cache[resid] = rn
                    resname = resname_cache[resid]
                    path_str_parts.append(f"{resname}{resid}")
                path_str = '[' + ', '.join(path_str_parts) + ']'

                # Get donor resname
                if donor_resid not in resname_cache:
                    resname_cache[donor_resid] = get_residue_type(pdb_file, donor_resid) or "UNK"
                donor_resname = resname_cache[donor_resid]

                # Calculate distance if not cached
                if donor_resid not in distance_cache:
                    donor_coords = get_donor_residue_coords(pdb_file, donor_resid)
                    if cu_coords is not None and donor_coords is not None:
                        distance_cache[donor_resid] = calculate_distance(cu_coords, donor_coords)
                    else:
                        distance_cache[donor_resid] = None
                
                distance = distance_cache[donor_resid]
                if distance is not None:
                    f.write(f"{idx:4d}. {path_str:70s} {time:.10e}  (donor: {donor_resname}{donor_resid}, dist: {distance:.2f} Å)\n")
                else:
                    f.write(f"{idx:4d}. {path_str:70s} {time:.10e}  (donor: {donor_resname}{donor_resid})\n")
            
            f.write("\n\n")
            
            # Write exact mean residence time paths
            f.write("Paths ranked according to exact mean residence time:\n")
            f.write("-" * 80 + "\n")
            
            # Sort by time
            sorted_exact = sorted(data['exact'], key=lambda x: x[1])
            
            for idx, (path, time, donor_resid) in enumerate(sorted_exact, 1):
                # Build path string with residue names
                path_str_parts = []
                for resid in path:
                    if resid not in resname_cache:
                        rn = get_residue_type(pdb_file, resid)
                        if rn is None:
                            rn = "CU" if resid == 999 else "UNK"
                        resname_cache[resid] = rn
                    resname = resname_cache[resid]
                    path_str_parts.append(f"{resname}{resid}")
                path_str = '[' + ', '.join(path_str_parts) + ']'

                # Get donor resname
                if donor_resid not in resname_cache:
                    resname_cache[donor_resid] = get_residue_type(pdb_file, donor_resid) or "UNK"
                donor_resname = resname_cache[donor_resid]

                # Use cached distance or calculate if not cached
                if donor_resid not in distance_cache:
                    donor_coords = get_donor_residue_coords(pdb_file, donor_resid)
                    if cu_coords is not None and donor_coords is not None:
                        distance_cache[donor_resid] = calculate_distance(cu_coords, donor_coords)
                    else:
                        distance_cache[donor_resid] = None
                
                distance = distance_cache[donor_resid]
                if distance is not None:
                    f.write(f"{idx:4d}. {path_str:70s} {time:.10e}  (donor: {donor_resname}{donor_resid}, dist: {distance:.2f} Å)\n")
                else:
                    f.write(f"{idx:4d}. {path_str:70s} {time:.10e}  (donor: {donor_resname}{donor_resid})\n")
        
        print(f"Written: {output_file} ({len(sorted_approx)} approx paths, {len(sorted_exact)} exact paths)")


if __name__ == "__main__":
    eh_results_dir = "EH_results"
    analyze_eh_results(eh_results_dir)
