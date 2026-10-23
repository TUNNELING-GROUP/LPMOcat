#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#   Copyright 2026 Katarzyna Szleper-Shendy <katarzyna.szleper@polsl.pl> 
# 


import os
import re


def load_protein(pdbin):
    """Generate PyMOL commands to load and style the protein."""
    output = "load ../{}, protein\n".format(pdbin)
    output += "hide cartoon\nshow ribbon, protein\ncolor gray50, protein\n"
    output += "color atomic, not elem C\n"
    output += "zoom protein\n"
    output += "copy cover, protein\nhide everything, cover\n"
    output += "color gray50, cover\nshow mesh, cover\n"
    return output


def add_path_resids(path_resids, path_id, donor_resid, donor_resname, time, distance, scene_name):
    """Generate PyMOL commands to visualize a path with informative scene name and labels."""
    # Use a safe name for selection
    safe_name = f"path_{path_id:0>2}"
    
    output = "sele {}, resi {}\n".format(safe_name, "+".join(map(str, path_resids)))
    output += "color atomic, {}\n".format(safe_name)
    output += "show sticks, {}\n".format(safe_name)
    output += "orient {}\n".format(safe_name)
    
    # Add labels for each residue in the path (one label per residue on CA atom)
    # First, hide all existing labels
    output += "hide labels, all\n"
    # Then show labels for path residues (label only CA atom to get one label per residue)
    for resid in path_resids:
        output += "label resi {} and {} and name CA, resn+resi\n".format(resid, safe_name)
    
    # Create scene with descriptive name (labels will be included in the scene)
    # Use simpler description without special characters to avoid PyMOL crashes
    output += 'scene {}, store, "Donor {} Time {:.2e} Dist {:.1f}A"\n'.format(
        scene_name, f"{donor_resname}{donor_resid}", time, distance)
    
    # After storing the scene, hide the path and labels for next iteration
    output += "color gray50, {}\nhide sticks, {}\n".format(safe_name, safe_name)
    output += "hide labels, {}\n".format(safe_name)
    return output


def parse_paths_file(paths_file):
    """
    Parse a paths.txt file and extract paths with metadata for approximate and exact residence times.
    Returns: (approx_paths, exact_paths)
    where each is a list of tuples: (path, donor_resid, donor_resname, time, distance)
    path is a list of residue IDs
    """
    approx_paths = []
    exact_paths = []
    
    current_section = None
    
    with open(paths_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Detect section headers
            if "approximate mean residence time" in line.lower():
                current_section = "approx"
                continue
            elif "exact mean residence time" in line.lower():
                current_section = "exact"
                continue
            
            # Skip separator lines
            if line.startswith('-') or line.startswith('=') or not line:
                continue
            
            # Parse path lines (format: "   1. [TRP133, HIS92, TYR151, HIS171]  2.1214165919e-05  (donor: TRP133, dist: 17.10 Å)")
            if line and current_section:
                # Extract path, time, donor, and distance
                path_match = re.search(r'\[([^\]]+)\]', line)
                time_match = re.search(r'(\d+\.?\d*[eE][+-]?\d+|\d+\.\d+)', line)
                donor_match = re.search(r'donor:\s*([A-Z]+)(\d+)', line)
                dist_match = re.search(r'dist:\s*([\d.]+)', line)
                
                if path_match:
                    path_str = path_match.group(1)
                    # Parse path as list of residue IDs, extracting numbers from resname+resid
                    path_items = [x.strip() for x in path_str.split(',') if x.strip()]
                    path = []
                    for item in path_items:
                        num_match = re.search(r'(\d+)', item)
                        if num_match:
                            path.append(int(num_match.group(1)))
                    
                    # Extract time, donor, and distance
                    time = float(time_match.group(1)) if time_match else 0.0
                    if donor_match:
                        donor_resname = donor_match.group(1)
                        donor_resid = int(donor_match.group(2))
                    else:
                        donor_resname = "UNK"
                        donor_resid = 0
                    distance = float(dist_match.group(1)) if dist_match else 0.0
                    
                    path_data = (path, donor_resid, donor_resname, time, distance)
                    
                    if current_section == "approx":
                        approx_paths.append(path_data)
                    elif current_section == "exact":
                        exact_paths.append(path_data)
    
    return approx_paths, exact_paths


def create_pymol_script(pdb_file, paths, output_file, time_type):
    """
    Create a PyMOL script for visualizing paths.
    
    Args:
        pdb_file: Path to PDB file
        paths: List of path data tuples: (path, donor_resid, donor_resname, time, distance)
        output_file: Output PyMOL script file
        time_type: "approximate" or "exact" (for naming)
    """
    pymscript = load_protein(pdb_file)
    
    # Add each path
    first_scene_name = None
    for idx, path_data in enumerate(paths, 1):
        path, donor_resid, donor_resname, time, distance = path_data
        # Create scene name: use simple format with donor and index
        # PyMOL scene names should be alphanumeric with underscores
        # We'll use a simple format: path_donor_index
        scene_name = f"path_d{donor_resid}_{idx:02d}"
        
        if first_scene_name is None:
            first_scene_name = scene_name
        
        pymscript += add_path_resids(path, idx, donor_resid, donor_resname, time, distance, scene_name)
    
    # Set initial scene to first path
    if paths and first_scene_name:
        pymscript += f"scene {first_scene_name}\n"
    
    with open(output_file, 'w') as outf:
        outf.write(pymscript)


def process_paths_files(eh_analysis_dir="EH_analysis", pdb_dir="1.PDB_structures", 
                        output_dir="pymol_paths"):
    """
    Process all paths.txt files and create PyMOL scripts for approximate and exact residence times.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all paths.txt files
    paths_files = [f for f in os.listdir(eh_analysis_dir) 
                   if f.endswith('_paths.txt') and os.path.isfile(os.path.join(eh_analysis_dir, f))]
    paths_files.sort()
    
    for paths_file in paths_files:
        paths_path = os.path.join(eh_analysis_dir, paths_file)
        
        # Extract protein ID from filename (e.g., "050_SmAA10A_8RRY_paths.txt" -> "050_SmAA10A_8RRY")
        protein_id = paths_file.replace('_paths.txt', '')
        
        # Find corresponding PDB file
        pdb_file = f"{protein_id}.pdb"
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        if not os.path.isfile(pdb_path):
            print(f"Warning: PDB file not found: {pdb_path}")
            continue
        
        # Parse paths file
        approx_paths, exact_paths = parse_paths_file(paths_path)
        
        # Full path to PDB file for PyMOL script (relative to output directory)
        pdb_path_for_pymol = os.path.join(pdb_dir, pdb_file)
        
        # Create PyMOL script for approximate residence time
        if approx_paths:
            approx_output = os.path.join(output_dir, f"{protein_id}_approx.pml")
            create_pymol_script(pdb_path_for_pymol, approx_paths, approx_output, "approximate")
            print(f"Created: {approx_output} ({len(approx_paths)} paths)")
        
        # Create PyMOL script for exact residence time
        if exact_paths:
            exact_output = os.path.join(output_dir, f"{protein_id}_exact.pml")
            create_pymol_script(pdb_path_for_pymol, exact_paths, exact_output, "exact")
            print(f"Created: {exact_output} ({len(exact_paths)} paths)")


if __name__ == "__main__":
    process_paths_files()
