#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2026 Carlos Eduardo Sequeiros Borja <csequeiros@polsl.pl>
#  


import os
import subprocess
import glob
import re


def get_donor_files(inpath, pdbin):
    """
    Find all individual donor CSV files for a given PDB file.
    Returns list of tuples: (donor_file_path, resname, resid)
    """
    pdb_base = pdbin.replace(".pdb", "")
    pattern = os.path.join(inpath, "{}_donor_*.csv".format(pdb_base))
    donor_files = glob.glob(pattern)
    
    donor_info = []
    for donor_file in donor_files:
        # Extract resname and resid from filename: {pdb}_donor_{resname}{resid}.csv
        basename = os.path.basename(donor_file)
        match = re.search(r'_donor_([A-Z]+)(\d+)\.csv$', basename)
        if match:
            resname = match.group(1)
            resid = match.group(2)
            donor_info.append((donor_file, resname, resid))
    
    return donor_info

def _execute_command(command, shell, cwd):
    with subprocess.Popen(
            command,
            cwd=cwd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE) as process:
        out, err = process.communicate()
        errcode = process.returncode
    if errcode == 0:
        return out.decode("utf-8")
    message = "Shell command failed: \n"
    message += " ".join(command)
    message += "\nerror code:"
    message += str(errcode)
    if out is not None:
        message += "\nstdout:\n"
        message += out.decode("utf-8")
    if err is not None:
        message += "\nstderr:\n"
        message += err.decode("utf-8")
    raise RuntimeError(message)

def run_EH(inpath, pdbin, don_csv, resname, resid, outpath):
    """
    Run EHPath for a specific donor file.
    
    Args:
        inpath: Path to input directory
        pdbin: PDB filename
        don_csv: Path to the specific donor CSV file
        resname: Residue name of the donor
        resid: Residue ID of the donor
        outpath: Path to output directory
    """
    ace_csv = os.path.join(inpath, pdbin.replace(".pdb", "_acceptor.csv"))
    bri_csv = os.path.join(inpath, pdbin.replace(".pdb", "_bridge.csv"))
    # Use the specific donor file instead of a generic one
    outfile = os.path.join(outpath, pdbin.replace(".pdb", "_donor_{}{}.log".format(resname, resid)))
    if os.path.isfile(outfile):
        return
    comms = ["python3", "EHPath/EHPath.py",
             don_csv, bri_csv, ace_csv, resid]
    print("Running:", " ".join(comms))
    output = _execute_command(comms, False, ".")
    os.makedirs(outpath, exist_ok=True)
    with open(outfile, "w") as outf:
        outf.write(output)


if __name__ == "__main__":
    lpmo_pdbs = "1.PDB_structures"
    pdb_files = os.listdir(lpmo_pdbs)
    pdb_files.sort()
    csv_lpmo = "EH_files"
    eh_results = "EH_results"
    # pdb_files = ["050_SmAA10A_8RRY.pdb", "087_MtAA9E_9BJQ.pdb", "037_BaAA10A_5IJU.pdb", "018_NcAA9C_4D7U.pdb"]
    for pdbf in pdb_files:
        donor_files = get_donor_files(csv_lpmo, pdbf)
        for don_csv, resname, resid in donor_files:
            try:
                run_EH(csv_lpmo, pdbf, don_csv, resname, resid, eh_results)
            except RuntimeError as e:
                print(f"Error running EHPath for {pdbf}, donor {resname}{resid}: {e}")
                exit(1)

