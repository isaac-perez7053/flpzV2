#TODO: I need to fix the imports as there is some confusion as to how to get the correct files
import subprocess 
import os
import re
from pathlib import Path
from shared import AbinitFile, AbinitUnitCell, UnitCell
import shared.smodes_postproc_abinit, shared.smodes_symmadapt_abinit
from .. import generate_boilerplate, run_loop_smodes_script

def parse_inputfile(input_file):
    global num_datapoints, name, genstruc, min_amp, max_amp, sbatch_preamble

    with open(input_file, 'r') as f: 
        lines = f.readlines()

        
    # Extract num_datapoints
    num_datapoints = None
    for i, line in enumerate(lines):
        if line.strip().startswith('num_datapoints'):
            match = re.search(r"\d+", line)
            if match:
                num_datapoints = int(match.group())
                break  # Exit the loop once 'num_datapoints' is found

    # If num_datapoints is not set, raise an exception
    if num_datapoints is None:
        raise Exception("Number of datapoints is missing in the input file!")


    # Extract name
    name = None
    for i, line in enumerate(lines):
        if line.strip().startswith('name'):
            match = re.search(r"\b[a-zA-Z0-9]+\b", line.strip())
            if match:
                name = match.group()
                break  

    if name is None:
        raise Exception("Name is missing in the input file!")
        
    # Extract genstruc abinit file name
    genstruc = None
    for i, line in enumerate(lines):
        if line.strip().startswith('genstruc'):
            match = re.search(r"\b[a-zA-Z]+\b", line.strip())
            if match:
                genstruc = match.group()
                break 
        
    if genstruc is None:
        raise Exception("The Abinit file (genstruc) is missing in the input file!")
        
    # Extract the minimum amplitude 
    min_amp = None
    for i, line in enumerate(lines):
        if line.strip().startswith('min'):
            match = re.search(r"^\s*[-+]?[0-9]*\.?[0-9]+", line.strip())
            if match:
                min_amp = match.group()
                break
        
    if min_amp is None:
        raise Exception("The min_amp (min) is missing in the input file!")
        
    # Extract the maximum amplitude 
    max_amp = None
    for i, line in enumerate(lines):
        if line.strip().startswith('max'):
            match = re.search(r"^\s*[-+]?[0-9]*\.?[0-9]+", line.strip())
            if match: 
                max_amp = match.group()
                break
        
    if max_amp is None: 
        raise Exception("the max_amp (max) is missing in the input file!")
        
    # Extract the sbatch preamble 
    sbatch_preamble = None
    for i, line in enumerate(lines):
        if line.strip().startswith('sbatch_preamble'):
            match = re.search(r"^[\w,\s-]+\.[A-Za-z]{2,}$", line.strip())
            if match: 
                sbatch_preamble = match.group()
                break
        
    if sbatch_preamble is None:
        raise Exception("sbatch preamble is missing in the input file!")


# TODO: The path the file isn't dynamic and msut be fixed. 
def run_smodes_symmadapt(smodes_input, irrep):
    """
    Runs the smodes_symmadapt_abinit script dynamically located relative to the shared directory,
    ensuring all output files are placed in the current working directory.

    Args:
        smodes_input (str): Path to the SMODES input file.
        irrep (str): Irreducible representation argument.

    Returns:
        str: The stdout of the script execution.
    
    Raises:
        RuntimeError: If the script execution fails.
    """
    # Locate the smodes_symmadapt_abinit script dynamically
    script_dir = Path(__file__).resolve().parent / "shared"
    smodes_script = script_dir / "smodes_symmadapt_abinit.py"

    if not smodes_script.exists():
        raise FileNotFoundError(f"Script not found at {smodes_script}")

    # Command to execute the smodes_symmadapt_abinit script
    cmd = [
        "python3",  # Or the appropriate Python executable
        str(smodes_script),
        str(smodes_input),
        str(irrep),
    ]

    # Run the script and ensure output files are placed in the current working directory
    try:
        result = subprocess.run(
            cmd,
            cwd=Path.cwd(),  # Set the working directory to the current directory
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        print(f"smodes_symmadapt_abinit completed successfully:\n{result.stdout}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        error_msg = f"Error running smodes_symmadapt_abinit:\n{e.stderr}"
        print(error_msg)
        raise RuntimeError(error_msg)


def energy_main(input_file, smodes_input, irrep):
    print("Energy program running")
    
    # Parse input file and store variables
    parse_inputfile(input_file=input_file)

    # Store variables in the Abinit file
    abinit_file = AbinitFile(filepath = genstruc)

    # Generate boilerplate
    generate_boilerplate(str(genstruc))

    # Run smodes_symmadapt_abinit
    try:
        run_smodes_symmadapt(smodes_input, irrep)
    except RuntimeError as e:
        print(f"Failed to run smodes_symmadapt_abinit: {e}")
        raise

    # Run loop_smodes.tcsh
    

    




    abinit_file.unit_cell



