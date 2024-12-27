#TODO: I need to fix the imports as there is some confusion as to how to get the correct files
import subprocess 
import os
import re
from pathlib import Path
from shared.boilerplate_generation import BoilerplateGenerator
from shared.smodes_processor import SmodesProcessor


def parse_inputfile(filepath):
    """
    Returns the variables contained in the flpz input file
    """
    with open(filepath, 'r') as f: 
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
    for line in lines:
        if line.strip().startswith('name'):
            match = re.search(r"name\s+([a-zA-Z0-9]+)", line)
            if match:
                name = match.group(1)  # Get the alphanumeric value after 'name'
                break

    if name is None:
        raise Exception("Name is missing in the input file!")
        
    # Extract genstruc abinit file name
    genstruc = None
    for line in lines:
        if line.strip().startswith('genstruc'):
            match = re.search(r"genstruc\s+([a-zA-Z0-9_.-]+)", line)
            if match:
                genstruc = match.group(1)  # Capture the filename after 'genstruc'
                break

    if genstruc is None:
        raise Exception("The Abinit file (genstruc) is missing in the input file!")

        
    # Extract the minimum amplitude
    min_amp = None
    for line in lines:
        if line.strip().startswith('min'):
            match = re.search(r"min\s+([-+]?\d*\.?\d+)", line)
            if match:
                min_amp = match.group(1)  # Get the numeric value after 'min'
                break

    if min_amp is None:
        raise Exception("The min_amp (min) is missing in the input file!")

    # Extract the maximum amplitude
    max_amp = None
    for line in lines:
        if line.strip().startswith('max'):
            match = re.search(r"max\s+([-+]?\d*\.?\d+)", line)
            if match:
                max_amp = match.group(1)  # Get the numeric value after 'max'
                break

    if max_amp is None:
        raise Exception("The max_amp (max) is missing in the input file!")
        
    # Extract the sbatch preamble
    sbatch_preamble = None
    for line in lines:
        if line.strip().startswith('sbatch_preamble'):
            match = re.search(r"sbatch_preamble\s+([a-zA-Z0-9_.-]+)", line)
            if match:
                sbatch_preamble = match.group(1)  # Capture the filename after 'sbatch_preamble'
                break

    if sbatch_preamble is None:
        raise Exception("sbatch preamble is missing in the input file!")

    
    return num_datapoints, str(name), str(genstruc), min_amp, max_amp, str(sbatch_preamble)


def energy_main(*args):
    print("Energy program running")
    from .. import generate_boilerplate, run_loop_smodes_script, AbinitFile   # Delayed import
    
    
    if len(args) < 3:
        raise ValueError("Missing required arguments: input_file, smodes_input, irrep")

    input_file, smodes_input, irrep = args[:3]
    run_piezo = args[3] if len(args) > 3 else False


    # Parse input file and store variables
    num_datapoints, name, genstruc, min_amp, max_amp, sbatch_preamble = parse_inputfile(input_file=input_file)

    # Store variables in the Abinit file
    smodes_file = SmodesProcessor()

    # Generate boilerplate


    # Detect instabilities
    pass

    # Input unit cell from abinit file with eigenvector and use a subclass to generate new instances of the class that will be used to write new files. 

    # Detect if unstable modes exist. 


if __name__ == "__main__":
    energy_main()    

    







