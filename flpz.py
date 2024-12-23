#!/usr/bin/env python3 
import sys 
import argparse
import subprocess
from pathlib import Path

# Add subdirectories to the Python path for importign submodules
current_dir = Path(__file__).parent
sys.path.append(str(current_dir / 'energy'))
sys.path.append(str(current_dir / 'perturbations'))
sys.path.append(str(current_dir / 'coupling'))

def main():
    parser = argparse.ArgumentParser(description='FLPZ: A tool to study flexo and piezoelectricity')
    parser.add_argument('inputs', nargs='*', help='Program type (energy, pert, copule) followed by their necessary arguments')

    args = parser.parse_args()
    inputs = args.inputs

     # Validate inputs
    if len(inputs) == 0:
        print("Error: No program type was specified.")
        sys.exit(1)

    program_type = inputs[0]
    additional_args = inputs[1:]  # The remaining arguments

    # Execute the corresponding program
    script_map = {
        'energy': str(current_dir / 'energy' / 'energy.py'),
        'pert': str(current_dir / 'perturbations' / 'perturbations.py'),
        'couple': str(current_dir / 'coupling' / 'coupling.py'),
    }

    if program_type in script_map:
        try:
            subprocess.run(['python3', script_map[program_type], *additional_args], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running {program_type}: {e}")
    else:
        print(f"Error: Invalid program type '{program_type}'. Available options are: energy, pert, couple.")
