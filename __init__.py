import subprocess
import os
from pathlib import Path
from .energy import energy_main
from .perturbations import perturbations_main
from .coupling import coupling_main
from .shared import AbinitUnitCell, UnitCell, AbinitFile

__all__ = [
    "energy_main",
    "pert_main",
    "couple_main",
    "AbinitUnitCell",
    "UnitCell",
    "AbinitFile",
]

def generate_boilerplate(input_file):
    """
    
    """
    # Get the current working directory
    current_dir = Path.cwd()

    # Path to the boilerplate script
    script_path = Path(__file__).parent / "boilerplate_generation.sh"

    # Execute the script in the current directory
    result = subprocess.run(
        ["bash", script_path, input_file],
        cwd=current_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    # Handle script output
    if result.returncode == 0:
        print(result.stdout)
    else:
        print(f"Error: {result.stderr}")

def run_loop_smodes_script():
    """
    
    """
    script_path = os.path.join(os.path.dirname(__file__), 'scripts', 'loop_smodes.tcsh')

    # Run script and wiat for it to finish
    result = subprocess.run(['tcsh', script_path], capture_output=True, text=True )



# This is where the program will be executed. Here, the user has 3 choices (so far) as to what they want
def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="FLPZ Program")
    parser.add_argument("program", choices=["energy", "pert", "couple"], help="Select the program to run")
    parser.add_argument("inputs", nargs="*", help="Input files or arguments for the selected program")
    args = parser.parse_args()

    # The main call of the energy program
    if args.program == "energy":
        # Ensure exactly 3 arguments are provided
        if len(args.inputs) != 3:
            print("Error: For 'energy' program, exactly 3 arguments are required: input_file, smodes_input, irrep.")
            sys.exit(1)

        # Assign the 3 arguments to the respective variables
        input_file, smodes_input, irrep = args.inputs

        energy_main(input_file=input_file, smodes_input=smodes_input, irrep=irrep)

    # The main call of the perturbations program
    elif args.program == "pert":
        # Ensure at least 3 argument are provided

        if len(args.inputs) < 3:
            print("Error: for 'perturbations' program, 3 are required: input_file, smodes_input, irrep, piezo (optional)")

        # Assign variables
        input_file, smodes_input, irrep, run_piezo=False = args.input

        perturbations_main(input_file=input_file, smodes_input=smodes_input, irrep=irrep, run_piezo=run_piezo)

    # The main call of the coupling program
    elif args.program == "couple":
        coupling_main(args.inputs)
    else:
        print("Error: Invalid program selected.")
        sys.exit(1)
