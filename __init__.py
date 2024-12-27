from pathlib import Path
from .energy.energy import energy_main
from .perturbations.perturbations import perturbations_main
from .coupling.coupling import coupling_main

__all__ = [
    "energy_main",
    "pert_main",
    "couple_main",
    "parse_inputfile"
    "run_smodes_symmadapt"
    "AbinitUnitCell",
    "UnitCell",
    "AbinitFile",
]



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

        # Ensure args.input has at least three elements
        if len(args.input) < 3:
            raise ValueError("At least three arguments are required: input_file, smodes_input, irrep")

        # Unpack the first three arguments
        input_file, smodes_input, irrep = args.input[:3]

        # Assign default value for run_piezo if not provided
        run_piezo = args.input[3] if len(args.input) > 3 else False

        energy_main(input_file=input_file, smodes_input=smodes_input, irrep=irrep)

    # The main call of the perturbations program
    elif args.program == "pert":
        # Ensure at least 3 argument are provided

        if len(args.inputs) < 3:
            print("Error: for 'perturbations' program, 3 are required: input_file, smodes_input, irrep, piezo (optional)")

        # Ensure args.input has at least three elements
        if len(args.input) < 3:
            raise ValueError("At least three arguments are required: input_file, smodes_input, irrep")

        # Unpack the first three arguments
        input_file, smodes_input, irrep = args.input[:3]

        # Assign default value for run_piezo if not provided
        run_piezo = args.input[3] if len(args.input) > 3 else False

        perturbations_main(input_file=input_file, smodes_input=smodes_input, irrep=irrep, run_piezo=run_piezo)

    # The main call of the coupling program
    elif args.program == "couple":
        coupling_main(args.inputs)
    else:
        print("Error: Invalid program selected.")
        sys.exit(1)
