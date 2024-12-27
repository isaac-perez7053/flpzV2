#!/usr/bin/env python3
import sys
import argparse
import subprocess
from pathlib import Path
from .shared import AbinitUnitCell, UnitCell, AbinitFile
from .energy import energy_main, parse_inputfile, run_smodes_symmadapt
from .perturbations import perturbations_main
from .coupling import coupling_main


# Add subdirectories to the Python path for importing submodules
current_dir = Path(__file__).parent
sys.path.append(str(current_dir / 'energy'))
sys.path.append(str(current_dir / 'perturbations'))
sys.path.append(str(current_dir / 'coupling'))


class flpz:
    """
    Main class for FLPZ package to study flexo and piezoelectricity.
    """

    @staticmethod
    def energy(*args):
        """Call the energy module with arguments."""
        energy_main(*args)

    @staticmethod
    def perturbations(*args):
        """Call the perturbations module with arguments."""
        perturbations_main(*args)

    @staticmethod
    def coupling(*args):
        """Call the coupling module with arguments."""
        coupling_main(*args)


def main():
    parser = argparse.ArgumentParser(description='FLPZ: A tool to study flexo and piezoelectricity')
    parser.add_argument('inputs', nargs='*', help='Program type (energy, pert, couple) followed by their necessary arguments')

    args = parser.parse_args()
    inputs = args.inputs

    # Validate inputs
    if len(inputs) == 0:
        print("Error: No program type was specified.")
        sys.exit(1)

    program_type = inputs[0]
    additional_args = inputs[1:]  # The remaining arguments

    # Instantiate FLPZ class
    flpz_instance = flpz()

    # Dispatch the appropriate method based on program type
    if program_type == 'energy':
        flpz_instance.energy(*additional_args)
    elif program_type == 'pert':
        flpz_instance.perturbations(*additional_args)
    elif program_type == 'couple':
        flpz_instance.coupling(*additional_args)
    else:
        print(f"Error: Invalid program type '{program_type}'. Available options are: energy, pert, couple.")


if __name__ == "__main__":
    main()