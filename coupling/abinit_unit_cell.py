from unit_cell import UnitCell
from pathlib import Path
import sys
import os
# To ensure the file is calling from the correct directory. 
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.append(str(current_dir))

class AbinitUnitCell(UnitCell):
    """
    Contains methods that allow the user to manipulate the unit cell and produce Abinit files

    Public Methods: 

    """

    def __init__(self, acell, rprim, coordinates, coord_type, num_atoms, atom_types, znucl, typat, header_path, convergence_path=None, batchScriptHeader_path=None):
        super().__init__(acell, rprim, coordinates, coord_type, num_atoms, atom_types, znucl, typat, header_path, convergence_path)
        self.convergence_path = convergence_path
        self.batchScriptHeader_path = batchScriptHeader_path

    def groundWFKFile(self, output_path, batch_script=False):
        """ 
        Creates the groundWFK.abi file with specified contents.

        Args: 
            output_path (str): path to directory to save the groundWFK.abi file.

        """
# TODO: I don't want the tolwfr to be hardcoded. I think an easy fix around is give users the ability to create their own file. 
        content = """
  ndtset 2

#Set 1 : Ground State Self-Consistent Calculation
#************************************************

  kptopt1 1
  tolvrs 1.0d-18

#Set 2 : Calculation of ddk wavefunctions
#************************************************
  kptopt2 2             # DDK can use only time reversal symmetry
  getwfk2 1             # require ground state wavefunctions from previous run
  rfelfd2 2             # activate DDK perturbation
  iscf2   -3            # this is a non-self-consistent calculation
  tolwfr2 1.0D-18       # tight convergence on wavefunction residuals
"""
        file_path = Path(output_path) / "groundWFK.abi"

        # Write inital content to the output file
        with open(file_path, 'w') as f: 
            f.write(content)

            # Include convergence file contents
            if os.path.isfile(self.convergence_path):
                with open(self.convergence_path, 'r') as conv_file:
                    f.writelines(conv_file.readlines())

            # Append unit cell details
            f.write("\n#Definition of unit cell\n")
            f.write(f"acell {' '.join(map(str, self.acell))}\n")
            f.write(f"rprim\n")
            for coord in self.rprim:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write(f"{self.coord_type}\n")
            for coord in self.coordinates:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write("\n#Definition of atoms\n")
            f.write(f"natom {self.num_atoms}\n")
            f.write(f"ntypat {len(self.atom_types)}\n")
            f.write(f"znucl {' '.join(map(str, self.znucl))}\n")
            f.write(f"typat {' '.join(map(str, self.typat))}\n")

        if batch_script == True: 

            file_path = Path(output_path) / "b-script-groundWFK"

            # Write new batch script file and copy header
            with open(file_path, 'w') as f:
                if os.path.isfile(self.batchScriptHeader_path):
                    with open(self.batchScriptHeader_path, 'r') as batch_file: 
                        f.writelines(batch_file.readlines())
                
        # TODO: Put in sbatch command


    def phononDispFile(self, output_path, batch_script=False):
        """ 
        Creates the phononDispCalc.abi file with specified contents.

        Args: 
            output_path (str): path to directory to save the groundWFK.abi file.

        """
# TODO: the ngqpt is currently hardcoded and needs to be manually calculated. 

        content = """

  ndtset 6

#Definition of q-point grid
#**************************

  nqpt 1     # One qpt for each dataset
  qptopt 1
  ngqpt 4 4 4
  nshiftq 1
  shiftq 0.0 0.0 0.0

iqpt: 5 iqpt+ 1   #automatically iterate through the q pts

#Set 1 : iqpt 1 is the gamma point, so Q=0 phonons and electric field pert.
#**************************************************************************

  getddk1   98         # d/dk wave functions
  kptopt1   2          # Use of time-reversal symmetry
  rfelfd1   3          # Electric-field perturbation response
                       # (in addition to default phonon)

#Sets 2-20 : Finite-wave-vector phonon calculations (defaults for all datasets)
#******************************************************************************

   getwfk  99           # Use GS wave functions
   kptopt  3
   rfphon  1          # Do phonon response
   tolvrs  1.0d-15    # Converge on potential residual

#******
#Flags*
#******

   prtwf 1
   prtden 1
   prtpot 1
   prteig 0
"""
        file_path = Path(output_path) / "phononDispCalc.abi"

        # Write inital content to the output file
        with open(file_path, 'w') as f: 
            f.write(content)

            # Include convergence file contents
            if os.path.isfile(self.convergence_path):
                with open(self.convergence_path, 'r') as conv_file:
                    f.writelines(conv_file.readlines())

            # Append unit cell details
            f.write("\n#Definition of unit cell\n")
            f.write(f"acell {' '.join(map(str, self.acell))}\n")
            f.write(f"rprim\n")
            for coord in self.rprim:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write(f"{self.coord_type}\n")
            for coord in self.coordinates:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write("\n#Definition of atoms\n")
            f.write(f"natom {self.num_atoms}\n")
            f.write(f"ntypat {len(self.atom_types)}\n")
            f.write(f"znucl {' '.join(map(str, self.znucl))}\n")
            f.write(f"typat {' '.join(map(str, self.typat))}\n")

        if batch_script == True: 

            file_path = Path(output_path) / "b-script-groundWFK"

            # Write new batch script file and copy header
            with open(file_path, 'w') as f:
                if os.path.isfile(self.batchScriptHeader_path):
                    with open(self.batchScriptHeader_path, 'r') as batch_file: 
                        f.writelines(batch_file.readlines())
                
        # TODO: Put in sbatch command

    def customFile(self, output_path, batch_script=False): 
        """ 
        Creates the groundWFK.abi file with specified contents. Note, this uses the classes header attribute
        and must be rewritten if you wish for it to be something different. 

        Args: 
            output_path (str): path to directory to save the groundWFK.abi file.

        """
        file_path = Path(output_path) / "groundWFK.abi"

        # Write inital content to the output file
        with open(file_path, 'w') as f: 
            f.write(self.header_path)

            # Include convergence file contents
            if os.path.isfile(self.convergence_path):
                with open(self.convergence_path, 'r') as conv_file:
                    f.writelines(conv_file.readlines())

            # Append unit cell details
            f.write("\n#Definition of unit cell\n")
            f.write(f"acell {' '.join(map(str, self.acell))}\n")
            f.write(f"rprim\n")
            for coord in self.rprim:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write(f"{self.coord_type}\n")
            for coord in self.coordinates:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write("\n#Definition of atoms\n")
            f.write(f"natom {self.num_atoms}\n")
            f.write(f"ntypat {len(self.atom_types)}\n")
            f.write(f"znucl {' '.join(map(str, self.znucl))}\n")
            f.write(f"typat {' '.join(map(str, self.typat))}\n")

        if batch_script == True: 

            file_path = Path(output_path) / "b-script-groundWFK"

            # Write new batch script file and copy header
            with open(file_path, 'w') as f:
                if os.path.isfile(self.batchScriptHeader_path):
                    with open(self.batchScriptHeader_path, 'r') as batch_file: 
                        f.writelines(batch_file.readlines())
                
        # TODO: Put in sbatch command

