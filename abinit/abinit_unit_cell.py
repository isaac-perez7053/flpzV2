from flpz.unit_cell_module import UnitCell
import numpy as np
import os
import re
import shutil
import sys
import tempfile
import subprocess
import time
import copy
from pathlib import Path

# To ensure the file is calling from the correct directory. 
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.append(str(current_dir))

class AbinitUnitCell(UnitCell):
    """
    Defines the AbinitUnitCell class, which encapsulates the necessary information for an Abinit simulation.

    This class extends the `UnitCell` class to provide additional functionality specifically for handling 
    Abinit input files and running simulations. 

    Attributes:
        abi_file (str): The path to the main Abinit file containing unit cell information.
        file_name (str): The name of the file derived from `abi_file`.
        ecut (int): Energy cutoff value.
        ecutsm (float): Smearing on the energy cutoff.
        nshiftk (int): Number of shifts applied to k-points.
        shiftk (list of float): Shift vectors for k-points.
        nstep (int): Number of steps for SCF calculation.
        diemac (float): Macroscopic dielectric constant.
        ixc (int): Exchange-correlation functional index.
        pp_dirpath (str): Directory path for pseudopotentials.
        pseudos (str): Pseudopotential filename(s).
        kptrlatt (list of list of int): Matrix defining reciprocal space vectors.
        nband (int): Number of bands for electronic structure calculation.
        toldfe (str): Tolerance on the total energy difference in calculations.
        convergence_path (str): Path to a file containing convergence parameters.
        batchScriptHeader_path (str): Path to a batch script header file.
        runningJobs (list of int): List of job IDs for currently running jobs.
        energy (float): Energy of the current configuration.

    Public Methods:
        findSpaceGroup(): Determines and returns the space group of the unit cell.
        convertToXcart(): Converts and returns Cartesian coordinates of the unit cell.
        convertToXred(): Converts and returns reduced coordinates of the unit cell.
        write_ground_workfunction_file(output_path): Creates an Abinit file for calculating work function.
        write_phonon_dispersion_file(output_path): Creates an Abinit file for phonon dispersion calculation.
        write_custom_abifile(output_file, header_file, toldfe=True): Writes a custom Abinit .abi file.
        all_jobs_finished(): Checks if all submitted jobs have finished.
        wait_for_jobs_to_finish(check_time=60): Waits until all submitted jobs are completed.
        run_abinit(input_file='abinit', batch_name='abinit_job', ...): Runs Abinit with specified settings.
        write_batch_script(batch_script_header_file='default_batch_file', ...): Writes a batch script for running simulations.
        perturbations(pert): Applies perturbations to unit cell coordinates and returns a new instance.
        copy_abinit_unit_cell(): Creates a deep copy of the current instance.
        run_energy_calculation(host_spec='mpirun -hosts=localhost -np 30'): Executes an energy calculation for the unit cell.
        grab_energy(abo_file=None): Retrieves the total energy from an Abinit output file.
        change_coordinates(new_coordinates, cartesian=False, reduced=False): Updates the coordinates of the unit cell.
    """
    # TODO: Make a method that will place all pseudopotentials into the pseudopotential folder in the program. All pp_dirpaths will then be the same

    def __init__(self, abi_file, convergence_path=None, batch_script_header_file=None):
        """
        Initializes an instance of AbinitUnitCell.

        Args:
            abi_file (str): Path to the Abinit input file containing unit cell details.
            convergence_path (str, optional): Path to a file with convergence parameters. Defaults to None.
            batch_script_header_file (str, optional): Path to a batch script header file. Defaults to None.
        """
        # Call the parent class's initializer with the keyword arguments
        self.abi_file = str(abi_file)
        self.file_name = self.abi_file.replace('.abi', '')

        super().__init__(abi_file=abi_file)

        # Convergence attributes
        self.ecut = None
        self.ecutsm = None
        self.ecutsm = None
        self.nshiftk = None
        self.shiftk = []
        self.nstep = None
        self.diemac = None
        self.ixc = None
        self.pp_dirpath = None
        self.pseudos = None
        self.kptrlatt = None
        self.nband = None
        self.toldfe = None

        # Initialize additional attributes specific to AbinitUnitCell
        if convergence_path is None: 
            self._initialize_convergence_from_file()
        self.convergence_path = convergence_path
        self.batchScriptHeader_path = batch_script_header_file
        self.runningJobs = []

        # Other attributes that need to be calculated after initialization
        self.energy = None
        self.piezo_tensor_clamped = None
        self.piezo_tensor_relaxed = None
        self.flexo_tensor = None
   

    #--------------------------
    # Initialization Methods
    #--------------------------

    def _initialize_convergence_from_file(self):
        """
        Extracts convergence features from the Abinit file by parsing the required parameters.
        Raises exceptions if critical parameters are missing from the file.
        """

        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            # Copy the contents of the original file to the temporary file
            shutil.copyfile(self.abi_file, temp_file.name)
            
            # Open the temporary file for reading
            with open(temp_file.name, 'r') as f:
                lines = f.readlines()

        # Extract ecut
        for i, line in enumerate(lines):
            if line.strip().startswith('ecut'):
                match = re.search(r"\d+", line)
                if match:
                    ecut = int(match.group())
                del lines[i]
                break

        
        self.ecut = ecut
        if self.ecut is None:
            raise Exception("ecut is missing in the Abinit file!")
        
        # Extract ecutsm
        for i, line in enumerate(lines):
            if line.strip().startswith('ecutsm'):
                match = re.search(r"\d+\.\d+|\d+", line)
                if match: 
                    ecustsm = float(match.group())
                del lines[i]
                break

        self.ecutsm = ecustsm

        # Extract nshiftk
        for i, line in enumerate(lines): 
            if line.strip().startswith('nshiftk'):
                match = re.search(r"\d+", line)
                if match:
                    nshiftk = int(match.group())
                del lines[i]
                break
        
        self.nshiftk = nshiftk
        if self.nshiftk is None:
            raise Exception("nshiftk is missing in the Abinit file!")

        # Extract shiftk
        for i, line in enumerate(lines):
            if line.strip().startswith('shiftk'):
                # Extract different features of the acell feature and map it to a float where it will be inserted into a list.
                match = re.search(r"(\d+)\*([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*)", line)
                if match:
                    count = int(match.group(1))
                    value = float(match.group(2))
                    shiftk = [value] * count
                else:
                    shiftk = list(map(float, re.findall(r"[-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*", line)))
                # Delete extracted lines in copy
                del lines[i]

        self.shiftk = shiftk
        if self.shiftk is None:
            raise Exception("shiftk is missing in the Abinit file!")
        
        # Extract nband
        for i, line in enumerate(lines):
            if line.strip().startswith('nband'):
                match = re.search(r"\d+", line)
                if match:
                    nband = int(match.group())
                del lines[i]
                break
        
        self.nband = nband
        if self.nband is None:
            raise Exception("nband is missing in the Abinit file!")
        
        # Extract nstep 
        for i, line in enumerate(lines):
            if line.strip().startswith('nstep'):
                match = re.search(r"\d+", line)
                if match:
                    nstep = int(match.group())
                del lines[i]
                break
        
        self.nstep = nstep
        if self.nstep is None:
            raise Exception("ecut is missing in the Abinit file!")
        
        # Extract diemac 
        for i, line in enumerate(lines):
            if line.strip().startswith('diemac'):
                match = re.search(r"\d+\.\d+|\d+", line)
                if match:
                    diemac = float(match.group())
                del lines[i]
                break
        
        self.diemac = diemac
        if self.diemac is None:
            self.diemac = 4.0

        # Extract toldfe
        for i, line in enumerate(lines):
            if line.strip().startswith('toldfe'):
                match = re.search(r"toldfe\s+(\d+\.\d+d[+-]?\d+)", line)
                if match:
                    toldfe = match.group(1)
                del lines[i]
                break
        
        self.toldfe = toldfe
        
        # Extract ixc
        for i, line in enumerate(lines):
            if line.strip().startswith('ixc'):
                match = re.search(r"[-+]?\d+", line)
                if match:
                    ixc = int(match.group())
                del lines[i]
                break
        
        self.ixc = ixc
        if self.ixc is None:
            raise Exception("ixc is missing in the Abinit file!")
        
        # Extract pp_dirpath
        for i, line in enumerate(lines):
            if line.strip().startswith('pp_dirpath'):
                match = re.search(r'pp_dirpath\s+"([^"]+)"', line)
                if match: 
                    pp_dirpath = str(match.group(1))
                del lines[i]
                break
    
        self.pp_dirpath = pp_dirpath
        if self.pp_dirpath is None:
            raise Exception("pp_dirpath is missing in the Abinit file!")
        
        # Extract pseudos
        for i, line in enumerate(lines):
            if line.strip().startswith('pseudos'):
                match = re.search(r'pseudos\s+"([^"]+)"', line)
                if match: 
                    pseudos = str(match.group(1))
                del lines[i]
                break
        
        self.pseudos = pseudos
        if self.pseudos is None:
            raise Exception("pseudos is missing in the Abinit file!")
        
        # Extract kptrlatt
        kptrlatt = []
        for i, line in enumerate(lines):
            if line.strip() == 'kptrlatt':
                del lines[i]
                j = i
                while j < len(lines) and re.match(r"^\s*[-+]?\d+", lines[j]):
                    kptrlatt.append(list(map(int, lines[j].split())))
                    del lines[j]
                break
        
        self.kptrlatt = kptrlatt

        os.remove(temp_file.name)


    def findSpaceGroup(self):
        """
        Determines and returns the space group of the unit cell using parent class functionality.

        Returns:
            str: The identifier or symbol of the space group.
        """
        return super().findSpaceGroup()
    
    def convertToXcart(self):
        """
        Converts and returns Cartesian coordinates of the unit cell using parent class functionality.

        Returns:
            np.ndarray: Cartesian coordinates of the unit cell.
        """
        return super().convertToXcart()
    
    def convertToXred(self):
        """
        Converts and returns reduced coordinates of the unit cell using parent class functionality.

        Returns:
            np.ndarray: Reduced coordinates of the unit cell.
        """
        return super().convertToXred()
    
    def change_coordinates(self, new_coordinates, cartesian=False, reduced=False):
        """
        Updates the coordinates of the unit cell to new values and resets the energy attribute.

        Args:
            new_coordinates (np.ndarray): New array of coordinates to set for the unit cell.
            cartesian (bool): If True, indicates the new coordinates are in Cartesian form. Default is False.
            reduced (bool): If True, indicates the new coordinates are in reduced form. Default is False.
        """
        return super().change_coordinates(new_coordinates, cartesian, reduced)


    #-----------------------------
    # File Writing Methods
    #-----------------------------

    def write_batch_script(self, batch_script_header_file='default_batch_file', input_file='input.in', batch_name='default_output', host_spec=None, log='log'):
        """
        Writes a batch script based on a predefined header, customizing it for Abinit execution.

        Args:
            batch_script_header_file (str): Path to the batch script header file. Defaults to 'default_batch_file'.
            input_file (str): Input file name used in Abinit execution. Defaults to 'input.in'.
            batch_name (str): Name of the resulting batch script. Defaults to 'default_output'.
            host_spec (str): Host specification line for distributed computing. Optional.
            log (str): Log file name for capturing output messages. Defaults to 'log'.
        
        Returns:
            bool: True if the batch script was successfully written; False otherwise.
        """

        # Read the content from the batch_script_header file
        try:
            with open(batch_script_header_file, 'r') as header_file:
                batch_script_header = header_file.read()
        except FileNotFoundError:
            print(f"Error: The file {batch_script_header_file} does not exist.")
            return False
        except Exception as e:
            print(f"An error occurred while reading the batch script header file: {e}")
            return False

        # Write to the output file
        try:
            with open(f"{batch_name}", 'w') as file:
                # Write the contents of the batch script header
                file.write("#!/bin/bash\n")
                file.write(batch_script_header)
                
                if host_spec is None:
                    file.write(f"\nmpirun -np 8 abinit < {input_file} > {log} \n")
                else:
                    file.write(f"\n{host_spec} abinit < {input_file} > {log} \n")

            print("Batch script was written successfully.")
            return True

        except Exception as e:
            print(f"An error occurred while writing the batch script: {e}")
            return False
        
    def write_ground_workfunction_file(self, output_path):
        """
        Creates an Abinit input file to calculate the ground work function of the unit cell.

        Args:
            output_path (str): Directory to save the generated Abinit input file.
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

        self.write_custom_abifile(output_file=output_path, header_file=content)


    def write_phonon_dispersion_file(self, output_path):
        """
        Generates an Abinit input file for calculating the phonon dispersion curve of the unit cell.

        Args:
            output_path (str): Directory to save the generated Abinit input file.
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
        self.write_custom_abifile(output_file=output_path, header_file=content)

    def write_custom_abifile(self, output_file, header_file, toldfe=True):
        """
        Writes a custom Abinit .abi file using user-defined or default parameters.

        Args:
            output_file (str): Path where the new Abinit file will be saved.
            header_file (str): Header content or path to a header file.
            toldfe (bool): Flag indicating whether to append toldfe parameter. Defaults to True.
        """
        # Determine if the header_file is actual content or a path to a file
        if "\n" in header_file or not os.path.exists(header_file):
            # If it's likely content due to newline characters or non-existent path
            header_content = header_file
        else:
            # It's a valid file path; read the content from the file
            with open(header_file, 'r') as hf:
                header_content = hf.read()

        # Write all content to the output file
        with open(f"{output_file}.abi", 'w') as outf:
            outf.write(header_content)
            
            # Append unit cell details
            outf.write("\n# Definition of unit cell")
            outf.write(f"\n#*********************************\n")
            outf.write(f"acell {' '.join(map(str, self.acell))}\n")
            outf.write(f"rprim\n")
            for coord in self.rprim:
                outf.write(f"  {'  '.join(map(str, coord))}\n")

            if self.coord_type == 'reduced':
                outf.write(f"xred\n")
            else:
                outf.write(f"xcart\n")

            for coord in self.coordinates:
                # Convert each numpy array to a flat list
                outf.write(f"  {'  '.join(map(str, coord))}\n")

            outf.write("\n# Definition of atoms")
            outf.write(f"\n#*********************************\n")
            outf.write(f"natom {self.num_atoms} \n")
            outf.write(f"ntypat {self.ntypat} \n")
            outf.write(f"znucl {' '.join(map(str, self.znucl))}\n")
            outf.write(f"typat {' '.join(map(str, self.typat))}\n")

            if self.convergence_path is None:
                outf.write("\n# Definition of the planewave basis set")
                outf.write(f"\n#*********************************\n")
                outf.write(f"ecut {self.ecut} \n")
                if self.ecutsm is not None:
                    outf.write(f"ecutsm {self.ecutsm} \n")

                outf.write(f"\n# Definition of the k-point grid")
                outf.write(f"\n#********************************* \n")
                outf.write(f"nshiftk {self.nshiftk} \n")
                if self.kptrlatt is not None:
                    for i in self.kptrlatt:
                        outf.write(f"  {' '.join(map(str, i))}\n")
                outf.write(f"shiftk {' '.join(map(str, self.shiftk))} \n")
                outf.write(f"nband {self.nband} \n")
                outf.write("\n# Definition of the SCF Procedure")
                outf.write(f"\n#********************************* \n")
                outf.write(f"nstep {self.nstep} \n")
                outf.write(f"diemac {self.diemac} \n")
                outf.write(f"ixc {self.ixc} \n")

                if toldfe == True:
                    outf.write(f"toldfe {self.toldfe}\n")

                outf.write(f"\npp_dirpath \"{self.pp_dirpath}\" \n")
                outf.write(f"pseudos \"{self.pseudos}\" \n")
                print(f"{output_file} was created successfully!")
            else:
                with open(self.convergence_path, 'r') as cf:
                    convergence_content = cf.read()
                outf.write(convergence_content)


        
    
    def create_mrgddb_file(self, content):
        pass

    def create_anaddb_file(self, content):
        pass

    #----------------------
    # File Execution Methods
    #----------------------

    def run_abinit(self, input_file='abinit', batch_name='abinit_job',
                batch_script_header_file=None, host_spec='mpirun -hosts=localhost -np 30', 
                delete_batch_script=True, log="log"):
        """
        Executes the Abinit program using a generated input file and specified settings.

        Args:
            input_file (str): The name of the Abinit input file to be executed. Defaults to 'abinit'.
            batch_name (str): The name of the batch script to be generated. Defaults to 'abinit_job'.
            batch_script_header_file (str): File path for the batch script header.
            host_spec (str): Command specifying options for parallel execution environment. Defaults as shown.
            delete_batch_script (bool): Whether to delete the batch script after execution. Default is True.
            log (str): Filename for logging output. Defaults to "log".
        """

        
        # Compile the content to be written in the file
        content = f"""{input_file}.abi
{input_file}.abo
{input_file}_generic_input_files
{input_file}_generic_output_files
{input_file}_generic_temp_files
    """

        if batch_script_header_file is not None:
            # Create a non-temporary file in the current directory
            file_path = f"{input_file}_abinit_input_data.txt"
            
            with open(file_path, 'w') as file:
                file.write(content)
            try:
                # Use the regular file's path in your operations
                script_created = self.write_batch_script(batch_script_header_file=batch_script_header_file, 
                                                        input_file=file_path, 
                                                        batch_name=f"{batch_name}.sh",
                                                        host_spec=host_spec, 
                                                        log=log)
                print(f"Was the batch script successfully created: {script_created}")

                # Submit the job using subprocess to capture output
                result = subprocess.run(['sbatch', f"{batch_name}.sh"], capture_output=True, text=True)

                if result.returncode == 0:
                    print("Batch job submitted using 'sbatch'.")
                    try:
                        job_number = int(result.stdout.strip().split()[-1])
                        self.runningJobs.append(job_number)
                        print(f"Job number {job_number} added to running jobs.")
                    except (ValueError, IndexError) as e:
                        print(f"Failed to parse job number: {e}")
                else:
                    print("Failed to submit batch job:", result.stderr)

            finally:
                # TODO: Write a method that cleans files for you
                print("Attempt to delete")

        else:
            command = f"abinit < {input_file} > {log}"
            os.system(command)
            print(f"Abinit executed directly. Output written to '{log}'.")

    def run_piezo_calculation(self, host_spec='mpirun -hosts=localhost -np 30'):
        """
        Runs a piezoelectricity calculation for the unit cell using default or provided host specifications.

        Args:
            host_spec (str): Specification for parallel execution environment. Defaults as shown.
        """
        content = f"""ndtset 4

# Set 1: Ground State Self-Consistency
#*************************************

getwfk1 0
kptopt1 1
tolvrs1 1.0d-18

# Set 2: Reponse function calculation of d/dk wave function
#**********************************************************

iscf2 -3
rfelfd2 2
tolwfr2 1.0d-20

# Set 3: Response function calculation of d2/dkdk wavefunction
#*************************************************************

getddk3 2
iscf3 -3
rf2_dkdk3 3
tolwfr3 1.0d-16
rf2_pert1_dir3 1 1 1
rf2_pert2_dir3 1 1 1

# Set 4: Response function calculation to q=0 phonons, electric field and strain
#*******************************************************************************
getddk4 2
rfelfd4 3
rfphon4 1
rfstrs4 3
rfstrs_ref4 1
tolvrs4 1.0d-8

# turn off various file outputs
prtpot 0
prteig 0
"""
        # Get the current working directory
        working_directory = os.getcwd()
        
        # Construct the full paths for the output and batch files
        output_file = os.path.join(working_directory, f"{self.file_name}_energy")
        batch_name = os.path.join(working_directory, f"{self.file_name}_bscript")
        
        # Use these paths in your methods
        self.write_custom_abifile(output_file=output_file, header_file=content, toldfe=False)
        self.run_abinit(input_file=output_file, batch_name=batch_name, host_spec=host_spec, batch_script_header_file=self.batchScriptHeader_path)

    def run_flexo_calculation(self, host_spec='mpirun -hosts=localhost -np 30'):
        """
        Runs a flexoelectricity calculation for the unit cell using default or provided host specifications.

        Args:
            host_spec (str): Specification for parallel execution environment. Defaults as shown.
        """
        content = f"""ndtset 5

# Set 1: Ground State Self-Consistency
#*************************************

getwfk1 0
kptopt1 1
tolvrs1 1.0d-18

# Set 2: Reponse function calculation of d/dk wave function
#**********************************************************

iscf2 -3
rfelfd2 2
tolwfr2 1.0d-20

# Set 3: Response function calculation of d2/dkdk wavefunction
#*************************************************************

getddk3 2
iscf3 -3
rf2_dkdk3 3
tolwfr3 1.0d-16
rf2_pert1_dir3 1 1 1
rf2_pert2_dir3 1 1 1

# Set 4: Response function calculation to q=0 phonons, electric field and strain
#*******************************************************************************
getddk4 2
rfelfd4 3
rfphon4 1
rfstrs4 3
rfstrs_ref4 1
tolvrs4 1.0d-8
prepalw4 1

# Set 5: Long-wave Calculations
#******************************

optdriver5 10
get1wf5 4
get1den5 4
getddk5 2
getdkdk5 3
lw_flexo5 1

# turn off various file outputs
prtpot 0
prteig 0
"""
        # Get the current working directory
        working_directory = os.getcwd()
        
        # Construct the full paths for the output and batch files
        output_file = os.path.join(working_directory, f"{self.file_name}_energy")
        batch_name = os.path.join(working_directory, f"{self.file_name}_bscript")
        
        # Use these paths in your methods
        self.write_custom_abifile(output_file=output_file, header_file=content, toldfe=False)
        self.run_abinit(input_file=output_file, batch_name=batch_name, host_spec=host_spec, batch_script_header_file=self.batchScriptHeader_path)

    def run_energy_calculation(self, host_spec='mpirun -hosts=localhost -np 30'):
        """
        Runs an energy calculation for the unit cell using default or provided host specifications.

        Args:
            host_spec (str): Specification for parallel execution environment. Defaults as shown.
        """
        content = f"""ndtset 1

# Ground State Self-Consistency
#*******************************

getwfk1 0
kptopt1 1
tolvrs 1.0d-18

# turn off various file outputs
prtpot 0
prteig 0


getwfk 1
useylm 1  # Use of spherical harmonics
kptopt 2  # Takes into account time-reversal symmetry. 

"""
        # Get the current working directory
        working_directory = os.getcwd()
        
        # Construct the full paths for the output and batch files
        output_file = os.path.join(working_directory, f"{self.file_name}_energy")
        batch_name = os.path.join(working_directory, f"{self.file_name}_bscript")
        
        # Use these paths in your methods
        self.write_custom_abifile(output_file=output_file, header_file=content, toldfe=False)
        self.run_abinit(input_file=output_file, batch_name=batch_name, host_spec=host_spec, batch_script_header_file=self.batchScriptHeader_path)


    def run_anaddb_file(self, content):
        pass

    def run_mrgddb_file(self, content):
        pass
        
    #---------------------------------
    # File Extraction methods
    #---------------------------------

    def grab_energy(self, abo_file=None):
        """
        Retrieves and assigns the total energy from a specified Abinit output file (`abo_file`).

        Args:
            abo_file (str, optional): The path to the Abinit output file. Defaults to auto-generated name.
        
        Raises:
            FileNotFoundError: If the specified `abo_file` does not exist.
        """
        if abo_file is None:
            abo_file = f"{self.file_name}_energy.abo"

        # Ensure total_energy_value is initialized
        total_energy_value = None
        
        try:
            with open(abo_file) as f:
                # Read all content as a single string
                abo_content = f.read()

            # Apply the regex pattern to the full content
            match = re.search(r"total_energy\s*:\s*(-?\d+\.\d+E?[+-]?\d*)", abo_content)
            
            if match:
                total_energy_value = match.group(1)
            else:
                print("Total energy not found.")
                
        except FileNotFoundError:
            print(f"The file {abo_file} was not found.")
        
        self.energy = float(total_energy_value)

    def grab_flexo_tensor(self, anaddb_file=None):
        """
        Retrieves and assigns the electronic flexoelectric coefficients from a specified Abinit output file (`abo_file`).

        Args:
            abo_file (str, optional): The path to the Abinit output file. Defaults to auto-generated name.
        
        Raises:
            FileNotFoundError: If the specified `abo_file` does not exist.
        """
        if anaddb_file is None:
            anaddb_file = f"{self.file_name}_energy.abo"
        
        try:
            with open(anaddb_file) as f:
                # Read all content as a single string
                abo_content = f.read()

            # Apply regex to find the flexoelectric tensor
            tensor_pattern = r"Type-II electronic \(clamped ion\) flexoelectric tensor.*?^([\s\S]+?)(?:^\s*$|#)"
            match = re.search(tensor_pattern, abo_content, re.MULTILINE)
            
            if match:
                tensor_str = match.group(1).strip()
                # Convert tensor_str to a numerical array
                self.flexo_tensor = self._parse_tensor(tensor_str)
                print("Flexo tensor successfully extracted.")
            else:
                print("Flexo tensor not found in the file.")
                
        except FileNotFoundError:
            print(f"The file {anaddb_file} was not found.")


    def _parse_tensor(self, tensor_str):
        """
        Parse a string representation of tensor into a nested list or numpy array.

        Args:
            tensor_str (str): String representation of the tensor.

        Returns:
            parsed_tensor (list or ndarray): Parsed tensor as a nested list or numpy array.
        """
        lines = tensor_str.splitlines()
        parsed_tensor = []
        
        for line in lines:
            numbers = [float(num) for num in line.split()]
            parsed_tensor.append(numbers)
        
        return parsed_tensor


    def grab_piezo_tensor(self, anaddb_file=None):
        """
        Retrieves and assigns the total energy from a specified Abinit output file (`abo_file`)
        along with the clamped and relaxed ion piezoelectric tensors.

        Args:
            abo_file (str, optional): The path to the Abinit output file. Defaults to auto-generated name.
            
        Raises:
            FileNotFoundError: If the specified `abo_file` does not exist.
        """
        if anaddb_file is None:
            anaddb_file = f"{self.file_name}_energy.abo"

        # Initialize arrays to store piezoelectric tensors
        piezo_tensor_clamped = None
        piezo_tensor_relaxed = None

        try:
            with open(anaddb_file) as f:
                # Read all content as a single string
                abo_content = f.read()

            # Extract total energy value
            energy_match = re.search(r"total_energy\s*:\s*(-?\d+\.\d+E?[+-]?\d*)", abo_content)
            if energy_match:
                self.energy = float(energy_match.group(1))
            else:
                print("Total energy not found.")

            # Extract clamped ion piezoelectric tensor
            clamped_match = re.search(
                r"Proper piezoelectric constants \(clamped ion\) \(unit:c/m\^2\)\s*\n((?:\s*-?\d+\.\d+\s+\n?)+)", 
                abo_content
            )
            if clamped_match:
                clamped_strings = clamped_match.group(1).strip().split('\n')
                piezo_tensor_clamped = np.array([list(map(float, line.split())) for line in clamped_strings])

            # Extract relaxed ion piezoelectric tensor
            relaxed_match = re.search(
                r"Proper piezoelectric constants \(relaxed ion\) \(unit:c/m\^2\)\s*\n((?:\s*-?\d+\.\d+\s+\n?)+)", 
                abo_content
            )
            if relaxed_match:
                relaxed_strings = relaxed_match.group(1).strip().split('\n')
                piezo_tensor_relaxed = np.array([list(map(float, line.split())) for line in relaxed_strings])
            
        except FileNotFoundError:
            print(f"The file {anaddb_file} was not found.")
        
        self.piezo_tensor_clamped = piezo_tensor_clamped
        self.piezo_tensor_relaxed = piezo_tensor_relaxed

    #--------------------------
    # Utilities
    #--------------------------

    def copy_abinit_unit_cell(self):
        """
        Creates a deep copy of the current AbinitUnitCell instance.

        Returns:
            AbinitUnitCell: A new instance that is a deep copy of the current instance.
        """
        # Perform a deep copy to ensure all nested objects are also copied
        copied_cell = copy.deepcopy(self)
        return copied_cell
    
    def all_jobs_finished(self):
        """
        Checks whether all submitted computational jobs have finished.

        Returns:
            bool: True if all jobs have completed; False otherwise.
        """
        for job_id in self.runningJobs:
            # Check the status of the specific job
            result = subprocess.run(['squeue', '-j', str(job_id)], capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"Error checking job {job_id} status:", result.stderr)
                continue
            
            # If the job is found in the queue, it means it's still running or pending
            if str(job_id) in result.stdout:
                return False
        
        # If none of the jobs were found in the queue, they have all finished
        return True
    
    def wait_for_jobs_to_finish(self, check_time=60):
        """
        Waits until all computational jobs are finished, checking their status periodically.

        Args:
            check_time (int): Time interval between checks, in seconds. Defaults to 60.
        """
        print("Waiting for jobs to finish...")
        while not self.all_jobs_finished():
            print(f"Jobs still running. Checking again in {check_time} seconds...")
            time.sleep(check_time)
        print("All jobs have finished.")
        self.runningJobs = []

    def perturbations(self, pert):
        """
        Applies a given perturbation to the unit cell's coordinates and returns a new UnitCell object.

        Args:
            pert (np.ndarray): Array representing the perturbation to be applied to current coordinates.

        Returns:
            UnitCell: A new instance of UnitCell with updated (perturbed) coordinates.
        """

        # Ensure the perturbation has the correct shape
        perturbation = np.array(pert, dtype=float)
        if perturbation.shape != self.coordinates.shape:
            raise ValueError("Perturbation must have the same shape as the coordinates.")

        copy_cell = self.copy_abinit_unit_cell()
        # Calculate new coordinates by adding the perturbation
        new_coordinates = self.coordinates + perturbation

        # Create a new UnitCell object with the updated coordinates
        # Assuming other properties remain the same; adjust as needed
        copy_cell.coordinates = new_coordinates

        return copy_cell
    