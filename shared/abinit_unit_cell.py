from shared.unit_cell_module import UnitCell
from pathlib import Path
import sys
import os
import shutil
import tempfile
import re
import subprocess
import time

# To ensure the file is calling from the correct directory. 
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.append(str(current_dir))

class AbinitUnitCell(UnitCell):
    """
    Contains methods that allow the user to manipulate the unit cell and produce Abinit files

    Public Methods: 

    """
    # TODO: Make a method that will place all pseudopotentials into the pseudopotential folder in the program. All pp_dirpaths will then be the same

    def __init__(self, abi_file, convergence_path=None, batch_script_header_file=None):
        # Call the parent class's initializer with the keyword arguments
        self.abi_file = abi_file
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

   


    def _initialize_convergence_from_file(self):
        """
        Extracts the convergence features from the Abinit file
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
        return super().findSpaceGroup()
    
    def convertToXcart(self):
        return super().convertToXcart()
    
    def convertToXred(self):
        return super().convertToXred()

    def groundWFKFile(self, output_path):
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

        self.write_custom_abifile(output_file=output_path, header_file=content)


    def phononDispFile(self, output_path):
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
        self.write_custom_abifile(output_file=output_path, header_file=content)

    def write_custom_abifile(self, output_file, header_file): 
        """ 
        Creates a custom Abinit file using the attributes of the class or custom files.
        
        Args: 
            output_file (str): Path to save new Abinit file
            header_file (str): The header content or path to a header file
        """
        # Determine if the header_file is actual content or a path to a file
        if "\n" in header_file or not os.path.exists(header_file):
            # If it's likely a content due to newline characters or non-existent path
            header_content = header_file
        else:
            # It's a valid file path; read the content from the file
            with open(header_file, 'r') as f:
                header_content = f.read()

        # Write initial content to the output file
        with open(output_file, 'w') as f: 
            f.write(header_content)

            # Append unit cell details
            f.write("\n# Definition of unit cell")
            f.write(f"\n #*********************************\n")
            f.write(f"acell {' '.join(map(str, self.acell))}\n")
            f.write(f"rprim\n")
            for coord in self.rprim:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            if self.coord_type == 'reduced':
                f.write(f"xred\n")
            else:
                f.write(f"xcart\n")

            for coord in self.coordinates:
                f.write(f"  {'  '.join(map(str, coord))}\n")

            f.write("\n# Definition of atoms")
            f.write(f"\n #*********************************\n")
            f.write(f"natom {self.num_atoms} \n")
            f.write(f"ntypat {self.ntypat} \n")
            f.write(f"znucl {' '.join(map(str, self.znucl))}\n")
            f.write(f"typat {' '.join(map(str, self.typat))}\n")

            if self.convergence_path is None:
                f.write("\n # Definition of the planewave basis set")
                f.write(f"\n #*********************************\n")
                f.write(f"ecut {self.ecut} \n")
                if self.ecutsm is not None:
                    f.write(f"ecutsm {self.ecutsm} \n")

                f.write(f"\n # Definition of the k-point grid")
                f.write(f"\n #********************************* \n")
                f.write(f"nshiftk {self.nshiftk} \n")
                if self.kptrlatt is not None:
                    for i in self.kptrlatt:
                        f.write(f"  {' '.join(map(str, i))}\n")
                f.write(f"shiftk {' '.join(map(str, self.shiftk))} \n")
                f.write(f"nband {self.nband} \n")
                f.write("\n # Definition of the SCF Procedure")
                f.write(f"\n #********************************* \n")
                f.write(f"nstep {self.nstep} \n")
                f.write(f"diemac {self.diemac} \n")
                f.write(f"ixc {self.ixc} \n")
                f.write(f"toldfe {self.toldfe}\n")
                f.write(f"\npp_dirpath \"{self.pp_dirpath}\" \n")
                f.write(f"pseudos \"{self.pseudos}\" \n")
            else: 
                with open(self.convergence_path, 'r') as f:
                    convergence_content = f.read()
                with open(output_file, 'w') as f: 
                    f.write(convergence_content)
                
    def all_jobs_finished(self):
            """
            Checks if all jobs in self.runningJobs have finished.

            Returns:
                bool: True if all jobs are finished, False otherwise.
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
        Waits until all jobs in self.runningJobs are finished.
        """
        print("Waiting for jobs to finish...")
        while not self.all_jobs_finished():
            print(f"Jobs still running. Checking again in {check_time} seconds...")
            time.sleep(check_time)
        print("All jobs have finished.")

    def run_abinit(self, input_file='abinit.in', output_file='abinit_job.sh', batch_script_header_file=None, host=None, delete_batch_script=True, log='log'):
        """
        Run the Abinit program using the generated input file.

        Args:
            input_file (str): The abinit file that will be executed
            output_file (str): Either the name of the created batch
                               script or the name of the log file if Abinit is ran directly
            batch_script_header_file (str): The header of the batch script
            host (str): The command to be run in the batch script
            delete_batch_scrtip (bool): Option to delete batch script after execution
        """

        if batch_script_header_file is not None:
            self.write_batch_script(batch_script_header_file=batch_script_header_file, output_file=output_file, host=host, log=log)
            
            # Submit the job using subprocess to capture output
            result = subprocess.run(['sbatch', output_file], capture_output=True, text=True)
            
            if result.returncode == 0:
                print("Batch job submitted using 'sbatch'.")
                
                # Extract job number (assuming it's in the output)
                # This assumes your system outputs something like: "Submitted batch job 12345"
                try:
                    job_number = int(result.stdout.strip().split()[-1])
                    self.runningJobs.append(job_number)
                    print(f"Job number {job_number} added to running jobs.")
                except (ValueError, IndexError) as e:
                    print(f"Failed to parse job number: {e}")
                    
            else:
                print("Failed to submit batch job:", result.stderr)

            # Delete the batch script after submitting the job
            if os.path.exists(output_file) and delete_batch_script:
                try:
                    os.remove(output_file)
                    print(f"Batch script '{output_file}' has been deleted.")
                except OSError as e:
                    print(f"Error deleting batch script: {e}")

        else:
            command = f"abinit < {input_file} > {output_file}"
            os.system(command)
            print(f"Abinit executed directly. Output written to '{output_file}'.")

    def write_batch_script(self, batch_script_header_file='default_batch_file', output_file='default_output', host=None, log='log'):
        """
        Writes a batch script using a prewritten header file 

        Args: 
            batch_script_header_file (str): The file containing the header of the batch script
            output_file (str): The name of the created batch script
            host (str): The command to be ran in the batch script
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
            with open(f"{output_file}", 'w') as file:
                # Write the contents of the batch script header
                file.write(batch_script_header)
                
                if host is None:
                    file.write(f"\nmpirun -np 8 abinit < temp_input_file > {log} \n")
                else:
                    file.write(f"\n {host} \n")

            print("Batch script was written successfully.")
            return True

        except Exception as e:
            print(f"An error occurred while writing the batch script: {e}")
            return False


