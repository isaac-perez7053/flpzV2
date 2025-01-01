import numpy as np
import os
import subprocess
import sys
from pathlib import Path
from scipy.sparse.linalg import cg
import warnings

from . import AbinitUnitCell
np.set_printoptions(precision=10)

import tracemalloc
tracemalloc.start()

class SmodesProcessor(AbinitUnitCell):
    """
    A class that processes symmetry modes (SMODES) to calculate phonon properties
    and analyzes them using Abinit simulations.

    Attributes:
        target_irrep (str): The irreducible representation targeted in the calculations.
        smodes_path (str): The path to the SMODES executable.
        symm_prec (float): Precision for recognizing symmetry operations.
        disp_mag (float): Magnitude of displacements used in calculations.
        host_spec (str): Specifications used for parallel execution.
        irrep (str): Irreducible representation identifier.
        num_sam (int): Number of systematic atomic modes (SAMs).
        sam_atom_label (list): Labels of atoms in SAMs.
        mass_list (list): List of atomic masses.
        dist_mat (np.ndarray): Displacement matrix for SAMs.
        pos_mat_cart (np.ndarray): Cartesian coordinates of atomic positions.
        type_count (list): Count of each atom type.
        isir (bool): Indicates if infrared active modes are present.
        israman (bool): Indicates if Raman active modes are present.
        transmodes (bool): Indicates if translational modes exist.
        crossdot_ispos (bool): If cross product of rprim vectors is positive.
        mass_matrix (np.ndarray): Matrix containing atomic masses.
        springs_constants_matrix (np.ndarray): Spring constants matrix.
        jobs_ran_abo (list): List of job identifiers processed by Abinit.
        dyn_freqs (list): List of dynamic frequencies.
        fc_evals (np.ndarray): Eigenvalues of the force constant matrix.
        phonon_vecs (np.ndarray): Phonon displacement vectors.
        red_mass (np.ndarray): Reduced masses for each mode.

    Public Methods:
        convertToXcart(): Converts and returns Cartesian coordinates of the unit cell.
        convertToXred(): Converts and returns reduced coordinates of the unit cell.
        findSpaceGroup(): Determines and returns the space group of the unit cell.
        write_custom_abifile(output_file, header_file): Writes a custom Abinit .abi file.
        run_smodes(smodes_input): Executes SMODES and processes its output.
        unstable_phonons(): Outputs unstable phonons or indicates their absence.
        symmadapt(): Runs the symmadapt sub-program for symmetry-related modes adaptation.
    """
    def __init__(self, abi_file, smodes_input, target_irrep, symm_prec=1e-5, disp_mag=0.001, smodes_path='../isobyu/smodes', host_spec='mpirun -hosts=localhost -np 30', b_script_header_file=None):
        """
        Initializes a SmodesProcessor with specified input file, SMODES parameters, and Abinit configurations.

        Args:
            abi_file (str): Path to the Abinit input file.
            smodes_input (str): Path to the SMODES input file.
            target_irrep (str): The target irreducible representation.
            symm_prec (float): Symmetry precision for identifying symmetry operations. Defaults to 1e-5.
            disp_mag (float): Displacement magnitude for calculations. Defaults to 0.001.
            smodes_path (str): Executable path for SMODES. Defaults to '../isobyu/smodes'.
            host_spec (str): Specification string for host setup in distributed computing. Defaults to 'mpirun -hosts=localhost -np 30'.
            b_script_header_file (str, optional): Path to batch script header file. Defaults to None.
        """
        super().__init__(abi_file=abi_file, batch_script_header_file=b_script_header_file)
        self.target_irrep = target_irrep
        self.smodes_path = smodes_path
        self.symm_prec = symm_prec
        self.disp_mag = disp_mag
        self.host_spec = host_spec

        # Initialize attributes to be populated by _create_header
        self.irrep = None
        self.num_sam = None
        self.sam_atom_label =[]
        self.mass_list = []
        self.dist_mat = None
        self.pos_mat_cart = None
        self.type_count = None

        self.isir = False
        self.israman = False
        self.transmodes = False
        self.crossdot_ispos = False

        self.mass_matrix = None
        self.springs_constants_matrix = None

        # Populate header information
        self._create_header(smodes_input)
        
        # Attributes for storing calculated data
        self.jobs_ran_abo = []
        self.dyn_freqs = []
        self.fc_evals = None
        self.phonon_vecs = None
        self.red_mass = None

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
    
    def findSpaceGroup(self):
        """
        Determines and returns the space group of the unit cell using parent class functionality.

        Returns:
            str: The identifier or symbol of the space group.
        """
        return super().findSpaceGroup()
    
    def write_custom_abifile(self, output_file, header_file):
        """
        Writes a custom Abinit .abi file using user-defined parameters.

        Args:
            output_file (str): Path where the new Abinit file will be saved.
            header_file (str): Header content or path to a header file.
        """
        return super().write_custom_abifile(output_file, header_file)

    def change_coordinates(self, new_coordinates, cartesian=False, reduced=False):
        """
        Updates the coordinates of the unit cell to new values and resets the energy attribute.

        Args:
            new_coordinates (np.ndarray): New array of coordinates to set for the unit cell.
            cartesian (bool): If True, indicates the new coordinates are in Cartesian form. Default is False.
            reduced (bool): If True, indicates the new coordinates are in reduced form. Default is False.
        """
        return super().change_coordinates(new_coordinates, cartesian, reduced)

    def _create_header(self, smodes_input):
        """
        Extract header information from SMODES input file and store it in class attributes.

        Args:
            smodes_input (str): Path to the SMODES input file.

        Raises:
            FileNotFoundError: If the SMODES executable is not found at the specified path.
        """
        atom_list = [
            ('1', 'H', '1.008'), ('2', 'He', '4.002'), ('3', 'Li', '6.94'), ('4', 'Be', '9.012'),
            ('5', 'B', '10.81'), ('6', 'C', '12.011'), ('7', 'N', '14.007'), ('8', 'O', '15.999'),
            ('9', 'F', '18.998'), ('10', 'Ne', '20.180'), ('11', 'Na', '22.990'), ('12', 'Mg', '24.305'),
            ('13', 'Al', '26.982'), ('14', 'Si', '28.085'), ('15', 'P', '30.974'), ('16', 'S', '32.06'),
            ('17', 'Cl', '35.45'), ('18', 'Ar', '39.948'), ('19', 'K', '39.098'), ('20', 'Ca', '40.078'),
            ('21', 'Sc', '44.956'), ('22', 'Ti', '47.867'), ('23', 'V', '50.942'), ('24', 'Cr', '51.996'),
            ('25', 'Mn', '54.938'), ('26', 'Fe', '55.845'), ('27', 'Co', '58.933'), ('28', 'Ni', '58.693'),
            ('29', 'Cu', '63.546'), ('30', 'Zn', '65.38'), ('31', 'Ga', '69.723'), ('32', 'Ge', '72.630'),
            ('33', 'As', '74.922'), ('34', 'Se', '78.971'), ('35', 'Br', '79.904'), ('36', 'Kr', '83.798'),
            ('37', 'Rb', '85.468'), ('38', 'Sr', '87.62'), ('39', 'Y', '88.906'), ('40', 'Zr', '91.224'),
            ('41', 'Nb', '92.906'), ('42', 'Mo', '95.95'), ('43', 'Tc', '98'), ('44', 'Ru', '101.07'),
            ('45', 'Rh', '102.91'), ('46', 'Pd', '106.42'), ('47', 'Ag', '107.87'), ('48', 'Cd', '112.41'),
            ('49', 'In', '114.82'), ('50', 'Sn', '118.71'), ('51', 'Sb', '121.76'), ('52', 'Te', '127.60'),
            ('53', 'I', '126.90'), ('54', 'Xe', '131.29'), ('55', 'Cs', '132.91'), ('56', 'Ba', '137.33'),
            ('57', 'La', '138.91'), ('58', 'Ce', '140.12'), ('59', 'Pr', '140.91'), ('60', 'Nd', '144.24'),
            ('61', 'Pm', '145'), ('62', 'Sm', '150.36'), ('63', 'Eu', '151.96'), ('64', 'Gd', '157.25'),
            ('65', 'Tb', '158.93'), ('66', 'Dy', '162.50'), ('67', 'Ho', '164.93'), ('68', 'Er', '167.26'),
            ('69', 'Tm', '168.93'), ('70', 'Yb', '173.05'), ('71', 'Lu', '174.97'), ('72', 'Hf', '178.49'),
            ('73', 'Ta', '180.95'), ('74', 'W', '183.84'), ('75', 'Re', '186.21'), ('76', 'Os', '190.23'),
            ('77', 'Ir', '192.22'), ('78', 'Pt', '195.08'), ('79', 'Au', '196.97'), ('80', 'Hg', '200.59'),
            ('81', 'Tl', '204.38'), ('82', 'Pb', '207.2'), ('83', 'Bi', '208.98'), ('84', 'Po', '209'),
            ('85', 'At', '210'), ('86', 'Rn', '222'), ('87', 'Fr', '223'), ('88', 'Ra', '226'),
            ('89', 'Ac', '227'), ('90', 'Th', '232.04'), ('91', 'Pa', '231.04'), ('92', 'U', '238.03'),
            ('93', 'Np', '237'), ('94', 'Pu', '244'), ('95', 'Am', '243'), ('96', 'Cm', '247'),
            ('97', 'Bk', '247'), ('98', 'Cf', '251'), ('99', 'Es', '252'), ('100', 'Fm', '257'),
            ('101', 'Md', '258'), ('102', 'No', '259'), ('103', 'Lr', '266'), ('104', 'Rf', '267'),
            ('105', 'Db', '268'), ('106', 'Sg', '269'), ('107', 'Bh', '270'), ('108', 'Hs', '277'),
            ('109', 'Mt', '278'), ('110', 'Ds', '281'), ('111', 'Rg', '282'), ('112', 'Cn', '285'),
            ('113', 'Nh', '286'), ('114', 'Fl', '289'), ('115', 'Mc', '290'), ('116', 'Lv', '293'),
            ('117', 'Ts', '294'), ('118', 'Og', '294')
        ]


        if not Path(self.smodes_path).is_file():
            raise FileNotFoundError(f"SMODES executable not found at: {self.smodes_path}. Current directory is {os.getcwd()}")

        # Open and read SMODES input file
        with open(smodes_input) as s:
            s_lines = s.readlines()
    

        # Parse lattice parameters
        prec_lat_param = [float(x) for x in s_lines[0].split()]

        print(f"Precision Lattice Parameters:\n {prec_lat_param}\n ")
        # I need to check if this is correct
        self.acell = [1, 1, 1]

        # Execute SMODES and process output
        command = f"{self.smodes_path} < {smodes_input}"
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

        # The computer doesn't like if I don't wait :3 
        proc.wait()
        output = proc.stdout.read().decode('ascii')
        proc.stdout.close()

        # Process the output from SMODES
        start_target = 999
        end_target = 0
        outlist = output.split("\n")

        for line in range(len(outlist)):
            line_content = outlist[line].split()
            if len(line_content) > 1 and line_content[0] == "Irrep" and line_content[1] == self.target_irrep:
                start_target = line

            if len(line_content) > 0 and start_target < 999:
                if line_content[0] == "***********************************************":
                    end_target = line
                    break

        target_output = outlist[start_target:end_target]

        if (target_output[3].split()[0]=='These'):
            self.transmodes = True
            del target_output[3] 

        if (target_output[3].split()[0]=='IR'):
            self.isir = True
            del target_output[3] 

        if (target_output[3].split()[0]=='Raman'):
            self.israman = True
            del target_output[3] 


        # Parse degeneracy and number of modes
        degeneracy = int(target_output[1].split()[-1])
        num_modes_without_degen = int(target_output[2].split()[-1])
        num_modes = num_modes_without_degen // degeneracy

        print(f"Degeneracy: {degeneracy}\n")
        print(f"Number of Modes: {num_modes_without_degen}")
        print(f"(Meaning {num_modes} modes to find) \n")

        # Process lattice vectors and atomic positions
        v1 = [float(i) for i in target_output[4].split()]
        v2 = [float(i) for i in target_output[5].split()]
        v3 = [float(i) for i in target_output[6].split()]
        shape_cell = np.array([v1, v2, v3])

        atom_names = []
        atom_positions_raw = []

        for l in range(8, len(target_output)):
            line_content = target_output[l].split()
            if target_output[l] == "Symmetry modes:":
                break
            if len(line_content) >= 4:
                atom_names.append(line_content[1])
                atom_positions_raw.append([float(line_content[2]), float(line_content[3]), float(line_content[4])])

        # Number of atoms (Begin updating the Abinit file with it's new attributes )
        self.natom = len(atom_names)

        # Find number of atom types and count them
        self.ntypat = len(list(set(atom_names)))

        # Dictionary to count the occurances of each element
        count_dict = {}
        result = []

        # Iterate over each element in the input array
        for item in list(atom_names):
            # Increment the count of the element
            if item in count_dict:
                count_dict[item] += 1
            else:
                count_dict[item] = 1
        
        # Create a list of multiplicities based on the original order
        multiplicity_list = []
        # Set to keep track of already processed elements
        seen = set()

        for item in list(atom_names):
            if item not in seen:
                multiplicity_list.append(count_dict[item])
                seen.add(item)

        self.type_count = multiplicity_list

        for index, value in enumerate(multiplicity_list):
            result.extend([(index + 1)] * value)

        self.typat = result 


        # Clean the cell for possible irrational numbers and clean fractions
        clean_list = self._generate_clean_list()
        shape_cell = self._clean_matrix(shape_cell, clean_list)

        # Extract the cell
        prec_lat_array = np.array([prec_lat_param, prec_lat_param, prec_lat_param])
        self.rprim = np.multiply(shape_cell, prec_lat_array)

        # Clean up atom positions
        atom_positions = self._clean_positions(atom_positions_raw, prec_lat_param, clean_list)
        print(f"Smodes Unit Cell Coordiantes:\n {atom_positions} \n")
        self.change_coordinates(new_coordinates=atom_positions, cartesian=True)
        pos_mat_cart = self.coordinates.copy()
        
        # Return list of atom names without duplicates
        seen = set()
        atom_names_nodup = []

        for item in atom_names:
            if item not in seen:
                atom_names_nodup.append(item)
                seen.add(item)
        self.type_list = atom_names_nodup

        # Find all the masses of the atoms
        atomic_num_list, atomic_mass_list = self._get_atomic_details(atom_names_nodup, atom_list)

        self.znucl = atomic_num_list

        # Store extracted data as class attributes
        self.irrep = self.target_irrep
        self.num_sam = num_modes
        self.mass_list = atomic_mass_list
        self.pos_mat_cart = pos_mat_cart

        # Calculate the displacement matrix
        start_line = self.natom + 11
        dist_mat = self._calculate_displacement_matrix(target_output, num_modes, self.natom, start_line)
        # Normalize and orthogonalize SAMs
        self.dist_mat = self._orthogonalize_sams(dist_mat, num_modes, self.natom)

        # Abinit requires this be positive
        crossDot = np.dot(np.cross(self.rprim[0,:], self.rprim[1,:]),np.transpose(self.rprim[2,:]))
        self.crossdot_ispos = crossDot > 0



    def _generate_clean_list(self):
        """
        Generate a list of rational approximations for cleanup.

        Returns:
            list: List of clean values for matrix approximation.
        """
        clean_list = [1.0 / 3.0, 2.0 / 3.0]
        for i in range(1, 10):
            for base in [np.sqrt(3), np.sqrt(2)]:
                clean_list.extend([
                    base / float(i), 2 * base / float(i), 3 * base / float(i),
                    4 * base / float(i), 5 * base / float(i), float(i) / 6.0, float(i) / 8.0
                ])
        return clean_list

    def _clean_matrix(self, matrix, clean_list):
        """
        Clean a matrix by replacing approximate values with exact ones using clean_list.

        Args:
            matrix (np.ndarray): Input matrix to be cleaned.
            clean_list (list): List of target values for cleaning.

        Returns:
            np.ndarray: Cleaned matrix.
        """
        for n in range(matrix.shape[0]):
            for i in range(matrix.shape[1]):
                for c in clean_list:
                    if abs(abs(matrix[n, i]) - abs(c)) < self.symm_prec:
                        matrix[n, i] = np.sign(matrix[n, i]) * c
        return matrix

    def _clean_positions(self, positions, prec_lat_param, clean_list):
        """
        Clean atomic positions and convert using lattice parameters.

        Args:
            positions (list): List of raw atomic positions.
            prec_lat_param (list): Lattice parameters for conversion.
            clean_list (list): List of values to use for cleaning.

        Returns:
            np.ndarray: Cleaned and converted atomic positions.
        """
        # Copy positions to avoid modifying the original input data
        cleaned_positions = positions.copy()

        for n, pos in enumerate(cleaned_positions):
            for i in range(3):
                for c in clean_list:
                    if abs(abs(pos[i]) - abs(c)) < self.symm_prec:
                        pos[i] = np.sign(pos[i]) * c
                pos[i] *= prec_lat_param[i]

        # Convert to a NumPy array to ensure consistent processing
        cleaned_positions = np.array(cleaned_positions)

        # Ensure dimensions are correct 
        if cleaned_positions.ndim != 2 or cleaned_positions.shape[1] != 3:
            raise ValueError(f"Cleaned positions do not have expected shape (n_atoms, 3): {cleaned_positions.shape}")

        return np.array(cleaned_positions)


    def _get_atomic_details(self, type_list, atom_list):
        """
        Retrieve atomic numbers and masses based on atom type.

        Args:
            type_list (list): List of atom types present in the unit cell.
            atom_list (list): Comprehensive list of known elements.

        Returns:
            tuple: Two lists containing atomic numbers and masses respectively.
        """
        atom_dict = {record[1]: (record[0], float(record[2])) for record in atom_list}

        atomic_num_list = []
        atomic_mass_list = []

        for atom_type in type_list:
            if atom_type in atom_dict:
                atomic_number, atomic_mass = atom_dict[atom_type]
                atomic_num_list.append(atomic_number)
                atomic_mass_list.append(atomic_mass)
            else:
                print(f"Unrecognized atom {atom_type}, manual entry required")
                atomic_num_list.append(0)
                atomic_mass_list.append(0)

        return atomic_num_list, atomic_mass_list

        
    def _calculate_displacement_matrix(self, target_output, num_modes, num_atoms, start_line):
        """
        Calculate the initial displacement matrix from SMODES output.

        Args:
            target_output (list): Parsed output lines from SMODES execution.
            num_modes (int): Number of modes considered in calculations.
            num_atoms (int): Total number of atoms present.
            start_line (int): Line number in output where parsing begins.

        Returns:
            np.ndarray: Calculated displacement matrix.
        """
        dist_mat = np.zeros((num_modes, num_atoms, 3))
        mode_index = -1
        self.sam_atom_label = [None] * num_modes

        for l in range(start_line, len(target_output)):
            line_content = target_output[l].split()
            if target_output[l] == '------------------------------------------':
                mode_index += 1
                
            else:
                atom = int(line_content[0]) - 1
                self.sam_atom_label[mode_index] = line_content[1]
                disp1, disp2, disp3 = map(float, line_content[2:5])
                dist_mat[mode_index, atom, 0] = disp1
                dist_mat[mode_index, atom, 1] = disp2
                dist_mat[mode_index, atom, 2] = disp3
        
        return dist_mat

    def _orthogonalize_sams(self, dist_mat, num_modes, num_atoms):
        """
        Normalize and orthogonalize the systematic atomic modes (SAMs).

        Args:
            dist_mat (np.ndarray): Initial displacement matrix.
            num_modes (int): Number of modes.
            num_atoms (int): Number of atoms.

        Returns:
            np.ndarray: Orthogonalized matrix of SAMs.
        """
        # Normalize the SAMs
        for m in range(0, num_modes):
            norm = np.linalg.norm(dist_mat[m, :, :])
            if norm == 0:
                raise ValueError(f"Zero norm encountered at index {m} during normalization.")
            dist_mat[m, :, :] /= norm

        # Orthogonalize the SAMs using a stable Gram-Schmidt Process
        orth_mat = np.zeros((num_modes, num_atoms, 3))
        for m in range(0, num_modes):
            sam = dist_mat[m, :, :]
            
            for n in range(m):
                proj = np.sum(np.multiply(sam, orth_mat[n, :, :])) * orth_mat[n, :, :]
                sam -= proj
            
            # Re-normalize
            norm = np.linalg.norm(sam)
            if norm > 0:
                orth_mat[m, :, :] = sam / norm
            else:
                # Handle the zero norm case, e.g., assigning a zero matrix or handling it differently
                orth_mat[m, :, :] = np.zeros_like(sam)
                print(f"Warning: Zero norm encountered at index {m} during orthogonalization.")

        return orth_mat

    def _loop_modes(self):
        """
        Creates the displaced cells and runs them through Abinit
        for the _perform_calculations method.
        """
        content = """
getwfk 1
useylm 1
kptopt 2
chkprim 0
"""
        # Run the original unchanged cell
        self.convertToXred()
        original_coords = self.coordinates.copy()
        abi_name = "dist_0"
        self.write_custom_abifile(abi_name, content)
        self.run_abinit(input_file=abi_name, batch_name='dist_0_sbatch', batch_script_header_file=self.batchScriptHeader_path, host_spec=self.host_spec, log='dist_0.log')
        self.jobs_ran_abo.append(f"dist_0.abo")

        # Displace each cell
        for i in range(self.num_sam):
            j = i + 1

            self.convertToXcart()
            # Calculate displacement
            perturbation = np.array(self.coordinates + (1.88973 * self.disp_mag * self.dist_mat[i]))
            self.change_coordinates(new_coordinates=perturbation, cartesian=True)
            self.convertToXred()
            abi_name = f"dist_{j}"
            batch_name = f"dist_{j}_sbatch"

            # Write abifile and run abinit
            self.write_custom_abifile(abi_name, content)
            self.run_abinit(input_file=abi_name, batch_name=batch_name, batch_script_header_file=self.batchScriptHeader_path, host_spec=self.host_spec, log=f"dist_{j}.log")

            # Change coordinates back to their original value
            self.change_coordinates(np.array(original_coords).copy(), reduced=True)
            self.jobs_ran_abo.append(f"dist_{j}.abo")

        # Originally set to 300
        self.wait_for_jobs_to_finish(60)


    def _perform_calculations(self, stabilize=False):
        """
        Calculates the eigen-frequencies associated with a particular representation.

        Raises:
            ValueError: If preconditions like positive cross-product are not met.
        """
        if not self.crossdot_ispos:
            raise ValueError("The cross (R1, R2)*R3 must be positive or Abinit will freak out")

        # Ensure force_mat_raw is float64
        force_mat_raw = np.zeros((self.num_sam + 1, self.natom, 3), dtype=np.float64)
        
        for sam, abo in enumerate(self.jobs_ran_abo):
            with open(abo) as f:
                abo_lines = f.readlines()

            line_start = 0
            atom_ind = 0
            
            for line_num, line in enumerate(abo_lines):
                words = line.split()
                if len(words) >= 1 and words[0] == "cartesian" and words[1] == "forces" and words[2] == "(eV/Angstrom)":
                    line_start = line_num + 1
                    break
                    
            for line_num in range(line_start, line_start + self.natom):
                words = abo_lines[line_num].split()
                force_mat_raw[sam, atom_ind, 0] = float(words[1])
                force_mat_raw[sam, atom_ind, 1] = float(words[2])
                force_mat_raw[sam, atom_ind, 2] = float(words[3])
                atom_ind += 1

        print(f"Printing force_mat_raw:\n \n {force_mat_raw} \n")

        # Create force_list with dtype float64
        force_list = np.zeros((self.num_sam, self.natom, 3), dtype=np.float64)
        
        # Subtract off the forces from the original cell
        for sam in range(self.num_sam):
            for i in range(self.natom):
                for j in range(3):
                    force_list[sam, i, j] = force_mat_raw[sam + 1, i, j] - force_mat_raw[0, i, j]

        print(f"Printing force list: \n \n {force_list} \n")
        
        # Initialize the force matrix with dtype float64
        force_matrix = np.zeros((self.num_sam, self.num_sam), dtype=np.float64)

        # Vectorized computation of the force matrix
        force_matrix = np.tensordot(force_list, self.dist_mat.astype(np.float64), axes=([1, 2], [1, 2]))

        # Print the initial force matrix
        print(f"Initial Force Matrix:\n{force_matrix}\n")
        cond_number = np.linalg.cond(force_matrix)
        print(f"Condition number of the force matrix: {cond_number}")

        if stabilize:
            # Stabilize the force matrix
            force_matrix = self.stabilize_matrix(force_matrix)

            # Print the stabilized force matrix
            print(f"Stabilized Force Matrix:\n{force_matrix}\n")

        # Store the force matrix
        self.springs_constants_matrix = np.array(force_matrix, dtype=np.float64)

        # Construct the mass matrix
        #############################

        # Initialize the mass vector
        mass_vector = np.zeros(self.num_sam, dtype=np.float64)

        # Build the mass vector
        for m in range(self.num_sam):
            this_mass = 0
            for n in range(self.ntypat):
                if self.sam_atom_label[m] == self.type_list[n]:
                    this_mass = self.mass_list[n]
                    mass_vector[m] = this_mass
            if this_mass == 0:
                raise ValueError("Problem with building mass matrix. Quitting...")

        # Print the mass vector for debugging
        print(f"Mass Vector:\n{mass_vector}\n")

        # Compute the square root of the mass vector
        sqrt_mass_vector = np.sqrt(mass_vector)

        # Fill the mass matrix using an outer product
        mass_matrix = np.outer(sqrt_mass_vector, sqrt_mass_vector)

        # Print the initial mass matrix
        print(f"Initial Mass Matrix:\n{mass_matrix}\n")
        cond_number = np.linalg.cond(mass_matrix)
        print(f"Condition number of the mass matrix: {cond_number}")

        if stabilize:
            # Stabilize the mass matrix
            mass_matrix = self.stabilize_matrix(mass_matrix, epsilon=1e-6, alpha=0.2)

        # Store the mass matrix
        self.mass_matrix = np.array(mass_matrix, dtype=np.float64)

        # Construct the fc_mat matrix
        ##############################
        fc_mat = (-force_matrix / self.disp_mag).astype(np.float64)
        fc_mat = (fc_mat + np.transpose(fc_mat)) / 2.0

        # Print the initial fc_mat matrix
        print(f"Initial fc_mat Matrix:\n{fc_mat}\n")
        cond_number = np.linalg.cond(fc_mat)
        print(f"Condition number of the fc_mat matrix: {cond_number}")

        if stabilize:
            # Stabilize the fc_mat matrix
            fc_mat = self.stabilize_matrix(fc_mat)

            # Print the stabilized fc_mat matrix
            print(f"Stabilized fc_mat Matrix:\n{fc_mat}\n")

        fc_evals, _ = np.linalg.eig(fc_mat)

        # Construct the dyn_mat matrix
        ###############################
        dyn_mat = np.divide(fc_mat, mass_matrix)

        # Print the initial dyn_mat
        print(f"Initial dyn_mat Matrix:\n{dyn_mat}\n")
        cond_number = np.linalg.cond(dyn_mat)
        print(f"Condition number of the dynamical matrix: {cond_number}")

        if stabilize:
            # Stabilize the dyn_mat matrix
            dyn_mat = self.stabilize_matrix(dyn_mat, epsilon=1e-9, alpha=0.01)

            # Print the stabilized dyn_mat matrix
            print(f"Stabilized dyn_mat Matrix:\n{dyn_mat}\n")

        dynevals, dynevecs_sam = np.linalg.eig(dyn_mat)

        print(f"DEBUG: Printing dynevecs_sam: \n {dynevecs_sam} \n")
        cond_number = np.linalg.cond(dynevecs_sam)
        print(f"Condition number of dynevecs_sam: {cond_number}")

        eV_to_J = 1.602177E-19
        ang_to_m = 1.0E-10
        AMU_to_kg = 1.66053E-27
        c = 2.9979458E10 

        freq_thz = np.multiply(np.sign(dynevals), np.sqrt(np.absolute(dynevals) * eV_to_J / (ang_to_m ** 2 * AMU_to_kg)) * 1.0E-12)
        fc_eval = np.multiply(np.sign(fc_evals), np.sqrt(np.absolute(fc_evals)))

        print(f"DEBUG: Printing freq_thz: \n {freq_thz} \n ")
        idx_dyn = np.flip(freq_thz.argsort()[::-1])
        print(f"DEBUG: Printing idx_dyn, \n {idx_dyn} \n")
        freq_thz = freq_thz[idx_dyn] / (2 * np.pi)
        dynevecs_sam = dynevecs_sam[:, idx_dyn]

        freq_cm = freq_thz * 1.0E12 / c

        self.dyn_freqs = [[freq_thz[i], freq_cm[i]] for i in range(self.num_sam)]
        self.fc_evals = fc_eval[idx_dyn]

        dynevecs = np.zeros((self.num_sam, self.natom, 3), dtype=np.float64)
        
        for evec in range(self.num_sam):
            real_dynevec = np.zeros((self.natom, 3), dtype=np.float64)
            for s in range(self.num_sam):
                real_dynevec += dynevecs_sam[s, evec] * self.dist_mat[s, :, :]
            dynevecs[evec, :, :] = real_dynevec

        print(f"DEGBUG: Printing Dynevecs: \n {dynevecs} \n")

        mass_col = np.zeros((self.natom, 3), dtype=np.float64)
        atomind = 0
        
        for atype in range(self.ntypat):
            for j in range(self.type_count[atype]):
                mass_col[atomind, :] = np.sqrt(self.mass_list[atype])
                atomind += 1

        phon_disp_eigs = np.zeros((self.num_sam, self.natom, 3), dtype=np.float64)
        redmass_vec = np.zeros((self.natom, 1), dtype=np.float64)

        for mode in range(self.num_sam):
            phon_disp_eigs[mode, :, :] = np.divide(dynevecs[mode, :, :], mass_col)
            mag_squared = np.sum(np.sum(np.multiply(phon_disp_eigs[mode, :, :], phon_disp_eigs[mode, :, :])))
            redmass_vec[mode] = 1.0 / mag_squared
            phon_disp_eigs[mode, :, :] /= np.sqrt(mag_squared)

        self.phonon_vecs = phon_disp_eigs 
        # Assign phonon vectors and reduced mass to object attributes
        self.red_mass = redmass_vec.astype(np.float64)
        
        print(f"DEBUG: Printing reduced mass vector: \n {self.red_mass} \n")

        # Store the phonon eigenvectors
        self.phonon_vecs = phon_disp_eigs.astype(np.float64)

        print("Computation completed. The resulting matrices and vectors are stored in the object's attributes.")

    def _perform_calculations_dfpt(self):
        pass
        
    def _imaginary_frequencies(self):
        """
        Detects whether imaginary frequencies exist in dynamic modes.

        Returns:
            list/bool: List of indices with imaginary frequencies, or False if none.
        """
        # List to store indicies with negative freq_thz values
        negative_indicies = []
        print(f"DEBUG: Printing phonons: \n {self.phonon_vecs} \n") 
        # Iterate over dyn_freqs to check the first element of each sublist
        for index, (_, freq_thz) in enumerate(self.dyn_freqs):
            if freq_thz < -20:
                negative_indicies.append(index)
        
        # Return the list of indicies or False if none are negative
        print(f"DEBUG: Printing unstable un normalized phonons: \n {negative_indicies} \n") 
        return negative_indicies if negative_indicies else False
    
    def stabilize_matrix(self, matrix, threshold=50000, epsilon=1e-12, alpha=0.001):
        """
        Stabilize a matrix by regularizing the diagonal, applying weighted symmetrization,
        and adjusting eigenvalues for numerical stability.

        Parameters:
            matrix (numpy.ndarray): The matrix to stabilize.
            threshold (float): Condition number threshold for stabilization.
            epsilon (float): Minimal regularization term for diagonal adjustment.
            alpha (float): Weight for symmetrization to preserve original values.

        Returns:
            numpy.ndarray: The stabilized matrix.
        """
        # Compute the initial condition number
        initial_cond_number = np.linalg.cond(matrix)
        print(f"Initial Condition Number of Matrix: {initial_cond_number}\n")

        if initial_cond_number > threshold:
            print("Condition number is too high; applying stabilization.")

            # Preserve diagonal values while improving stability
            initial_diagonal = np.diag(matrix).copy()

            # Regularize the diagonal minimally
            for i in range(matrix.shape[0]):
                row_sum = np.sum(np.abs(matrix[i, :])) - matrix[i, i]
                if matrix[i, i] < row_sum:
                    matrix[i, i] = (1 - epsilon) * initial_diagonal[i] + epsilon * row_sum

            # Weighted symmetrization to balance original values and numerical stability
            symmetrized_matrix = (matrix + matrix.T) / 2
            matrix = (1 - alpha) * matrix + alpha * symmetrized_matrix

        # Compute the stabilized condition number
        stabilized_cond_number = np.linalg.cond(matrix)
        print(f"Stabilized Condition Number of Matrix: {stabilized_cond_number}\n")

        return matrix

    def run_smodes(self, smodes_input):
        """
        Run the SMODES executable and process its output.

        Args:
            smodes_input (str): Path to SMODES input file.

        Returns:
            str: Output captured from SMODES execution.

        Raises:
            FileNotFoundError: If SMODES executable is not found.
            RuntimeError: If SMODES execution fails.
        """
        if not Path(self.smodes_path).is_file():
            raise FileNotFoundError(f"SMODES executable not found at: {self.smodes_path}")

        # Redirect output to the designated output directory
        command = f"{self.smodes_path} < {smodes_input} > output.log"
        process = subprocess.run(command, shell=True, capture_output=True, text=True)

        if process.returncode != 0:
            raise RuntimeError(f"SMODES execution failed: {process.stderr}")

        return process.stdout

    # TODO: check the vectors are normalized correctly
    def unstable_phonons(self):
        """
        Outputs the unstable phonons and false if none are present.

        Returns:
            list/bool: List of normalized matrices representing unstable phonons, 
                       or False if none are detected.
        """
        unstable_phonons_normalized = []
        if self._imaginary_frequencies() == False:
            print("No unstable phonons are present")
            return False
        else:
            for i in self._imaginary_frequencies():
                # Normalize the matrix as if it were a single vector.
                flattened = self.phonon_vecs[i].flatten()
                euclidean_norm = np.linalg.norm(flattened)
                normalized_flattened = flattened / euclidean_norm
                normalized_matrix = normalized_flattened.reshape(self.phonon_vecs[i].shape)
                unstable_phonons_normalized.append(normalized_matrix)
            print(f"DEBUG: Printing normalized unstable phonons: \n {unstable_phonons_normalized}")
            return unstable_phonons_normalized

    def symmadapt(self):
        """
        Runs the symmadapt sub-program to determine and adapt symmetry-related modes.

        Returns:
            list/bool: Result of unstable phonons calculation.
        """
        self._loop_modes() 
        self._perform_calculations()
        return self.unstable_phonons()      

def main():
    input_file = str(sys.argv[1])
    smodesInput = str(sys.argv[2])
    target_irrep = str(sys.argv[3])  # Example usage with command-line argument
    calculator = SmodesProcessor(input_file, smodesInput, target_irrep)
    
    # print("Dyn Frequencies (THz, cm^-1):", calculator.dyn_freqs)
    # print("Force Constant Evals:", calculator.fc_evals)
    # print("Phonon Vecs Shape:", calculator.phonon_vecs.shape)
    # print("Reduced Masses:", calculator.red_mass)


# Usage example
if __name__ == "__main__":
    main()
