import numpy as np
import os 
from pathlib import Path
import sys
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import re
import shutil

# To ensure the file is calling from the correct directory. 
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.append(str(current_dir))


# TODO: Test this for a bunch of cases
class UnitCell:
    """
    Defines the UnitCell class which contains the all the necessary information of a UnitCell from an abinit file

    Public Methods: 
        find_spaceGroup(): Returns space group of the UnitCell
        xred(): Returns the reduced coordinates of the UnitCell
        xcart(): Returns the cartesian coordinats of the UnitCell
        groundWFKFile(): Generate an Abinit groundWFK file 
        phononDispFile(): Generate an Abinit file used to calculate the phonon dispersion curve

    """
    
    def __init__(self, **kwargs):
        """
        Initialize the class. Either provide the variables directly or specify an abifile for automatic configuration.

        Acceptable keyword arguments are:
        - acell, rprim, coordinates, coord_type, num_atoms, atom_types, znucl, typat, header_path, abi_file
        """

        # Set default values if needed
        self.acell = None
        self.rprim = None
        self.coordinates = None
        self.coord_type = None
        self.num_atoms = None
        self.ntypat = None
        self.znucl = None
        self.typat = None
        self.header_path = kwargs.get('header_path')

        # If 'abi_file' is provided and exists, initialize from file
        abi_file = kwargs.get('abi_file')
        if abi_file and os.path.isfile(abi_file):
            self.abi_file = abi_file
            self._initialize_from_abifile()
        else:
            # Use the provided keyword arguments
            required_fields = ['acell', 'rprim', 'coordinates', 'coord_type', 'num_atoms', 'atom_types', 'znucl', 'typat']
            missing_fields = [field for field in required_fields if field not in kwargs]

            if missing_fields:
                raise ValueError(f"Missing required fields: {', '.join(missing_fields)}")

            self.acell = np.array(kwargs['acell'], dtype=float)
            self.rprim = np.array(kwargs['rprim'], dtype=float)
            self.coordinates = np.array(kwargs['coordinates'], dtype=float)
            self.coord_type = kwargs['coord_type']
            self.num_atoms = kwargs['num_atoms']
            self.atom_types = np.array(kwargs['atom_types'], dtype=int)
            self.znucl = np.array(kwargs['znucl'], dtype=int)
            self.typat = np.array(kwargs['typat'], dtype=int)

    def _initialize_from_abifile(self):
        """
        Extracts all lines from the Abinit file
        """

        temp_filepath = self.abi_file + ".temp"
        shutil.copy(self.abi_file, temp_filepath)

        with open(temp_filepath, 'r') as f: 
            lines = f.readlines()

        # Extract acell
        acell = []
        for i, line in enumerate(lines):
            if line.strip().startswith('acell'):
                # Extract different features of the acell feature and map it to a float where it will be inserted into a list.
                match = re.search(r"(\d+)\*([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*)", line)
                if match:
                    count = int(match.group(1))
                    value = float(match.group(2))
                    acell = [value] * count
                else:
                    acell = list(map(float, re.findall(r"[-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*", line)))
                # Delete extracted lines in copy
                del lines[i]

        self.acell = acell
        if not acell:  # Check if acell is still an empty list
            raise Exception("acell is missing in the Abinit file!")

        # Extract primitive vectors
        rprim = []
        for i, line in enumerate(lines):
            if line.strip() == 'rprim':
                del lines[i]
                j = i
                while j < len(lines) and re.match(r"^\s*[-+]?[0-9]*\.?[0-9]+", lines[j]):
                    rprim.append(list(map(float, lines[j].split())))
                    del lines[j]
                break
        else:
            # Default rprim value if not specified
            rprim = [[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]]
            
        self.rprim = rprim

        # Extract coordinates (xred or cartesian)
        coord_type = None
        coordinates = []
        for i, line in enumerate(lines):
            if line.strip() == 'xred':
                coord_type = 'reduced'
                del lines[i]
                j = i
                while j < len(lines) and re.match(r"^\s*[-+]?[0-9]*\.?[0-9]+", lines[j]):
                    coordinates.append(list(map(float, lines[j].split())))
                    del lines[j]
                break
            elif line.strip() == 'xcart':
                coord_type = 'cartesian'
                del lines[i]
                j = i
                while j < len(lines) and re.match(r"^\s*[-+]?[0-9]*\.?[0-9]+", lines[j]):
                    coordinates.append(list(map(float, lines[j].split())))
                    del lines[j]
                break

        self.coordinates = coordinates
        self.coord_type = coord_type
        if not coordinates:  # Check if coordinates list is empty
            raise Exception("coordinates are missing in the Abinit file!")

        # Extract natom
        num_atoms = None
        for i, line in enumerate(lines):
            if line.strip().startswith('natom'):
                match = re.search(r"\d+", line)
                if match:
                    num_atoms = int(match.group())
                del lines[i]
                break

        self.num_atoms = num_atoms
        if num_atoms is None:
            raise Exception("natom is missing in the Abinit file!")

        # Extract ntypat
        ntypat = None
        for i, line in enumerate(lines):
            if line.strip().startswith('ntypat'):
                match = re.search(r"\d+", line)
                if match:
                    ntypat = int(match.group())
                del lines[i]
                break

        self.ntypat = ntypat
        if self.ntypat is None:
            raise Exception("ntypat is missing in the Abinit file!")

        # Extract znucl
        znucl = []
        for i, line in enumerate(lines):
            if line.strip().startswith('znucl'):
                znucl = list(map(int, re.findall(r"\d+", line)))  # Fixed typo in re.finall to re.findall
                del lines[i]
                break
        
        self.znucl = znucl
        if not znucl:  # Check if znucl is still an empty list
            raise Exception("znucl is missing in the Abinit file!")

        # Extract typat
        typat = []
        for i, line in enumerate(lines):
            if line.strip().startswith('typat'):
                # TODO: Length may be hardcoded
                typat_tokens = re.findall(r"(\d+\*\d+|\d+)", line)
                for token in typat_tokens:
                    if '*' in token:
                        count, value = map(int, token.split('*'))
                        typat.extend([value] * count)
                    else:
                        typat.append(int(token))
                del lines[i]
                break
        
        self.typat = typat
        if not typat:  # Check if typat is still an empty list
            raise Exception("typat is missing in the Abinit file!")

        # Save the remaining content as the header
        header_path = self.abi_file + ".header"
        self.header_path = header_path

        with open(header_path, "w") as header_file:
            header_file.writelines(lines)


    def _close_to_zero(self, value, threshold=1e-10):
        """Helper function to treat small numbers as zero."""
        return np.abs(value) < threshold

    def convertToXred(self):
        if self.coord_type == 'reduced':
            print("The unit cell is already expressed in reduced coordinates")
        else:
            self.coord_type = 'reduced'
            # Calculate the reduced coordinates of the system
            rprim_scaled = self.rprim.copy()
            # Scaling the primitive vectors by the lattice constants
            for i in range(3):
                rprim_scaled[i] *= self.acell[i]
            
            # Convert Cartesian to reduced coordinates
            xred = np.dot(self.coordinates, np.linalg.inv(rprim_scaled))
            
            # Set small values to zero
            xred[np.abs(xred) < 1e-10] = 0
            
            self.coordinates = xred
            print("Converted to reduced coordinates.")

    def convertToXcart(self):
        if self.coord_type == 'cartesian':
            print("The unit cell is already expressed in cartesian coordinates")
        else:
            self.coord_type = 'cartesian'
            # Calculate the cartesian coordinates of the system
            rprim_scaled = self.rprim.copy()
            # Scaling the primitive vectors by the lattice constants
            for i in range(3):
                rprim_scaled[i] *= self.acell[i]
            
            # Convert reduced to Cartesian coordinates
            xcart = np.dot(self.coordinates, rprim_scaled)
            
            # Set small values to zero
            xcart[np.abs(xcart) < 1e-10] = 0
            
            self.coordinates = xcart
            print("Converted to Cartesian coordinates.")
    
    def findSpaceGroup(self):
        """
        Calculates and returns the space group of the unit cell. Ensure coordinates are in cartesian

        Returns: 
            str: The space group symbol (e.g. 'Pm-3m') and international number.
        """

        # Scale the lattice vectors by acell 
        lattice = Lattice(self.rprim * self.acell[:, None])

        # Map typat and znucl to get the atomic numbers
        species = [self.znucl[typ - 1] for typ in self.typat]

        # Create the pymatgen Structure
        structure = Structure(
            lattice=lattice,
            species=species, 
            coords=self.coordinates, 
            coords_are_cartesian=(self.coord_type=="cartesian")

        )

        analyzer = SpacegroupAnalyzer(structure)
        space_group = analyzer.get_space_group_symbol()
        space_group_number = analyzer.get_space_group_number()

        return f"{space_group} ({space_group_number})"

    def __repr__(self):
        return (f"UnitCell(acell={self.acell}, coord_type='{self.coord_type}', rprim={self.rprim}, coordinates={self.coordinates}, "
                f"num_atoms={self.num_atoms}, atom_types={self.atom_types}, znucl={self.znucl}, typat={self.typat}, "
                f"header_path='{self.header_path}')")


