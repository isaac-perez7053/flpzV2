import shutil
import re
import os
from unit_cell import UnitCell

# TODO: I may want to place this class into its own file
class AbinitFile: 
    """
    Class that primarily deals with the extraction of the molecule from the abinit file and initializes the UnitCell

    Public Methods: 

    """
    def __init__(self, filepath):
        self.filepath = filepath 
        self.unit_cell = None
        self.parse_file()

    """
    Extract all lines from the Abinit file
    """
    def parse_file(self): 
        temp_filepath = self.filepath + ".temp"
        shutil.copy(self.filepath, temp_filepath)

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

        if num_atoms is None:
            raise Exception("natom is missing in the Abinit file!")

        # Extract ntypat
        atom_types = None
        for i, line in enumerate(lines):
            if line.strip().startswith('ntypat'):
                match = re.search(r"\d+", line)
                if match:
                    atom_types = int(match.group())
                del lines[i]
                break

        if atom_types is None:
            raise Exception("ntypat is missing in the Abinit file!")

        # Extract znucl
        znucl = []
        for i, line in enumerate(lines):
            if line.strip().startswith('znucl'):
                znucl = list(map(int, re.findall(r"\d+", line)))  # Fixed typo in re.finall to re.findall
                del lines[i]
                break

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
        else:
            raise Exception("typat is missing from the Abinit file!")

        if not typat:  # Check if typat is still an empty list
            raise Exception("typat is missing in the Abinit file!")

        # Save the remaining content as the header
        header_path = self.filepath + ".header"
        with open(header_path, "w") as header_file:
            header_file.writelines(lines)
        
        # Create UnitCell instance 
        self.unit_cell = UnitCell(acell, rprim, coordinates, coord_type, num_atoms, atom_types, znucl, typat)

        # Remove temporary file
        os.remove(temp_filepath)
    
    def __repr__(self):
        return f"Molecule(unit_cell={self.unit_cell})"
    

if __name__ == "__main__":
    # For standalone execution and debugging
    abinit_file = AbinitFile(file_path = "example_file.abi")
    print(abinit_file.unit_cell)
    