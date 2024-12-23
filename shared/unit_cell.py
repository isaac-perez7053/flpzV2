import numpy as np
import os 
from pathlib import Path
import sys
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from abinit_file import AbinitFile
from abinit_unit_cell import AbinitUnitCell

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
    # TODO: extract rprim 
    def __init__(self, acell, rprim, coordinates, coord_type, num_atoms, atom_types, znucl, typat, header_path): 
        self.acell = np.array(acell, dtype=float)
        self.rprim = np.array(rprim, dtype=float)
        self.coordiantes = np.array(coordinates, dtype=float)
        self.coord_type = coord_type
        self.num_atoms = num_atoms
        self.atom_types = np.array(atom_types, dtype=int)
        self.znucl = np.array(znucl, dtype=int)
        self.typat = np.array(typat, dtype=int)
        self.header_path = header_path

    def convertToXred(self):
        if self.coord_type == 'reduced': 
            print("The unit cell is already expressed in reduced coordiantes")
        else: 
            self.coord_type = 'reduced'
            # Calculate the reduced coordinates of the system
            rprim_scaled = self.rprim.copy()
            xred = np.empty((self.num_atoms, 3))
            rprim_scaled[0] *= self.acell[0]
            rprim_scaled[1] *= self.acell[1]
            rprim_scaled[2] *= self.acell[2]
            xred = np.dot(self.coordiantes, np.linalg.inv(rprim_scaled))
            self.coordiantes = xred
    
    def convertToXcart(self):
        if self.coord_type == 'cartesian':
            print("The unit cell is already expressed in cartesian coordinates")
        else: 
            self.coord_type = 'cartesian'
            # Caclulate the cartesian coordinates of the system
            rprim_scaled = self.rprim.copy()
            xcart = np.empty((self.num_atoms, 3))
            rprim_scaled[0] *= self.acell[0]
            rprim_scaled[1] *= self.acell[1]
            rprim_scaled[2] *= self.acell[2]
            xcart = np.dot(self.coordiantes, np.linalg.inv(rprim_scaled))
            self.coordiantes = xcart


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
            coords=self.coordiantes, 
            coords_are_cartesian=(self.coord_type=="cartesian")

        )

        analyzer = SpacegroupAnalyzer(structure)
        space_group = analyzer.get_space_group_symbol()
        space_group_number = analyzer.get_space_group_number()

        return f"{space_group} ({space_group_number})"
        

    def createAbinitUnitCell(self):
        """
        Initialize an Abinit Unit Cell object 
        """
        AbinitUnitCell(self)


    def __repr__(self):
        return (f"UnitCell(acell={self.acell}, coord_type='{self.coord_type}', coordinates={self.coordinates}, "
                f"num_atoms={self.num_atoms}, atom_types={self.atom_types}, znucl={self.znucl}, typat={self.typat}, "
                f"header_path='{self.header_path}')")


# Main
if __name__ == "__main__": 
    filepath = input("Enter the path to the Abinit file: ").strip()
    if not os.path.isfile(filepath):
        print(f"File is not found: {filepath}")
    else: 
        molecule = AbinitFile(filepath)
        print(AbinitFile)

            
# An example of how this might be called is, 
# Assuming you have an instance of UnitCell
# unit_cell = molecule.unit_cell  # Example: retrieved from the Molecule class
# space_group = unit_cell.findSpaceGroup()
# print(f"The space group of the unit cell is: {space_group}")

