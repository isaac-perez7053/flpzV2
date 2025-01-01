from . import AbinitUnitCell
import numpy as np

class DfptProcessor(AbinitUnitCell):    
    """
    
    
    """
    def __init__(self, abi_file='abi_file.abi', convergence_file='convergence_file.abi'):
        # Initialize the input variables from the Abinit File
        super().__init__(abi_file=abi_file, convergence_path=convergence_file)

    def write_pertrubation_dynamical_file(self):
        pass


