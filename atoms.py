#Use to collect, visualise atomic structures
from ase.io import read

#Load atoms into graphing function
def load_atoms(path, repeat=(1,1,1)):
    return read(path).repeat(repeat)