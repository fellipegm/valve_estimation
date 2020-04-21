'''
This package is used to load data saved as .sim into python

'''


import numpy as np
from struct import unpack

class WrongFileHeader(Exception):
    pass

def load_sim(filename):

    with open(filename, "rb") as file:
        file.seek(0,0)
        header = file.read(5)
        expected_header = b'\x73\x69\x6d\x00\x00'
        if header != expected_header:
            print("Wrong file header - the file might be corrupted")
            raise WrongFileHeader
        
        n_variables = file.read(4)
        size_vectors = file.read(4)

        n_variables = int.from_bytes(n_variables, byteorder='little', signed=False)
        size_vectors = int.from_bytes(size_vectors, byteorder='little', signed=False)

        variables = dict()
        for i in range(n_variables):
            name_aux = file.read(20)
            name_aux = name_aux.decode("ascii", errors="ignore").replace("\x00", "")
            variables[name_aux] = []
        
        for var in variables.keys():
            aux_list = []
            for i in range(size_vectors):
                aux_list.extend( unpack('d', file.read(8)) )
            variables[var] = np.array(aux_list, dtype='float64')

    return variables




