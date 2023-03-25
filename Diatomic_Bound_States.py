import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs, eigsh
from numpy.linalg import eig

def Hydrogen_Atom_Bound_States(l):

    def Coulomb(grid, l, z=1):
        """ This is a simple coulomb potential with nuclear charge z and l specifying the quantum number 
        which is related to the centrifugal term in the potential"""
        return -z*np.power(grid, -1.0) + 0.5*l*(l+1)*np.power(grid, -2.0)

    potential = Coulomb(grid, l, z=1)

    matrix_size = grid.size 

    total_matrix = []
    matrix_row = []

    diag = np.array(range(matrix_size))
    row = diag
    col = diag
    data = np.zeros(len(row)) + potential + (15.0/ 12.0)/h2


    col = np.concatenate((col, diag[1:]))
    row = np.concatenate((row, diag[:-1]))
    data_right_one =  np.zeros(len(diag[:-1]))  + (-2.0/3.0)/h2
    data = np.concatenate((data, data_right_one))


    col = np.concatenate((col, diag[2:]))
    row = np.concatenate((row, diag[:-2]))
    data_right_two =  np.zeros(len(diag[:-2]))  + (1.0/24.0)/h2
    data = np.concatenate((data, data_right_two))

    col = np.concatenate((col, diag[:-1]))
    row = np.concatenate((row, diag[1:]))
    data_left_one =  np.zeros(len(diag[:-1]))  + (-2.0/3.0)/h2
    data = np.concatenate((data, data_left_one))


    col = np.concatenate((col, diag[:-2]))
    row = np.concatenate((row, diag[2:]))
    data_left_two =  np.zeros(len(diag[:-2]))  + (1.0/24.0)/h2
    data = np.concatenate((data, data_left_two))
    

    sparseMatrix = csr_matrix((data, (row, col)), 
                          shape = (matrix_size, matrix_size)).toarray()

    
    sparseMatrix[0,0] = potential[0]  + (20.0/24.0)/h2
    sparseMatrix[0,1] =  (-6.0/24.0)/h2
    sparseMatrix[0,1] =  (-4.0/24.0)/h2
    sparseMatrix[0,1] =  (1.0/24.0)/h2

    j = grid.size - 1
    sparseMatrix[j,j] = potential[j]  + (20.0/24.0)/h2
    sparseMatrix[0,j-1] =  (-6.0/24.0)/h2
    sparseMatrix[0,j-2] =  (-4.0/24.0)/h2
    sparseMatrix[0,j-3] =  (1.0/24.0)/h2

    # Hamiltonian.setValue(0,0, potential[0]  + (20.0/24.0)/h2)
    # Hamiltonian.setValue(0,1, (-6.0/24.0)/h2)
    # Hamiltonian.setValue(0,2, (-4.0/24.0)/h2)
    # Hamiltonian.setValue(0,3, (1.0/24.0)/h2)
    
    # j = Grid.size - 1
    # Hamiltonian.setValue(j,j, potential[j] + (20.0/24.0)/h2)
    # Hamiltonian.setValue(j,j - 1, (-6.0/24.0)/h2)
    # Hamiltonian.setValue(j,j - 2, (-4.0/24.0)/h2)
    # Hamiltonian.setValue(j,j - 3, (1.0/24.0)/h2)

    vals, vecs = eigs(sparseMatrix, k=3, sigma=-1)

    vals, vecs=eig(sparseMatrix)

    print(vals)

if __name__=="__main__":

    h = 0.0001
    h2 = h*h
    grid_size = 10
    grid = np.arange(h, grid_size, h)

    Hydrogen_Atom_Bound_States(0)