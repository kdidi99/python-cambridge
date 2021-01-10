import numpy as np
import sys 
import csv

#different methods for getting Hückel matrix (platonic solids are hard-coded)
def hueckel_linear(n):
    M = np.zeros((n,n))
    for i in range(n-1):
        M[i,i+1]=1
        M[i+1,i]=1
    
    return M

def hueckel_cyclic(n):
     M = np.zeros((n,n))
     for i in range(n-1):
         M[i,i+1]=1
         M[i+1,i]=1
         M[0,n-1]=1
         M[n-1,0]=1
    
     return M
 
def hueckel_tetrahedron():
    M = np.matrix([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]])
    
    return M

def hueckel_cube():
    M = np.matrix([[0,1,1,0,0,0,1,0],
                   [1,0,0,1,0,0,0,1],
                   [1,0,0,1,1,0,0,0],
                   [0,1,1,0,0,1,0,0],
                   [0,0,1,0,0,1,1,0],
                   [0,0,0,1,1,0,0,1],
                   [1,0,0,0,1,0,0,1],
                   [0,1,0,0,0,1,1,0],])
    
    return M

def hueckel_dodecahedron():
    #to save time and work, here a publicly available adjacency matrix for dodecahedrons was used as csv
    #downloaded from: https://www.distanceregular.org/graphs/dodecahedron.html
    reader = csv.reader(open("dodecahedron.csv"), delimiter=",")
    x = list(reader)
    result = np.array(x).astype("float")
    M = np.matrix(result)
    
    return M

#evaluation of eigenvalues of Hückel matrix
def get_evals(M):
    evals, evecs = np.linalg.eigh(M)
    eigvals = evals
    eigvals = np.array(eigvals)
    return sorted(eigvals)

#calculate the degeneracy of the eigenvalues
def degeneracy(eigvals):
    #round eigvals so that they can be compared
    eigvals_rounded = np.round(eigvals,2)
    eigvals_rounded = eigvals_rounded.tolist()
    #compare evals and identify degeneracies, 
    #save unique energies in energies and degeneracy in degen
    degen = []
    energies = [] 
    for i in eigvals_rounded: 
        if i not in energies: 
            energies.append(i) 
    l = len(energies)
    for i in range(0,l):
        degen.append(eigvals_rounded.count(energies[i]))       
        
    return energies, degen


#main method that defines the interaction with the user
def main():
    system = input("Please enter the type of system you want to calculate (linear, cyclic, tetrahedron, cube, dodecahedron):\n")

    if system == "linear":
        n1 = input("Please enter the number of carbon atoms your system has:\n")
        n = int(n1)
        M = hueckel_linear(n)
    elif system == "cyclic":
        n1 = input("Please enter the number of carbon atoms your system has:\n")
        n = int(n1)
        M = hueckel_cyclic(n)
    elif system == "tetrahedron":
        M = hueckel_tetrahedron()
    elif system == "cube":
        M = hueckel_cube()
    elif system == "dodecahedron":
        M = hueckel_dodecahedron()
    else:
        sys.exit("The system you provided cannot be calculated by this program. Check for typos.")
    
    evals = get_evals(M)
    energies, degen = degeneracy(evals)
    #pack corresponding energies and degeneracies in tuple and print them
    for en, de in zip(energies, degen): 
        print ("Energy:", en, " Degeneracy:", de)



