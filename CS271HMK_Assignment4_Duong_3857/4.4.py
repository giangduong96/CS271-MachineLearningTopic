"""
Name: Giang Duong
ID : 014533857
# Homework 4.1 - PCA
"""
import numpy as np
A = [[1,2], [3,1], [2, 3]]
AT = np.transpose(A)
C = 1/2 * np.matmul(A, AT)
print(C)
eigenvalue, eigenvector = np.linalg.eig(C)
eigenvalueNEW = []
array = []
min = min(eigenvalue)
for i in range(len(eigenvalue)):
    if eigenvalue[i] > min:
        eigenvalueNEW.append(eigenvalue[i])
        array.append(i)


print("Eigenvalues of C = \n",eigenvalueNEW)
for i in array:
    print(eigenvector[i])