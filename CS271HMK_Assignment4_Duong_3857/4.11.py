"""
Name: Giang Duong
ID : 014533857
# Homework 4.11 - PCA
"""
import numpy as np
import math
from scipy.spatial import distance
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

print("11.a")
x1 = [2, -1, 0, 1, 1, -3, 5, 2]
x2 = [-2, 3, 2, 3, 0, 2, -1, 1]
x3 = [-1, 3, 3, 1, -1, 4, 5, 2]
x4 = [3, -1, 0, 3, 2, -1, 3, 0]

# join 4 vector into 1 matrix
B = np.array([x1,x2,x3,x4], dtype=float)
B = np.matrix.transpose(B)
# Note B.T = B.np.matrix.transpose(B)
print("B matrix = \n", B);
rowmean= []
for i in range(8):
    sum = np.sum(B[i]) / 4
    rowmean.append(sum)
    for j in range(4):
        B[i][j] = B[i][j] - sum
print(B)
print("rowmean" ,rowmean)
# assign B matrix to A

A = B
print("B matrix after mean = 0 \n", A)
# create A transpose (AT)
AT = np.matrix.transpose(A)
# C = 1/n (A AT)
C = 1 / (len(x1)) * np.matmul(A, AT)
print("Covariance matrix = \n", C)
# eigenvalues and eigenvector of C
eigenvalue, eigenvector = np.linalg.eig(C)
eigenvalue= eigenvalue.real
eigenvector= eigenvector.real
print("Eigenvalues of C = \n", eigenvalue)
print("Eigenvector of C = \n", eigenvector)

# SVD of C
u, s, vh = np.linalg.svd(B,full_matrices=True)
print("U : \n", u)
print("S : \n", s)
print("V : \n", vh)

# scoring matrix based on the 3 most significant eigenvector of C
significantEigenvector = []
print("The 3 most significant eigenvector of C is : \n")
# minEigenvector = sorted(eigenvalue)
# minEigenvector = minEigenvector[len(A)-((len(A)-4))]
# for i in range(len(eigenvalue)):
#    if eigenvalue[i] > minEigenvector: significantEigenvector.append(i)
# print(significantEigenvector)
array = [0,1,2]
newdeltaScoreMatrix = []
deltaScoreMatrix = np.matmul(np.matrix.transpose(u), A)
for i in array:
    newdeltaScoreMatrix.append(deltaScoreMatrix[i])
# print(newdeltaScoreMatrix)
#newdeltaScoreMatrix[0] = newdeltaScoreMatrix[0] * -1
#newdeltaScoreMatrix[2] = newdeltaScoreMatrix[2] * -1
#newdeltaScoreMatrix = np.multiply(newdeltaScoreMatrix,  -1)
print(np.array(newdeltaScoreMatrix))

print("11.b")
Y1 = [1, 5, 1, 5, 5, 1, 1, 3]
Y2 = [-2, 3, 2, 3, 0, 2, -1, 1]
Y3 = [2, -3, 2, 3, 0, 0, 2, -1]
Y4 = [2, -2, 2, 2, -1, 1, 2, 2]


def scorePCA(Y, uvectors, newdelta,rowmean):
    # convert Y into array and transpose matrix
    print("Scoring phrase:")
    Y = np.array(Y, dtype=float)
    Y = np.matrix.transpose(Y)
    # Note: Use row mean from the matrix A to score
    for i in range(len(Y)):
        Y[i] = Y[i] - rowmean[i]
    print("Y~: \n",Y)

    newdelta =np.array(newdelta) # convert newdelta into array list
    # W = Y.u(i)
    W = np.array([np.dot(Y,u) for u in uvectors])
    print("W = \n", W)
    # # only take 3
    # W2 = W[:3]
    # print("W2 = \n", W2)
    # # change to array, transpose and take only 3 from delta
    # newdelta = np.array(newdelta)
    # newdelta = np.matrix.transpose(newdelta)
    # # newdelta = newdelta[:4]
    # #print(newdelta)
    score = []
    for i in range(newdelta.shape[1]):
        distance = math.sqrt(np.sum((W - newdelta[:,i])**2))
        score.append(distance)
    print(score)
    k = min(score)
    print("Emin = ", k)
    return k


Uvectors = []
for i in array:Uvectors.append(u[:,i])
Uvectors= np.array(Uvectors)
print("UVector choose: \n", Uvectors)

#
scorePCA(Y1, Uvectors, newdeltaScoreMatrix,rowmean)
scorePCA(Y2, Uvectors, newdeltaScoreMatrix,rowmean)
scorePCA(Y3, Uvectors, newdeltaScoreMatrix,rowmean)
scorePCA(Y4, Uvectors, newdeltaScoreMatrix,rowmean)
