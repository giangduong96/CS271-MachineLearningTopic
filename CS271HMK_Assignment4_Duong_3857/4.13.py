"""
Name: Giang Duong
ID : 014533857
# Homework 4.13 - PCA
"""
import numpy as np
import math
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)
# mean for each row from 4.8
meanrow4_8 = np.array([7.0/4, 7.0/4, 5.0/4, 2, 2, 1])
# u1, u1 rows from 4.10
u1 = np.array([0.1641, 0.6278, -0.2604, -0.5389, 0.4637, 0.0752])
u2 = np.array([0.2443, 0.1070, -0.8017, 0.4277, -0.1373, -0.2904])
Uvectors = [u1, u2]

print("question 13.a")
M1 = np.array([1, -1, 1, -1, -1, 1])
M2 = np.array([-2, 2, 2, -1, -2, 2])
M3 = np.array([1, 3, 0, 1, 3, 1])
M4 = np.array([2, 3, 1, 1, -2, 0])

B = np.array([M1,M2,M3,M4], dtype=float)
B = np.matrix.transpose(B)
# Note B.T = B.np.matrix.transpose(B)
print("B matrix = \n", B);

for i in range(6):
    for j in range(4):
        B[i][j] = B[i][j] - meanrow4_8[i]

# assign B matrix to A

A = B
print("B matrix after mean = 0 \n", A)
# create A transpose (AT)
AT = np.matrix.transpose(A)
# C = 1/n (A AT)
C = 1 / (len(M1)) * np.matmul(A, AT)
print("Covariance matrix = \n", C)
# eigenvalues and eigenvector of C
#eigenvalue, eigenvector = np.linalg.eig(C)
#eigenvalue= eigenvalue.real
#eigenvector= eigenvector.real
#print("Eigenvalues of C = \n", eigenvalue)
#print("Eigenvector of C = \n", eigenvector)

# SVD of C
u, s, vh = np.linalg.svd(B,full_matrices=True)
#print("U : \n", u)
#print("S : \n", s)
#print("V : \n", vh)

# scoring matrix based on the 3 most significant eigenvector of C
print("deltaScoreMatrix 13A : \n")
newdeltaScoreMatrixA = np.concatenate([u.reshape((1,-1)) for u in Uvectors], axis=0)
newdeltaScoreMatrixA = np.dot(newdeltaScoreMatrixA, A)
print(np.array(newdeltaScoreMatrixA))

print("13b")
D1 = np.array([-1, 2, 1, 2, -1, 0])
D2 = np.array([-2, 1, 2, 3, 2, 1])
D3 = np.array([-1, 3, 0, 1, 3, -1])
D4 = np.array([0, 2, 3, 1, 1, -2])

B1 = np.array([D1,D2,D3,D4], dtype=float)
B1 = np.matrix.transpose(B1)
# Note B.T = B.np.matrix.transpose(B)
print("B1 matrix = \n", B1);

for i in range(6):
    for j in range(4):
        B1[i][j] = B1[i][j] - meanrow4_8[i]

# assign B matrix to A1
A1 = B1
print("B1 matrix after mean = 0 \n", A1)
# create A transpose (AT)
A1T = np.matrix.transpose(A1)
# C = 1/n (A AT)
C1 = 1 / (len(D1)) * np.matmul(A1, A1T)
print("Covariance matrix = \n", C1)


# SVD of C
u, s, vh = np.linalg.svd(B,full_matrices=True)


# scoring matrix based on the 3 most significant eigenvector of C
print("deltaScoreMatrix 13B: \n")
newdeltaScoreMatrixB = np.concatenate([u.reshape((1,-1)) for u in Uvectors], axis=0)
newdeltaScoreMatrixB = np.dot(newdeltaScoreMatrixB, A1)
print(np.array(newdeltaScoreMatrixB))

def scorePCA(Y, uvectors, newdelta,rowmean):
    # convert Y into array and transpose matrix
    #print("Scoring phrase:")
    Y = np.array(Y, dtype=float)
    Y = np.matrix.transpose(Y)
    # Note: Use row mean from the matrix A to score
    for i in range(len(Y)):
        Y[i] = Y[i] - rowmean[i]
    #print("Y~: \n",Y)

    newdelta =np.array(newdelta) # convert newdelta into array list
    # W = Y.u(i)
    W = np.array([np.dot(Y,u) for u in uvectors])
    #print("W = \n", W)
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
    k = min(score)
    #print("Emin = ", k )
    return k



# 13 c
"""Input : Test set. Output: The min when scoring test set vs 2 scoring matrix from 13a abd 13 b"""
def detectMalware(Y, uvectors, newdeltaA, newdeltaB,rowmean):
    score13a = scorePCA(Y, uvectors, newdeltaA,rowmean)
    print("Score from 13a:\n", score13a)
    score13b = scorePCA(Y, uvectors, newdeltaB,rowmean)
    print("Score from 13b:\n", score13b)
    if (score13a <= score13b):
        print("Malware\n")
    else: print("Benign\n")

Y1 = np.array([1, 5, 1, 5, 5, 1])
Y2 = np.array([-2, 3, 2, 3, 0, 2])
Y3 = np.array([2, -3, 2, 3, 0 ,0])
Y4 = np.array([2, -2, 2, 2, -1, 1])

print("Y1")
detectMalware(Y1, Uvectors, newdeltaScoreMatrixA, newdeltaScoreMatrixB, meanrow4_8)
print("Y2")
detectMalware(Y2, Uvectors, newdeltaScoreMatrixA, newdeltaScoreMatrixB, meanrow4_8)
print("Y3")
detectMalware(Y3, Uvectors, newdeltaScoreMatrixA, newdeltaScoreMatrixB, meanrow4_8)
print("Y4")
detectMalware(Y4, Uvectors, newdeltaScoreMatrixA, newdeltaScoreMatrixB, meanrow4_8)

