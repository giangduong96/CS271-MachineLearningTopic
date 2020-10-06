"""
Name: Giang Duong
ID : 014533857
# Homework 4.1 - PCA
"""
# np.linalg is library for linear algebra
import numpy as np

A = [[1, 0, -2],
     [-2, 3, 1]]
B = [[1, -1],
     [2, -2],
     [3, 1]]
C = [[3, 2],
     [-1, -3],
     [2, -1]]
print("a. 2A = ")
aa = np.array(A) * 2
print(aa)

print("b. B +C = ")
b_and_c = np.array(B) + np.array(C)
print(b_and_c)

print("c. A+B is undefined because A 2x3 and B is 3x2. To add matrix with another matrix we need nxm = sxt ")
#a_and_b = np.array(A) + np.array(B)
#print(a_and_b) will create an error

print("d. AB = ")
axb = np.matmul(A,B)
print(axb)

print("e. BA = ")
bxa = np.matmul(B,A)
print(bxa)

print("f. BC is undefined because B 3x2 and C is 3x2. \nTo get product of 2 matrices we need m = s where B nxm = C sxt ")

