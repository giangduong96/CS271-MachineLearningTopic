"""
Name: Giang Duong
# Homework 1:
3. Write a HMM program for the English text problem in Section 9.2 of Chapter 9.
Test your each program on each of the following cases.

"""
import re
import string
import math
from sys import maxsize
import numpy as np
import random

def FABS(x):
    if x < 0.0:
        x = -x
    else:
        x = x
    return x

def setupfilename():
    # remove all characters except alphabets and word-space, save the new file to new_brown
    with open("brown.txt", "r") as filename, open("newbrown.txt", "w") as newline:
        x = filename.read()
        result = re.sub("[^a-z\s]", "", x, 0, re.IGNORECASE | re.MULTILINE)
        result = result.lower()
        newline.write(result)
        filename.close()

    occurrences = {}  # count the alphabet
    with open('newbrown.txt', 'r') as file:
        letters = string.ascii_lowercase
        for i in file:
            text_lower = i.lower()
            for letter in letters:
                if letter in text_lower:
                    occurrences[letter] = occurrences.get(letter, 0) + 1
        sorted(occurrences)
    # for word in occurrences:
    # print(word, ":", occurrences[word], "times.")
    # stringAZ = list(string.ascii_lowercase +' ')
    # print(stringAZ)

    count = 0
    with open('newbrown.txt', 'r') as file:
        while (count < countto):
            c = file.read(1)
            count += 1
            if not c:
                print("End of File")
                break
            char_to_num = ord(c) - 97
            if (char_to_num < 0):
                char_to_num = 26  # for white space (DEC = 32); 32-97 < 0 then put it to the last : a-z then space
            stringObservation.append(char_to_num)
    # print(stringObservation)
    return stringObservation


def getRandom(N):
    r = 1.0 / N
    diff = 0.1 * r
    t = r + (-diff + random.random() * 2 * diff)
    t = np.round(t, decimals=5)
    return t

def rowstocastic(numlist):
    s = sum(numlist)
    norm = [float(i)/s for i in numlist]
    norm = np.round(norm, decimals=5)
    return norm

def createTableA(number_of_state):
    A = []
    for i in range(number_of_state):
        A.append([])
        for j in range (number_of_state):
            A[i].append(getRandom(number_of_state))
        A[i] = rowstocastic(A[i])

    print(np.matrix(A))
    #A = [[intialA for x in range(number_of_state)] for x in range(number_of_state)]

    return A
def createTableB(number_of_state, symbols):
    B = []
    for i in range(number_of_state):
        B.append([])
        for j in range(symbols):
            B[i].append(getRandom(number_of_state))
        B[i] = rowstocastic(B[i])

    print(np.matrix(B))
    #B = [[intialB] * symbols] * number_of_state  # B = NxM matrix
    return B

def createTablepi(number_of_state):
    pi = []
    for i in range(number_of_state):
        pi.append(getRandom(number_of_state))
    pi = rowstocastic(pi)
    print(pi)
    #pi = [intialA, intialA]  # Pi = 1xN
    return pi

"""
Forward algorithm
"""


def forward():
    c[0] = 0.0
    for i in range(number_of_state):
        alpha[0][i] = (pi[i] * B[i][stringObservation[0]])
        c[0] += alpha[0][i]
    """scale alpha """
    c[0] = (1 / c[0])
    for i in range(number_of_state):
        alpha[0][i] = c[0] * alpha[0][i]
        # alpha[0][i] = alpha[0][i]/ c[0]

    """compute alpha (i)"""
    for t in range(1, T):
        c[t] = 0.0
        for i in range(number_of_state):
            alpha[t][i] = 0.0
            for j in range(number_of_state):
                alpha[t][i] = alpha[t][i] + alpha[t - 1][j] * A[j][i]
            alpha[t][i] = alpha[t][i] * B[i][stringObservation[t]]
            c[t] = c[t] + alpha[t][i]
        c[t] = 1 / c[t]
        # scale A[t][i]
        for i in range(number_of_state):
            alpha[t][i] = c[t] * alpha[t][i]
            # alpha[t][i] = alpha[t][i] / c[t]


""" Backward algorithm"""


def backward():
    Tbackward = T - 1
    for i in range(number_of_state):
        beta[Tbackward][i] = 1.0 * c[Tbackward]
    # compute beta (i)
    for t in reversed(range(T - 1)):
        for i in range(number_of_state):
            beta[t][i] = 0.0
            for j in range(number_of_state):
                beta[t][i] += (A[i][j] * B[j][stringObservation[t + 1]] * beta[t + 1][j])
            beta[t][i] = c[t] * beta[t][i]


""" Compute the gammas and di-gammas"""


def gammasAnddigammas():
    for t in range(0, T - 1):
        denom = 0.0
        temp2 = 0.0
        for i in range(number_of_state):
            gammas[t][i] = 0
            for j in range(number_of_state):
                digammas[t][i][j] = (alpha[t][i] * A[i][j] * B[j][stringObservation[t + 1]] * beta[t + 1][j])  # /denom
                gammas[t][i] += digammas[t][i][j]
            # check gammas, check if gamm[i] == alpha[i] *beta[i] / sum(alpha[j] *beta[j]
            temp2 += gammas[t][i]
            temp = 0.0
            for j in range(number_of_state):
                temp += alpha[t][j] * beta[t][j]
            temp = (alpha[t][i] * beta[t][i]) / temp
            if ((FABS(temp - gammas[t][i])) > EPSILON):
                print("gammas ", i, "=", gammas[t][i], temp, "Error!!!")
        if (FABS(1.0 - temp2) > EPSILON):
            print("Sum of gammas's = ", temp2, "should sum to 1.0. \n")

    # special case for gammas[T-1](i)
    for i in range(number_of_state):
        gammas[Tbackward][i] = alpha[Tbackward][i]


""" Re-estimate the model π"""


def reestimate():
    for i in range(number_of_state):
        pi[i] = gammas[0][i]

    # Re-estimate matrix A: Note: follow the row stocastic
    for i in range(number_of_state):
        for j in range(number_of_state):
            numer = 0.0
            denom = 0.0
            for t in range(Tbackward):
                numer += digammas[t][i][j]
                denom += gammas[t][i]
            # if numer == 0:
            #    A[i][j] = 0
            # else:
            A[i][j] = numer / denom

    # Re-estimate matrix B: Note: follow the row stocastic
    for i in range(number_of_state):
        for j in range(symbols):
            numer = 0.0
            denom = 0.0
            for t in range(T):
                if (stringObservation[t] == j):
                    numer += gammas[t][i]
                denom += gammas[t][i]
            # if numer == 0:
            #    B[i][j] = 0
            # else:
            B[i][j] = numer / denom
    return A, B, pi


def computeLog(c):
    """
    Compute log P(O|lambda)
    """
    logProb = 0
    for t in range(T):
        logProb += math.log(c[t])
    return -logProb

EPSILON = 0.001
symbols = 27  # 26 characters and word-space; M = 27
number_of_state = 2  # N = 2 because of vowel and consonant
intialA = 1 / number_of_state
intialB = 1 / symbols
T = 50000  # length of observation sequences
Tbackward = T - 1
countfrom = 0
countto = 50000
maxIters = 200
minIters = 20
c = np.zeros((T,), dtype=np.float32)
alpha = np.zeros((T, number_of_state), dtype=np.float32)
beta = np.zeros((T, number_of_state), dtype=np.float32)
gammas = np.zeros((T, number_of_state), dtype=np.float32)
digammas = np.zeros((T, number_of_state, number_of_state), dtype=np.float32)

stringObservation = []
stringObservation = setupfilename()
# observationsequence = symbols ** T
# observation sequence = M^T
A = createTableA(number_of_state)
B = createTableB(number_of_state, symbols)
pi = createTablepi(number_of_state)

iters = 0
logProb = -1.0
newLogProb = 0.0
diff = maxsize
while ((iters < maxIters) and (newLogProb > logProb)):
    print("Iteration: ", iters)
    logProb = newLogProb
    forward()
   #print("Alpha pass. Done!")
    backward()
    #print("Beta pass. Done!")
    gammasAnddigammas()
    #print("Gammas, Digammas pass. Done!")
    reestimate()
    #print("Reestimate pass. Done!")

    logProb = computeLog(c)
    # trick so no initial logProb is requires
    if iters == 0:
        logProb = newLogProb - 1.0
    #print("A=,", A)
    #print("B=,", B)
    #print("pi=,", pi)

    print("Iteration = %d, score log [P(observation |lambda)] = %s" % (iters, logProb))
    iters += 1

# HMM training finished
print("HMM training finished. Model:")
A = np.round(A, decimals=3)
B = np.round(B, decimals=3)
pi = np.round(pi, decimals=3)
print("Final  ==A==")
print(np.matrix(A))
print("Final  ==B== ")
print(np.matrix(B))
print("Final  ==pi== %s")
print(np.matrix(pi))

