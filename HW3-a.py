import math
A = [[0.7, 0.3],
     [0.4, 0.6]]
B = [[0.1, 0.4, 0.5],
     [0.7, 0.2, 0.1]]
pi= [0.6, 0.4]

S = 0;M = 1;L = 2
Set1 = [0, 1, 2]

O = []
sum = []
totalsum = 0

K = []
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                O.append([a, b, c, d])
print(O)



for i in range(0, 81):
    prob = []
    #HHHH
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][0]*B[0][O[i][1]])*(A[0][0]*B[0][O[i][2]])*(A[0][0]*B[0][O[i][3]]))
    #HHHC
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][0]*B[0][O[i][1]])*(A[0][0]*B[0][O[i][2]])*(A[0][1]*B[1][O[i][3]]))
    #HHCH
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][0]*B[0][O[i][1]])*(A[0][1]*B[1][O[i][2]])*(A[1][0]*B[0][O[i][3]]))
    #HHCC
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][0]*B[0][O[i][1]])*(A[0][1]*B[1][O[i][2]])*(A[1][1]*B[1][O[i][3]]))
    #HCHH
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][1]*B[1][O[i][1]])*(A[1][0]*B[0][O[i][2]])*(A[0][0]*B[0][O[i][3]]))
    #HCHC
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][1]*B[1][O[i][1]])*(A[1][0]*B[0][O[i][2]])*(A[0][1]*B[1][O[i][3]]))
    #HCCH
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][1]*B[1][O[i][1]])*(A[1][1]*B[1][O[i][2]])*(A[1][0]*B[0][O[i][3]]))
    #HCCC
    prob.append((pi[0]*B[0][O[i][0]])*(A[0][1]*B[1][O[i][1]])*(A[1][1]*B[1][O[i][2]])*(A[1][1]*B[1][O[i][3]]))
    #CHHH
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][0]*B[0][O[i][1]])*(A[0][0]*B[0][O[i][2]])*(A[0][0]*B[0][O[i][3]]))
    #CHHC
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][0]*B[0][O[i][1]])*(A[0][0]*B[0][O[i][2]])*(A[0][1]*B[1][O[i][3]]))
    #CHCH
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][0]*B[0][O[i][1]])*(A[0][1]*B[1][O[i][2]])*(A[1][0]*B[0][O[i][3]]))
    #CHCC
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][0]*B[0][O[i][1]])*(A[0][1]*B[1][O[i][2]])*(A[1][1]*B[1][O[i][3]]))
    #CCHH
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][1]*B[1][O[i][1]])*(A[1][0]*B[0][O[i][2]])*(A[0][0]*B[0][O[i][3]]))
    #CCHC
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][1]*B[1][O[i][1]])*(A[1][0]*B[0][O[i][2]])*(A[0][1]*B[1][O[i][3]]))
    #CCCH
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][1]*B[1][O[i][1]])*(A[1][1]*B[1][O[i][2]])*(A[1][0]*B[0][O[i][3]]))
    #CCCC
    prob.append((pi[1]*B[1][O[i][0]])*(A[1][1]*B[1][O[i][1]])*(A[1][1]*B[1][O[i][2]])*(A[1][1]*B[1][O[i][3]]))
    print("This is the observation set :", O[i])
    #print(prob)

    sum = math.fsum(prob)
    print("Sum = : ", sum, "\n")
    totalsum += sum

print("Total Sum = " , totalsum)

