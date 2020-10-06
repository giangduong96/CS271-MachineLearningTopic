
A = [[0.7, 0.3],
     [0.4, 0.6]]
B = [[0.1, 0.4, 0.5],
     [0.7, 0.2, 0.1]]
pi = [0.6, 0.4]

S = 0;
M = 1;
L = 2
Set1 = [0, 1, 2]

O = []
totalsum = 0

for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                O.append([a, b, c, d])

N = len(pi)  # total states
M = len(O[0])  # total observations
#print(N, M)
#print(O)
print("M = ",M)
# compute alpha zero
for step in range(81):
    alpha = [[0.0 for i in range(N)] for _ in range(M)]
    sum = 0.0
    # print(alpha)
    for i in range(N):
        alpha[0][i] = (pi[i] * B[i][O[step][0]])
    # compute alpha (i)
    for t in range(1, M):
        for i in range(N):
            alpha[t][i] = 0
            for j in range(0, N):
                alpha[t][i] = alpha[t][i] + alpha[t - 1][j] * A[j][i]
            alpha[t][i] = alpha[t][i] * B[i][O[step][t]]
            sum += alpha[M-1][i]
        totalsum += sum
    print("Observation : ", O[step])
    print("Alpha = ",alpha)
    print("Sum = ", sum)


print("Total Sum = ",totalsum)
