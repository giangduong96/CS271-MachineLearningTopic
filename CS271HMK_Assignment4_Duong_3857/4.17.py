"""
Name: Giang Duong
ID : 014533857
# Homework 4.17 - PCA
"""
import numpy as np
import math
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)
# u1, u1 rows from 4.10
u1 = np.array([0.1641, 0.6278, -0.2604, -0.5389, 0.4637, 0.0752])
u2 = np.array([0.2443, 0.1070, -0.8017, 0.4277, -0.1373, -0.2904])
u1 = abs(u1)
u2 = abs(u2)
u3 = list(u1)
u4 = list(u2)
Uvectors = [u1, u2]
print("U1:", u1, "U2:", u2)
maxu1 = max(u1)
maxu2 = max(u2)
index1 = u3.index(maxu1)
index2 = u4.index(maxu2)
#17.a
print("a) Vector u1 has the greatest positive impact on the projection space is  ",maxu1, "at index", index1 + 1 )
print("a) Vector u2 has the greatest negative impact on the projection space is  ",maxu2, "at index", index2 + 1 )

print("17.b")
u1 = np.array([0.1641, 0.6278, -0.2604, -0.5389, 0.4637, 0.0752])
u2 = np.array([0.2443, 0.1070, -0.8017, 0.4277, -0.1373, -0.2904])
u3 = np.array([-0.0710, 0.2934, 0.3952, 0.3439, 0.3644, -0.7083])

#eigenvalues from (4.9)
e1 = 4.0833
e2 = 1.2364
e3 = 0.7428

l1 = math.sqrt(e1)*u1   # eigenvalue * u
l2 = math.sqrt(e2)*u2   # eigenvalue * u
l3 = math.sqrt(e3)*u3   # eigenvalue * u
print("L1 = %s" % l1)
print("L2 = %s" % l2)
print("L3 = %s" % l3)
print("L1^2 + L2^2 = %s" % ((l1**2)+(l2**2),))
print("L2^2 + L3^2 = %s" % ((l2**2)+(l3**2),))
print("L1^2 + L3^2 = %s" % ((l1**2)+(l3**2),))

print("17.c")
# LL1 = np.array([(l1**2)+(l2**2),])
# LL2 = np.array([(l2**2)+(l3**2),])
# LL3 = np.array([(l1**2)+(l3**2),])
# arrayLL = np.array(LL1, LL2, LL3, dtype=float)
# print("a",arrayLL)

imp = np.array((l1**2)+(l2**2), dtype=float)
imp2 = list(imp)
maximp = max(imp2)
minimp = min(imp2)
print(imp)
indeximpmax = imp2.index(maximp)
indeximpmin = imp2.index(minimp)
print("L1^2 + L2^2 = %s" % ((l1**2)+(l2**2),))
print("The most important feature is at index ",indeximpmax +1, "with value ", maximp)
print("The less important feature is at index ",indeximpmin +1, "with value ", minimp)
