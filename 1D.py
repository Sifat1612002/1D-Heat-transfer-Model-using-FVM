import numpy as np
import matplotlib.pyplot as plt

# This program will solve 1D Heat Transfer through a bar

# inputs
T_A = 100  # Boundary condition
T_B = 300  # Boundary condition
k = 100  # Conductivity
node = 10  # Number of nodes
L = 500  # Length of the wire
S = 1000  # Heat generation rate in the wire

# properties
d = L / node  # length of each node
D = k / d
A = 0.1
V = A * d

coe = np.zeros((node, node))  # coefficient matrix

# solving to get the temperature at nodes
for n in range(1, node + 1):
    if n == 1:
        a_L = 0
        a_R = D * A
        S_p = -2 * D * A
        S_u = T_A * (2 * D * A) + (S * V)
        a_p = a_L + a_R - S_p

        coe[n - 1, n - 1] = a_p
        coe[n - 1, n] = -a_R

        ans = np.array([S_u])

    elif n == node:
        a_R = 0
        a_L = D * A
        S_p = -2 * D * A
        S_u = T_B * (2 * D * A) + (S * V)
        a_p = a_L + a_R - S_p

        coe[n - 1, n - 1] = a_p
        coe[n - 1, n - 2] = -a_L
        ans = np.append(ans, S_u)

    elif 1 < n < node:
        S_p = 0
        S_u = S * V
        a_R = A * D
        a_L = A * D
        a_p = a_L + a_R - S_p

        coe[n - 1, n - 2] = -a_L
        coe[n - 1, n - 1] = a_p
        coe[n - 1, n] = -a_R

        ans = np.append(ans, S_u)

Temperature_Dist = np.linalg.solve(coe, ans)

# plotting discrete points and analytical results for comparison
x_value = np.zeros(node)
y_value = np.zeros(node)

for i in range(node):
    x1 = (d / 2) + (i) * (d)
    y1 = Temperature_Dist[i]
    x_value[i] = x1
    y_value[i] = y1

plt.plot(x_value, y_value, '*', markersize=10)
x = np.arange(0, L + 0.05, 0.05)
analytical_T = T_A + (x * (T_B - T_A) / L) + (S * x * (L - x)) / (2 * k)
plt.plot(x, analytical_T)

# comparing the results
matrix_of_analytical_T = np.zeros(node)
for i in range(node):
    x = x_value[i]
    matrix_of_analytical_T[i] = T_A + (x * (T_B - T_A) / L) + (S * x * (L - x)) / (2 * k)

table_data = np.column_stack((x_value, y_value, matrix_of_analytical_T, matrix_of_analytical_T - y_value))
table_header = ['X', 'TEMP', 'Analytical_Temp', 'ERROR']
print("{:<10} {:<10} {:<20} {:<20}".format(*table_header))
for row in table_data:
    print("{:<10} {:<10} {:<20} {:<20}".format(*row))
plt.show()
