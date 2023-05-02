from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from tools import *

# LaTeX output
init_printing()
fname = current_dir + r'/chaplygin_sleigh_v1v2.tex'

# Define variables
t = symbols('t')

m = symbols('m')
J1, J2, J3 = symbols('J1, J2, J3')

b = symbols('b')
g = symbols('g')

u1 = Function('u_1')(t)
u2 = Function('u_2')(t)

# pseudo velocities
qdot = Matrix([u1, u2])
qddot = diff(qdot, t)

# Angular velocity 
omega01_1 = Matrix([[0],[0],[u2/b]])
Omega01_11 = skew(omega01_1)

# COM velocity
vC_1 = Matrix([[u1],[u2],[0]])

# Jacobian 
Xdot = Matrix([vC_1, omega01_1])
B = Xdot.jacobian(qdot)
Bdot = diff(B, t)

# No applied forces
G = zeros(6, 1)

# Setup matrices for projection
M = zeros(6,6)
M[:3,:3] = m*eye(3)
#print(type(J))
M[3:,3:] = diag(J1, J2, J3)
D = zeros(6,6)
D[0:3,0:3] = Omega01_11
D[3:6,3:6] = Omega01_11

Mstar = simplify(B.T * M * B)
Nstar = simplify(B.T * (D*M*B + M*Bdot))
Gstar = simplify(B.T * G)
eom = simplify(Mstar.inv() * (Gstar - Nstar*qdot))

original_stdout = sys.stdout # Save a reference to the original standard output
with open(fname, 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    # LaTeX output
    print(r'\documentclass{article}')
    print(r'\usepackage[margin=0.7in]{geometry}')
    print(r'\usepackage{lscape}')
    print(r'\usepackage[parfill]{parskip}')
    print(r'\usepackage[utf8]{inputenc}')
    print(r'\usepackage{amsmath,amssymb,amsfonts,amsthm}')
    print(r'\begin{document}')
    # print(r'\begin{landscape}')
    print(r'\begin{align}')
    print(r'B &= ', replace_values_in_string(latex(B)))
    print(r'\end{align}')
    print(r'\\')
    print(r'\begin{align}')
    print(r'eom &= ', replace_values_in_string(latex(eom)))
    print(r'\end{align}')
    print('\n', r'\end{document}')

    sys.stdout = original_stdout # Reset the standard output to its original value

