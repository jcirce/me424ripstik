from sympy import *
import matplotlib.pyplot as plt
from tools import * 
import sys
import os

# LaTeX file output
fname = os.path.dirname(os.path.realpath(__file__)) + r'/pend_cart.tex'

# Define variables
t = symbols('t')
x2 = Function('x_2')(t)
theta = Function('theta')(t)

k = symbols('k')
J = symbols('J')
l = symbols('l')
m1, m2, g = symbols('m1 m2 g')

init_printing()
def skew(a):
    return Matrix([
        [0, -a[2], a[1]],
        [a[2], 0, -a[0]],
        [-a[1], a[0], 0]])
def unskew(a):
    return Matrix([
        [a[2,1]],
        [a[0,2]],
        [a[1,0]]])

def rot1(t):
    return Matrix([
        [1, 0, 0],
        [0, cos(t), -sin(t)],
        [0, sin(t), cos(t)]])

def rot2(t):
    return Matrix([
        [cos(t), 0, sin(t)],
        [0, 1, 0],
        [-sin(t), 0, cos(t)]])

def rot3(t):
    return Matrix([
        [cos(t), -sin(t), 0],
        [sin(t), cos(t), 0],
        [0, 0, 1]])

# Set derivatives
thetadot = diff(theta, t)
x2dot = diff(x2, t)
q = Matrix([x2, theta])
qdot = diff(q, t)
qddot = diff(qdot, t)

# Rotation matrix
R01 = rot3(theta)

# Angular velocity 
Omega01_11 = simplify(R01.T * diff(R01,t))
omega01_1 = unskew(Omega01_11)

# Body 1 velocity
v1C_0 = Matrix([[0],[x2dot],[0]])
# Body 2 velocity
v2C_0 = v1C_0 + R01 *  Omega01_11 * Matrix([[l],[0],[0]])

# Jacobians /  partial velocities
B1 = v1C_0.jacobian([x2dot, thetadot])
B2 = v2C_0.jacobian([x2dot, thetadot])
# Jacobian matrix and derivative
B = Matrix([[B1],[B2]])
Bdot = diff(B, t)

# External (non-constraint) forces 
Fs = -k*x2 # in negative 2 dir.
G1_0 = Matrix([
    [m1*g],
    [Fs],
    [0]])
G2_0 = Matrix([
    [m2*g],
    [0],
    [0]])

# Setup matrices for projection
M = zeros(6,6)
M[:3,:3] = m1*eye(3)
M[3:,3:] = m2*eye(3)
G = zeros(6,1)
G[:3,:] = G1_0
G[3:,:] = G2_0

Mstar = simplify(B.T * M * B)
Nstar = simplify(B.T * M * Bdot)
Gstar = simplify(B.T * G)
eom = Mstar * qddot + Nstar * qdot - Gstar # = 0

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
    print(r'\tilde{\omega}_{1/0}^{(1,1)} &= ', replace_values_in_string(latex(Omega01_11)))
    print(r'\end{align}')
    print(r'\\')
    print(r'\begin{align}')
    print(r'{\omega}_{1/0}^{(1)} &= ', replace_values_in_string(latex(omega01_1)))
    print(r'\end{align}')
    print(r'\\')    
    print(r'\begin{align}')
    print(r'B &= ', replace_values_in_string(latex(B)))
    print(r'\end{align}')
    print(r'\\')
    print(r'\begin{align}')
    print(r'M^{*} &= ', replace_values_in_string(latex(Mstar)))
    print(r'\end{align}')
    print(r'\\')
    print(r'\begin{align}')
    print(r'N^{*} &= ', replace_values_in_string(latex(Nstar)))
    print(r'\end{align}')    
    print(r'\\')
    print(r'\begin{align}')
    print(r'G^{*} &= ', replace_values_in_string(latex(Gstar)))
    print(r'\end{align}')     
    print(r'\\')
    print(r'\begin{align}')
    print(r'\textrm{EOM: } 0 &= ', replace_values_in_string(latex(eom)))
    print(r'\end{align}') 
    print('\n', r'\end{document}')

    sys.stdout = original_stdout # Reset the standard output to its original value


