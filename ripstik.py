from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import animation
from scipy.integrate import solve_ivp
from tools import * 
import sys
import os

fname = os.path.dirname(os.path.realpath(__file__)) + r'/ripstik.tex'
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

#Define variables
t = symbols('t')

#for rear footpad (body 2)
x_2 = Function('x_2')(t)
y_2 = Function('y_2')(t)
z_2 = Function('z_2')(t)
alpha_2 = Function('alpha_2')(t)
beta_2 = Function('beta_2')(t)
gamma_2 = Function('gamma_2')(t)

#twist forward deck to rear deck
theta_52 = Function('theta_52')(t)

#relative angle caster to deck
theta_32 = Function('theta_32')(t)
theta_65 = Function('theta_65')(t)
#fixed caster angle
theta_c = symbols('theta_c')

#rotation of wheel
theta_43 = Function('theta_43')(t)
theta_76 = Function('theta_76')(t)

#rotation from inertial to rear footpad
R21 = rot1(alpha_2)*rot2(beta_2)*rot3(gamma_2)
R32 = rot2(-theta_c)*rot3(theta_32)
R31 = R21*R32 #following paper even tho its backwards:/
R41 = R31*rot3(theta_43)

#front footpad rotation
R52 = rot1(theta_52)
R51 = R21*R52
R61 = R51*rot2(-theta_c)*rot3(theta_65) 
R62 = R52*rot2(-theta_c)*rot3(theta_65) 
R71 = R61*rot3(theta_76)
R72 = R62*rot3(theta_76)

#location of com of rear footpad in inertial frame (absolute position)
r_G2_1 = Matrix([
    [x_2],
    [y_2],
    [z_2]])

l_z, l_cx, l_cz, l_x1, l_x2, l_x3 = symbols('l_z l_cx l_cz l_x1 l_x2 l_x3')
#lz = vertical difference between footpad and caster
#lcx = x distance between caster com and wheel com
#lcz = z ""
#lx1 = distance between two footpad's coms
#lx2 = distance between G2 and G3
#lx3 = distance between G5 and G6

#position of rear caster (body3) in 2-frame
r_G3_2 = Matrix([
    [-l_x2],
    [0],
    [-l_z]])

#position of rear caster in inertial frame
r_G3_1 =  r_G2_1 + R21*r_G3_2 

#position rear wheel in caster frame
r_G4_3 = Matrix([
    [-l_cx],
    [0],
    [-l_cz]])

#position of rear wheel in 2 frame
r_G4_2 = R32*r_G4_3

#position of rear wheel in inertial frame
r_G4_1 = r_G3_1 + R31*r_G4_3

#position of front footpad(5) in rear(2)frame
r_G5_2 = Matrix([
    [l_x1],
    [0],
    [0]])

#front footpad(5) in inertial frame
r_G5_1 =  r_G2_1 + R21*r_G5_2 

#front caster(6) in front footpad 5-frame
r_G6_5 = Matrix([
    [l_x3],
    [0],
    [-l_z]])

#front caster in inertial
r_G6_2 = R52*r_G6_5
r_G6_1 = r_G5_1 + R51*r_G6_5 

#front wheel in caster frame
r_G7_6 = Matrix([
    [-l_cx],
    [0],
    [-l_cz]])

#front wheel in inertial
r_G7_1 = r_G6_1 + R61*r_G7_6
r_G7_2 = R62*r_G7_6


original_stdout = sys.stdout # Save a reference to the original standard output
fname = os.path.dirname(os.path.realpath(__file__)) + r'/ripstik.tex'
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
    print(r'R^{(2,1)} &= ', replace_values_in_string(latex(R21)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{2}^{(1)} &= ', replace_values_in_string(latex(r_G2_1)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(3,1)} &= ', replace_values_in_string(latex(R31)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{3}^{(1)} &= ', replace_values_in_string(latex(r_G3_1)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{3}^{(2)} &= ', replace_values_in_string(latex(r_G3_2)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(4,1)} &= ', replace_values_in_string(latex(R41)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{4}^{(1)} &= ', replace_values_in_string(latex(r_G4_1)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{4}^{(2)} &= ', replace_values_in_string(latex(r_G4_2)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(5,1)} &= ', replace_values_in_string(latex(R51)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{5}^{(1)} &= ', replace_values_in_string(latex(r_G5_1)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{5}^{(2)} &= ', replace_values_in_string(latex(r_G5_2)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(6,1)} &= ', replace_values_in_string(latex(R61)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{6}^{(1)} &= ', replace_values_in_string(latex(r_G6_1)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{6}^{(2)} &= ', replace_values_in_string(latex(r_G6_2)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(7,1)} &= ', replace_values_in_string(latex(R71)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{7}^{(1)} &= ', replace_values_in_string(latex(r_G7_1)))
    print(r'\end{align}')

    print(r'\begin{align}')
    print(r'r_{7}^{(2)} &= ', replace_values_in_string(latex(r_G7_2)))
    print(r'\end{align}')
    print(r'\\')

    
    print('\n', r'\end{document}')

    sys.stdout = original_stdout # Reset the standard output to its original value