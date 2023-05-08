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

def rotx(t):
    return Matrix([
        [1, 0, 0],
        [0, cos(t), -sin(t)],
        [0, sin(t), cos(t)]])

def roty(t):
    return Matrix([
        [cos(t), 0, sin(t)],
        [0, 1, 0],
        [-sin(t), 0, cos(t)]])

def rotz(t):
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

# x_2dot = diff(x_2, t)
# y_2dot = diff(y_2, t)
# z_2dot = diff(z_2, t)

alpha_2 = Function('alpha_2')(t)
beta_2 = Function('beta_2')(t)
gamma_2 = Function('gamma_2')(t)

# alpha_2dot = diff(alpha_2, t)
# beta_2dot = diff(beta_2, t)
# gamma_2dot = diff(gamma_2, t)

#twist forward deck to rear deck
theta_52 = Function('theta_52')(t)
# theta_52dot = diff(theta_52, t)

#relative angle caster to deck
theta_32 = Function('theta_32')(t)
theta_65 = Function('theta_65')(t)
# theta_32dot = diff(theta_32, t)
# theta_65dot = diff(theta_65, t)

#fixed caster angle
theta_c = symbols('theta_c')

#rotation of wheel
theta_43 = Function('theta_43')(t)
theta_76 = Function('theta_76')(t)
# theta_43dot = diff(theta_43, t)
# theta_76dot = diff(theta_76, t)

#non-rotating wheel angle
xi_R = Function('xi_R')(t)
xi_F = Function('xi_F')(t)
# xi_Rdot = diff(xi_R, t)
# xi_Fdot = diff(xi_F, t)

q = Matrix([x_2,y_2,z_2,alpha_2,beta_2,gamma_2,theta_52,theta_32,theta_65,theta_43,theta_76,xi_R,xi_F])
qdot = diff(q,t)
# qdot = Matrix([x_2dot,
#                y_2dot,
#                z_2dot,
#                alpha_2dot,
#                beta_2dot,
#                gamma_2dot,
#                theta_52dot,
#                theta_32dot,
#                theta_65dot,
#                theta_43dot,
#                theta_76dot,
#                xi_Rdot,
#                xi_Fdot])


#rotation from inertial to rear footpad
#from 1 to 2
R12 = rotx(alpha_2)*roty(beta_2)*rotz(gamma_2)

#rear caster
#from 2 to 3
R23 = roty(-theta_c)*rotz(theta_32)

R13 = R12*R23  #2's cancel

#rear wheel
R34 = rotz(theta_43)
R14 = R13*R34 

#front footpad rotation
R25 = rotx(theta_52)
R15 = R12*R25 #2's cancel

#front caster
R56 = roty(-theta_c)*rotz(theta_65)
R16 = R15*R56 

R26 = R25*R56 

#front wheel
R67 = rotz(theta_76)
R17 = R16*R67

R27 = R26*R67

#non-rotating wheel frames
Rxi_R = roty(xi_R)
Rxi_F = roty(xi_F)

# B = Matrix([[B2],[B3],[B4],[B5],[B6],[B7]])
# print(B3)

#location of com of rear footpad in inertial frame (absolute position)
r_G2_1 = Matrix([
    [x_2],
    [y_2],
    [z_2]])

v_G2_1 = diff(r_G2_1, t)
# print(v_G2_1)

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
r_G3_1 =  r_G2_1 + R12*r_G3_2 

#position rear wheel in caster frame
r_G4_3 = Matrix([
    [-l_cx],
    [0],
    [-l_cz]])

#position of rear wheel in 2 frame
r_G4_2 = R23*r_G4_3

#position of rear wheel in inertial frame
r_G4_1 = r_G3_1 + R13*r_G4_3

#position of front footpad(5) in rear(2)frame
r_G5_2 = Matrix([
    [l_x1],
    [0],
    [0]])

#front footpad(5) in inertial frame
r_G5_1 =  r_G2_1 + R12*r_G5_2 

#front caster(6) in front footpad 5-frame
r_G6_5 = Matrix([
    [l_x3],
    [0],
    [-l_z]])

#front caster in inertial
r_G6_2 = R25*r_G6_5
r_G6_1 = r_G5_1 + R15*r_G6_5 

#front wheel in caster frame
r_G7_6 = Matrix([
    [-l_cx],
    [0],
    [-l_cz]])

#front wheel in inertial
r_G7_1 = r_G6_1 + R16*r_G7_6
r_G7_2 = R26*r_G7_6
 
#location of wheel contact with ground
R = symbols('R')
wheel_radius = Matrix([
    [R],
    [0],
    [0]])

#position of contact point on rear and front wheels   
r_R = r_G4_1 + R14*Rxi_R*wheel_radius
r_F = r_G7_1 + R17*Rxi_F*wheel_radius

#tangent vectors of wheels at R and F
t_R = R14* diff(Rxi_R*wheel_radius,xi_R)
t_F = R17 *diff(Rxi_F*wheel_radius,xi_F)

#holonomic constraints on the system Ch(q) = [CH1,CH2,CH3,CH4].T
CH1 = r_R[2]
CH2 = r_F[2]
normVect = Matrix([0,0,1]) #normal vector to ground in 1 frame
CH3 = t_R.dot(normVect)
CH4 = t_F.dot(normVect)

CH = Matrix([CH1,CH2,CH3,CH4]) #holonomic constraints  CH = 0

#non holonomic constriants on the system Cnh(q) = [CNH1,CNH2,CNH3,CNH4].T
#velocity vectors of the contact points F and R
v_R = r_R.diff(t)
v_F = r_F.diff(t)

CNH1= v_R[0]
CNH2 = v_R[1]
CNH3 = v_F[0]
CNH4 = v_F[1]

CNH = Matrix([CNH1,CNH2,CNH3,CNH4]) #non holonomic constraints CNH=0

#finding the holonomic velocity of lagranges equation?? not sure what its called 
q = Matrix([x_2,y_2,z_2,alpha_2,beta_2,gamma_2,theta_52,theta_32,theta_65,theta_43,theta_76,xi_R,xi_F])

CH_dt = CH.jacobian(q) #CH_dt * qdot = 0
B = CNH.jacobian(q) #B * qdot = 0 #is the jacobian of the system

D = Matrix([[CH_dt,B]]) # it should be vertical [[CH_dt],[B]].T, plese double check   #D(q)qdot = 0


mp = 1.14 #kg mass of decks
mh=  70 #kg mass of human

m2 = mp+mh/2
m5 = m2
m3 = 0.25 #mass of casters
m6 = m3
m4 = 0.11 #mass of wheels
m7 = m4 

M_m = diag(m2,m3,m4,m5,m6,m7) #mass matrix
# M_I = []#MOI Tensor
# M_q = B*M*B


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
    print(r'R^{(1,2)} &= ', replace_values_in_string(latex(R12)))
    print(r'\end{align}')
    print(r'\\')


    # print(r'\begin{align}')
    # print(r'\tilde{\omega}_{1/2}^{(1,1)} &= ', replace_values_in_string(latex(Omega12_11)))
    # print(r'\end{align}')
    # print(r'\\')

    print(r'\begin{align}')
    print(r'r_{2}^{(1)} &= ', replace_values_in_string(latex(r_G2_1)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(3,1)} &= ', replace_values_in_string(latex(R13)))
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
    print(r'R^{(4,1)} &= ', replace_values_in_string(latex(R14)))
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
    print(r'R^{(5,1)} &= ', replace_values_in_string(latex(R15)))
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
    print(r'R^{(6,1)} &= ', replace_values_in_string(latex(R16)))
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
    print(r'R^{(7,1)} &= ', replace_values_in_string(latex(R17)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{7}^{(1)} &= ', replace_values_in_string(latex(r_G7_1)))
    print(r'\end{align}')

    print(r'\begin{align}')
    print(r'r_{7}^{(2)} &= ', replace_values_in_string(latex(r_G7_2)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(\xi_{R},4)} &= ', replace_values_in_string(latex(Rxi_R)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'R^{(\xi_{F},7)} &= ', replace_values_in_string(latex(Rxi_F)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{R}^{(1)} &= ', replace_values_in_string(latex(r_R)))
    print(r'\end{align}')
    print(r'\\')

    print(r'\begin{align}')
    print(r'r_{L}^{(1)} &= ', replace_values_in_string(latex(r_F)))
    print(r'\end{align}')
    print(r'\\')

    print('\n', r'\end{document}')

    sys.stdout = original_stdout # Reset the standard output to its original value