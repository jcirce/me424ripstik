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
y_2 = Function('x_2')(t)
z_2 = Function('x_2')(t)
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

#caster rotation 
R32 = rot2(theta_c)*rot3(theta_32)
R31 = R21*R32 #following paper even tho its backwards:/

#location of com of rear footpad in inertial frame (absolute position)
r_G2_1 = Matrix([
    [x_2],
    [y_2],
    [z_2]])

lz, lcx, lcz, lx1, lx2, lx3 = symbols('lz lcx lcz lx1 lx2 lx3')


#position of rear caster (body3) in 2-frame
r_G3_2 = Matrix([
    [-lx2],
    [0],
    [-lz]])

#position of rear caster in inertial frame
r_G3_1 =  r_G2_1 + R21*r_G3_2 


k = symbols('k')
J = symbols('J')
l = symbols('l')
m1, m2, g = symbols('m1 m2 g')