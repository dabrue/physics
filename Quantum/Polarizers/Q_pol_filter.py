#!/usr/bin/python
"""
    Quantum Polarizers in Series

    This script looks at two questions: 

    1. Three polarizing filters in series, with the first at 0 degrees and the
        last at 45 degrees, what angle should the middle one be to maximize throughput? 

    2. In the same setup with question 1, where the input is light polarized to 0
        degrees and the output is at 45 degrees, what is the ideal number of 
        polarizers to put between them to maximize output, and at what relative angles?

    3. For an arbitrary angle difference between input and output, what is the ideal
        number of filters and relative angles for maximum output? 

    Questions 2 and 3 can probably be answered simultaneously if coded correctly. 
"""

__author__ = "Daniel A. Brue"
__status__ = "Development"

import numpy
import math
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc("text",usetex=True)

pi=math.pi

Text1=["\n\n",
" Question 1: If filter 1 is at 0 degrees and filter 3 is at 45 degres, what angle",
"   should filter 2 be set to maximize light transmission?",
"",
" Answer 1: The probability of transmission of light with a given linear polarization",
"   through a polarizing lense that is not at the same angle is given by",
"",
"   P = (cos(theta))^2",
"",
"   where theta is the angle between the polarizer and the known initial polarization",
"   of the light.",
"",
"   If the light does pass through in this new basis, then it assumes the new",
"   polarization direction. ",
"",
"   In this case, phi will be the angle of polarizer 2, and the probability of ",
"   transmission between all three filters will be given by ",
"",
"   P3 = (cos(phi))^2 * (cos(45-phi))^2",
"",
"   Which can be plotted as a function of phi... \n\n"]

for line in Text1:
    print line

phimin=0.0
phimax=pi
nphi = 101

phi=numpy.linspace(phimin, phimax, nphi)
P3 = numpy.zeros_like(phi)
piovr4=pi/4.0
for i in range(nphi):
    P3[i] = math.cos(phi[i])**2  * math.cos((piovr4 - phi[i]))**2

xtix=[0.0, pi/8.0, pi/4.0, 3*pi/8.0, pi/2.0, 5*pi/8.0, 3*pi/4.0, 7*pi/8.0, pi]
xlbl=["0", r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$", r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$", r"$\frac{5\pi}{8}$", r"$\frac{3\pi}{4}$", r"$\frac{7\pi}{8}$", r"$\pi$"]


fig=plt.figure()
ax=fig.add_subplot(1,1,1)
plt.plot(phi,P3)
plt.xticks(xtix,xlbl)
plt.ylabel("Transmission Intensity")
plt.xlabel(r"\Phi")
plt.ylim(0,1.1)
plt.suptitle("Transmission as a function of inner filter rotation")
plt.grid()
plt.show()
    

Text2=["\n\n",
" Question 2: What number of intermediate polarizers will maximize the total throughput?",
"\n",
" Answer 2: We saw in question 1 that when the middle polarizer was oriented to an angle",
"   half way between the first and third polarizers, then transmission had the highest",
"   probability for a single photon and the transmission intensity was highest for an ",
"   aggregate light beam.",
"\n",
"   Assuming each filter is angled half way between the filters before and after it,",
"   find the best possible number for a 45 degree difference between initial and final.\n\n"
]

for line in Text2:
    print line

MaxFilters = 50  # Total number of filters in between the first and last one. 
PhiStart = 0.0
PhiEnd   = pi/4.0

N=[]
P=[]
for i in range(MaxFilters+1):
    N.append(i)
    dphi = (PhiEnd-PhiStart)/(i+1)
    P.append((math.cos(dphi)**2)**(i+1))

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
plt.plot(N,P,"bo-")
plt.suptitle("Polarization from 0 to 45 degrees")
plt.xlabel("Number of inner filters")
plt.ylabel("Transmission Intensity")
plt.ylim(0,1.1)
plt.grid()
plt.show()


Text3=["\n\n",
" Question 3: How does the ideal number of filters change with respect to final angle?\n"
]

for line in Text3:
    print line

PhiEndRange = [15,30,45,60,75,90]



Ps=[]
for PE in PhiEndRange:
    PERad = PE*pi/180.0
    P=[]
    for i in range(MaxFilters+1):
        dphi=(PERad-PhiStart)/(i+1)
        P.append((math.cos(dphi)**2)**(i+1))
    Ps.append(P)

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
for i in range(len(Ps)):
    plt.plot(N,Ps[i],label=str(PhiEndRange[i]))
plt.legend(loc=4)
plt.xlabel("Number of inner filters")
plt.ylabel("Transmission Intensity")
plt.ylim(0,1.1)
plt.suptitle("Number of filters for different final angles")
plt.grid()
plt.show()
