# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:35:38 2016

@author: nabraham
"""

oplsaa = '/home/nabraham/tinker/params/oplsaa.prm'

import os
import sys
import subprocess
import numpy as np
import itertools as it

#################################################################################
#################             Tinker Outputs          ###########################
#################################################################################

def eig_wvn(fil,key):
    #Calculates the eigenvalues and wavenumbers for a particular structre using Tinker
    #fil is the Tinker .xyz file for the desired output
    #key is the corresponding Tinker .key file for
    eigwvn = subprocess.check_output("vibrate %s -k %s  CR |  grep -oP '[-+]*[0-9]*\.[0-9]{2,9}'"%(fil,key),shell=True)
    eigwvn = eigwvn.split('\n')
    eigwvn2 = []
    for i in eigwvn:
      if i == '':
        pass
      else:
        eigwvn2.append(float(i))
    eig = np.array(eigwvn2[:len(eigwvn2)/2])
    wvn = np.array(eigwvn2[len(eigwvn2)/2:])
    return eig, wvn

def lat_ener(fil,key):
    #Extracts the lattice energy of the lattice
    #fil is the Tinker .xyz file for the desired output
    #key is the corresponding Tinker .key file for
    lat_ener = subprocess.check_output("analyze %s -k %s E | grep 'Total'| grep -oP '[-+]*[0-9]*\.[0-9]*'"%(fil,key),shell=True)
    lat_ener = float(lat_ener)
    return lat_ener

def num_atm(fil,nmol):
    ###Takes a Tinker .xyz file and determines the number of atoms per molecule
    #fil is the Tinker .xyz file for the desired output
    #nmol is the number of molecules per cell specified
    with open('%s'%(fil),'r') as l:
      atm = [lines.split() for lines in l]
      atm = np.array(list(it.izip_longest(*atm,fillvalue=' '))).T
    atm = int(atm[0,0])/nmol
    return atm

#################################################################################
#################          Force Field Inputs         ###########################
#################################################################################
def mass(fil,nmol):
    ###Outputs an array of the masses of the molecules, in order of that in the .xyz file
    #fil is the Tinker .xyz file for the desired output
    with open('%s'%(fil),'r') as l:
      cord = [lines.split() for lines in l]
      cord = np.array(list(it.izip_longest(*cord[2:],fillvalue=' '))).T
    mass = np.zeros(len(cord[:,0])/nmol)
    for i in np.arange(0,len(mass),1):
      num = int(cord[i,5])
      if num < 10:
        m = subprocess.check_output("less %s | grep 'atom          %s' | grep -oP '[0-9]*\.[0-9]*'"%(oplsaa,num),shell=True)
      elif 10 <= num < 100:
        m = subprocess.check_output("less %s | grep 'atom         %s' | grep -oP '[0-9]*\.[0-9]*'"%(oplsaa,num),shell=True)
      elif num >= 100:
        m = subprocess.check_output("less %s | grep 'atom        %s' | grep -oP '[0-9]*\.[0-9]*'"%(oplsaa,num),shell=True)
      mass[i] = m
    return mass

#################################################################################
#################          Thermo Properties          ###########################
#################################################################################

def vib_helm(T,wvn,mole):
    #Calculates the vibrational free energy
    #Output is quantum and classical vibrational free energy
    #T is the input temperature
    #wvn is the input wavenumbers
    #mole is the number of molecules in the lattice
    c = 2.998 * 10**10 #Speed of light in cm/s
    h = 2.520 * 10**(-35) #Reduced Plank's constant in cal*s
    k = 3.2998 * 10**(-24) #Boltzmann constant in cal*K
    Na = 6.022 * 10**(23) #Avogadro's number
    B = 1/(k*T)
    mol = float(mole)
    wvn = np.sort(wvn)
    fq = []
    fc = []
    for i in wvn[3:]:
      if i > 0:
        q = ((h*i*c*np.pi) + (1/B)* np.log(1 - np.exp(-B*h*i*c*2*np.pi)))*Na/1000
        fq.append(q)
        cl = (1/B)*np.log(B*h*i*c*2*np.pi)*Na/1000
        fc.append(cl)
      else:
        pass
    Fq = sum(fq)/mol
    Fc = sum(fc)/mol
    F =[Fq, Fc]
    return F

def vib_ent(T,wvn,mole):
    #Calculates the vibrational entropy
    #Output is quantum and classical vibrational free energy
    #T is the input temperature
    #wvn is the input wavenumbers
    #mole is the number of molecules in the lattice
    c = 2.998 * 10**10 #Speed of light in cm/s
    h = 2.520 * 10**(-35) #Reduced Plank's constant in cal*s
    k = 3.2998 * 10**(-24) #Boltzmann constant in cal*K
    Na = 6.022 * 10**(23) #Avogadro's number
    B = 1/(k*T)
    mol = float(mole)
    wvn = np.sort(wvn)
    sq = []
    sc = []
    for i in wvn[3:]:
      if i > 0:
        q = (h*i*c*2*np.pi/(T*(np.exp(B*h*i*c*2*np.pi)-1)) - k*np.log(1 - np.exp(-B*h*i*c*2*np.pi)))*Na/1000
        sq.append(q)
        cl = k*(1 - np.log(B*h*i*c*2*np.pi))*Na/1000
        sc.append(cl)
      else:
        pass
    Sq = sum(sq)/mol
    Sc = sum(sc)/mol
    S =np.array([Sq, Sc])
    return S

def vib_ener(T,wvn,mole):
    #Calculates the vibrational entropy
    #Output is quantum and classical vibrational free energy
    #T is the input temperature
    #wvn is the input wavenumbers
    #mole is the number of molecules in the lattice
    c = 2.998 * 10**10 #Speed of light in cm/s
    h = 2.520 * 10**(-35) #Reduced Plank's constant in cal*s
    k = 3.2998 * 10**(-24) #Boltzmann constant in cal*K
    Na = 6.022 * 10**(23) #Avogadro's number
    B = 1/(k*T)
    mol = float(mole)
    wvn = np.sort(wvn)
    sq = []
    sc = []
    for i in wvn[3:]:
      if i > 0:
        q = (0.5*h*i + 0.5*h*i/(np.exp(0.5*B*h*i)-1))*Na/1000
        sq.append(q)
        cl = (1/B)*Na/1000
        sc.append(cl)
      else:
        pass
    Uq = sum(sq)/mol
    Uc = sum(sc)/mol
    Uv =np.array([Uq, Uc])
    return Uv

def PV(fil,P):
    #Calculates the energy contributing to the PV term
    #fil is the Tinker .xyz file for the desired output
    #P is the desired pressure
    Na = 6.022 * 10**(23) #Avogadro's number
    convert = 2.390 * 10**(-29) #Converting bar*A^3 to kcal
    with open('%s'%(fil),'r') as l:
      vec = [lines.split() for lines in l]
      vec = np.array(list(it.izip_longest(*vec,fillvalue=' '))).T
      vec = vec[1,:6].astype(float)
    V = np.prod(vec[:3])*np.sqrt(1-np.cos(np.radians(vec[3]))**2 - np.cos(np.radians(vec[4]))**2 - np.cos(np.radians(vec[5]))**2 + 2*np.cos(np.radians(vec[3]))*np.cos(np.radians(vec[4]))*np.cos(np.radians(vec[5])))
    PV = P*V*convert*Na
    return V, PV, vec

def Vol(a,b,c,alph,bet,gam):
    V = a*b*c*np.sqrt(1-np.cos(np.radians(alph))**2 - np.cos(np.radians(bet))**2 - np.cos(np.radians(gam))**2 + 2*np.cos(np.radians(alph))*np.cos(np.radians(bet))*np.cos(np.radians(gam)))
    return V

def Gru_param(wvn,Vol,wvn2,Vol2):
    #Calcualtes all of the Gruneisen parameters for a particular system
    #wvn is the input wavenumbers of the reference structure
    #Vol is the input volume of the reference strucutre
    gru = -(np.log(wvn) - np.log(wvn2))/(np.log(Vol) - np.log(Vol2))
    gru[:3] = 0
    return gru
