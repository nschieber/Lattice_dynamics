#!/usr/bin/env python
###Harmonic approximation

###Created by: Nathan Abraham
###Shirts Group
###University of Colorado, Boulder
###Department of Chemical and Biological Engineering

###Update on April 29, 2016
###Update on August 12, 2016 - Additon of mode movements for negative wavenumbers
###Update on August 31, 2016 - Removed mode movements for negative wavenumbers added in entropy  terms
###Update on September 20, 2016 - Changed some input options and modified thermal properties ###Update on October 27, 2016 - Increased efficiency and added in enthalpy calculation

import os
import sys
import numpy as np
import properties as prop
import NumHes as NH
from optparse import OptionParser
import time
start = time.time()

###Input options
parser = OptionParser()
parser.add_option('-f', dest = 'fil', help = '.xyz file to be converted')
parser.add_option('-k', dest = 'key', help = 'key file for parameters')
parser.add_option('-n', dest = 'nmol', help = 'number of molecules in .xyz file')
parser.add_option('-p', dest = 'ply', help = 'polymorph characteristic')
       #Optional inputs#
parser.add_option('-P', dest = 'prs', help = 'Pressure of the system in atm', default = 1)
parser.add_option('-T', dest = 'temp', help = '.npy file containing desired temperature points', default = False)
parser.add_option('-H', dest = 'Hes', help = 'Analytical(a) or numerical(n) hessian',default = 'a')
parser.add_option('-s', dest = 'step', help = 'Analytical(a) or numerical(n) hessian',default = 0.00001)

(options, args) = parser.parse_args()
fil = options.fil
key = options.key
nmol = int(options.nmol)
ply = int(options.ply)
prs = float(options.prs)
temp = options.temp
Hes = options.Hes
step = float(options.step)

###Temperature range
if temp == False:
  T = np.arange(10,301,10)
  T = np.insert(T,0,0)
else:
  T = np.load(temp)

###Seeing if the wavenumbers have been calculated already
if Hes == 'a':
  if os.path.isfile('wave_HA_%s.npy'%(ply)) == True:
    #If so, unpacking the pre-existing wavenumbers
    wave = np.load('wave_HA_%s.npy'%(ply))
  else:
    #If not, calculation of the wavenumbers
    wave = prop.eig_wvn(fil,key)[1]
    np.save('wave_HA_%s'%(ply),wave)

if Hes == 'n':
  if os.path.isfile('wave_N_HA_%s.npy'%(ply)) == True:
    #If so, unpacking the pre-existing wavenumbers
    wave = np.load('wave_N_HA_%s.npy'%(ply))
  else:
    #If not, calculation of the wavenumbers
    wave = NH.num_wvn(fil,key,nmol,step)
    np.save('wave_N_HA_%s'%(ply),wave)

###Ending run if any negative wavenumbers exist
if any(wave < -1.0) == True:
  print "Structure is not properly minimized"
  os.system("rm hessian.txt analyze.txt")
  sys.exit()

###Calculating Thermal properties
proper = np.zeros((len(T),6,2))
  #Quantum in column 0
  #Classical in column 1
  #Rows: 0 U, 1 Fv, 2 PV, 3 G, 4 Sv, 5 H

proper[:,0,:] = prop.lat_ener('%s'%(fil),'%s'%(key))/nmol

for i in np.arange(len(T)):
  if T[i] != 0:
    proper[i,1,:] = prop.vib_helm(T[i],wave,nmol)
    proper[i,4,:] = prop.vib_ent(T[i],wave,nmol)
  proper[i,2,:] = prop.PV('%s'%(fil),prs)[1]/nmol
  proper[i,3,:] = np.sum(proper[i,:3,:],axis=0)
  proper[i,5,:] = proper[0,0,0] + proper[0,2,0]

end = time.time()

##Saving files
if Hes == 'a':
  #np.save('U_HA_%s_%s'%(ply,nmol),proper[:,0,0])
  #np.save('H_HA_%s_%s'%(ply,nmol),proper[:,5,0])
  #np.save('FqHA_%s_%s'%(ply,nmol),proper[:,1,0])
  #np.save('FcHA_%s_%s'%(ply,nmol),proper[:,1,1])
  #np.save('SqHA_%s_%s'%(ply,nmol),proper[:,4,0])
  #np.save('ScHA_%s_%s'%(ply,nmol),proper[:,4,1])
  #np.save('GqHA_%s_%s'%(ply,nmol),proper[:,3,0])
  np.save('GcHA_%s_%s'%(ply,nmol),proper[:,3,1])
  np.save('T_HA',T)
  #np.save('Time_HA_%s'%(ply),end-start)
if Hes == 'n':
  #np.save('U_NHA_%s_%s'%(ply,nmol),proper[:,0,0])
  #np.save('H_NHA_%s_%s'%(ply,nmol),proper[:,5,0])
  #np.save('FqNHA_%s_%s'%(ply,nmol),proper[:,1,0])
  #np.save('FcNHA_%s_%s'%(ply,nmol),proper[:,1,1])
  #np.save('SqNHA_%s_%s'%(ply,nmol),proper[:,4,0])
  #np.save('ScNHA_%s_%s'%(ply,nmol),proper[:,4,1])
  #np.save('GqNHA_%s_%s'%(ply,nmol),proper[:,3,0])
  np.save('GcNHA_%s_%s'%(ply,nmol),proper[:,3,1])
  np.save('T_NHA',T)
  #np.save('TimeN_HA_%s'%(ply),end-start)
