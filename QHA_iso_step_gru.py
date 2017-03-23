#!/usr/bin/env python
###QHA based off of center of mass movement using the Gruneisen parameter

###Created by: Nathan Abraham
###Shirts Group
###University of Colorado, Boulder
###Department of Chemical and Biological Engineering

import os
import sys
import os.path
import numpy as np
import properties as prop
import QHAsub as Qs
import NumHes as NH
from optparse import OptionParser

###Input options
parser = OptionParser()
parser.add_option('-f', dest = 'fil', help = '.xyz file to be converted')
parser.add_option('-k', dest = 'key', help = 'key file for parameters')
parser.add_option('-n', dest = 'nmol', help = 'number of molecules in .xyz file')
parser.add_option('-p', dest = 'ply', help = 'polymorph characteristic marker')
       #Optional inputs#
parser.add_option('-P', dest = 'prs', help = 'Pressure of the system in atm', default = 1)
parser.add_option('-T', dest = 'temp', help = '.npy file containing desired temperature points', default = False)
parser.add_option('-H', dest = 'Hes', help = 'Analytical(a) or numerical(n) hessian',default = 'a')

(options, args) = parser.parse_args()
fil = options.fil
key = options.key
nmol = int(options.nmol)
ply = int(options.ply)
prs = float(options.prs)
temp = options.temp
Hes = options.Hes

###Checking if cordinate and wavenumber directories exist
if os.path.isdir('cords') != True:
  os.system('mkdir cords')
if os.path.isdir('wave') != True:
  os.system('mkdir wave')

###Temperature range
if temp == False:
  T = np.arange(25,301.,25)
  T = np.insert(T,0,0.1)
elif os.path.isfile(temp) == True:
  T = np.load(temp).astype(float)
elif len(temp.split(',')) >= 2:
  T = np.array(temp.split(',')).astype(float)
elif len(temp.split('X')) == 2:
  temp = np.array(temp.split('X')).astype(float)
  T = np.arange(0,temp[1],temp[0])
else:
  print "Not an appropriate temperature input"
  sys.exit()

###Vector fraction for second point to find the Gruneisen parameter
Vec_n = 1.0005

###Finding isotropically expanded strucutre based off of Vec_n
if os.path.isfile('cords/p%s_%s.xyz'%(ply,Vec_n)) != True:
  Qs.iso_expand(fil,key,Vec_n,1.0,nmol,ply)
else:
  #If the structure has already been created, moving it up to here
  os.system('cp cords/p%s_%s.xyz ./'%(ply,Vec_n))

fil2 = 'p%s_%s.xyz'%(ply,Vec_n)

#Setting the volumes of the two unit cells
Vol = prop.PV(fil,prs)[0]
Vol2 = prop.PV(fil2,1)[0]


###Finding the Gruneisen parameter and Reference wavenumbers
print "Finding the Gruneisen paramters"
#Analytically
if Hes == 'a':
  if os.path.isfile('wave/wave_QHA_p%s_v%s.npy'%(ply,1.0)) != True:
    wvn = prop.eig_wvn(fil,key)[1]
    np.save('wave/wave_QHA_p%s_v%s.npy'%(ply,1.0),wvn)
  else:
    wvn = np.load('wave/wave_QHA_p%s_v%s.npy'%(ply,1.0))
  if os.path.isfile('wave/wave_QHA_p%s_v%s.npy'%(ply,Vec_n)) != True:
    wvn2 = prop.eig_wvn(fil2,key)[1]
    np.save('wave/wave_QHA_p%s_v%s.npy'%(ply,Vec_n),wvn2)
  else:
    wvn2 = np.load('wave/wave_QHA_p%s_v%s.npy'%(ply,Vec_n))
  if (any(wvn) or any(wvn2)) < -1.0:
    print "Unacceptable negative wavenumbers detected"
    print wvn[:10]
    print wvn2[:10]
    sys.exit()
  gru = prop.Gru_param(wvn,Vol,wvn2,Vol2)

#Numerically
elif Hes == 'n':
  if os.path.isfile('wave/wave_N_QHA_p%s_v%s.npy'%(ply,1.0)) != True:
    wvn = NH.num_wvn(fil,key,nmol)
    np.save('wave/wave_N_QHA_p%s_v%s.npy'%(ply,1.0),wvn)
  else:
    wvn = np.load('wave/wave_N_QHA_p%s_v%s.npy'%(ply,1.0))
  if os.path.isfile('wave/wave_N_QHA_p%s_v%s.npy'%(ply,Vec_n)) != True:
    wvn2 = NH.num_wvn(fil2,key,nmol)
    np.save('wave/wave_N_QHA_p%s_v%s.npy'%(ply,Vec_n),wvn2)
  else:
    wvn2 = np.load('wave/wave_N_QHA_p%s_v%s.npy'%(ply,Vec_n))
  if (any(wvn) or any(wvn2)) < -1.0:
    print "Unacceptable negative wavenumbers detected"
    print wvn[:10]
    print wvn2[:10]
    sys.exit()
  gru = prop.Gru_param(wvn,Vol,wvn2,Vol2)

#Leaving if any other input
else:
  print "That is not a method type for Hessian, choose either Numerical(n) or Analytical(a)"
  sys.exit()

os.system('mv p%s_%s.xyz cords'%(ply,Vec_n))


###Calculating Thermal properties
Vec_s = 0.0005 #Step size for lattice vector multiplier
Vec_l = 0.97  #Lower bound of lattice vector multiplier
Vec_h = 1.04  #Uppder bound of lattice vector multiplier
Vec = np.arange(Vec_l,Vec_h,Vec_s)

print "Calculating thermal properties"

properR = np.zeros((len(Vec),len(T),8,2))
  #Quantum in column 0
  #Classical in column 1
  #Rows: 0 U, 1 Fv, 2 PV, 3 G, 4 Sv, 5 H, 6 V, 7 Uv

k = '%s'%(fil)
hold = 1.0
for j in np.arange(len(Vec)):
  if os.path.isfile('cords/p%s_%s.xyz'%(ply,Vec[j])) == True:
    os.system('cp cords/p%s_%s.xyz temp.xyz'%(ply,Vec[j]))
  else:
    Qs.iso_expand(k,key,Vec[j],hold,nmol,ply)
    os.system('cp p%s_%s.xyz temp.xyz'%(ply,Vec[j]))
    os.system('mv p%s_%s.xyz cords/'%(ply,Vec[j]))
  k = 'temp.xyz'
  properR[j,:,0,:] = prop.lat_ener('%s'%(k),'%s'%(key))/nmol
  Vol_new = Vol*(Vec[j])**(3)
  wvn_new = Qs.wvn_gru(wvn,Vol,Vol_new,gru)
  hold = Vec[j]
  for i in np.arange(len(T)):
    if T[i] != 0:
      properR[j,i,1,:] = prop.vib_helm(T[i],wvn_new,nmol)
      properR[j,i,4,:] = prop.vib_ent(T[i],wvn_new,nmol)
      properR[j,i,7,:] = prop.vib_ener(T[i],wvn,nmol)
    properR[j,i,2,:] = prop.PV('%s'%(fil),prs)[1]*Vec[j]**3/nmol
    properR[j,i,3,:] = np.sum(properR[j,i,:3,:],axis=0)
    properR[j,i,5,:] = properR[j,0,0,0] + properR[j,0,2,0] + properR[j,0,7,:]
    properR[j,i,6,:] = Vol_new

os.system('rm temp.xyz')

###Saving raw data
if Hes == 'a':
  np.save('Raw_gQ_%s_%s'%(ply,nmol),properR)
elif Hes == 'n':
  np.save('Raw_gNQ_%s_%s'%(ply,nmol),properR)

###Finding the minimum Gibbs free energy path
proper = np.zeros((len(T),8,2))

for i in np.arange(len(T)):
  for j in np.arange(len(Vec)):
    if properR[j,i,3,0] == min(properR[:,i,3,0]):
      proper[i,:,0] = properR[j,i,:,0]
    if properR[j,i,3,1] == min(properR[:,i,3,1]):
      proper[i,:,1] = properR[j,i,:,1]

#Saving minimum gibbs free energy data
if Hes == 'a':
#  np.save('UqgQ_%s_%s'%(ply,nmol),proper[:,0,0])
#  np.save('UcgQ_%s_%s'%(ply,nmol),proper[:,0,1])
#  np.save('FqgQ_%s_%s'%(ply,nmol),proper[:,1,0])
#  np.save('FcgQ_%s_%s'%(ply,nmol),proper[:,1,1])
#  np.save('GqgQ_%s_%s'%(ply,nmol),proper[:,3,0])
  np.save('GcgQ_%s_%s'%(ply,nmol),proper[:,3,1])
#  np.save('SqgQ_%s_%s'%(ply,nmol),proper[:,4,0])
#  np.save('ScgQ_%s_%s'%(ply,nmol),proper[:,4,1])
#  np.save('HqgQ_%s_%s'%(ply,nmol),proper[:,5,0])
#  np.save('HcgQ_%s_%s'%(ply,nmol),proper[:,5,1])
#  np.save('VqgQ_%s_%s'%(ply,nmol),proper[:,6,0])
  np.save('VcgQ_%s_%s'%(ply,nmol),proper[:,6,1])
  np.save('T_gQ',T)
elif Hes == 'n':
#  np.save('UqgNQ_%s_%s'%(ply,nmol),proper[:,0,0])
#  np.save('UcgNQ_%s_%s'%(ply,nmol),proper[:,0,1])
#  np.save('FqgNQ_%s_%s'%(ply,nmol),proper[:,1,0])
#  np.save('FcgNQ_%s_%s'%(ply,nmol),proper[:,1,1])
#  np.save('GqgNQ_%s_%s'%(ply,nmol),proper[:,3,0])
  np.save('GcgNQ_%s_%s'%(ply,nmol),proper[:,3,1])
#  np.save('SqgNQ_%s_%s'%(ply,nmol),proper[:,4,0])
#  np.save('ScgNQ_%s_%s'%(ply,nmol),proper[:,4,1])
#  np.save('HqgNQ_%s_%s'%(ply,nmol),proper[:,5,0])
#  np.save('HcgNQ_%s_%s'%(ply,nmol),proper[:,5,1])
#  np.save('VqgNQ_%s_%s'%(ply,nmol),proper[:,6,0])
  np.save('VcgNQ_%s_%s'%(ply,nmol),proper[:,6,1])
  np.save('T_gNQ',T)

