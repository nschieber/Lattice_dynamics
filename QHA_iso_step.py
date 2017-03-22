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

(options, args) = parser.parse_args()
fil = options.fil
key = options.key
nmol = int(options.nmol)
ply = int(options.ply)
prs = float(options.prs)
temp = options.temp

###Temperature range
if temp == False:
  T = np.arange(5,301.,5)
  T = np.insert(T,0,0.1)
else:
  T = np.load(temp)

###Creating expanded and compressed structures
Vec_s = 0.0005
Vec_l = 0.99
Vec_h = 1.04
Vec = np.arange(Vec_l,Vec_h,Vec_s)

print "Searching for expanded and compressed structure"

#Seeing if a directory is made to save expanded strucutres
if os.path.isdir('cords') != True:
  os.system('mkdir cords')
  Qs.iso_expand('%s'%(fil),key,1.0,1.0,nmol,ply)
  os.system('mv p%s_1.0.xyz cords/'%(ply))
if os.path.isdir('wave') != True:
  os.system('mkdir wave')


hold = 1.0
for j in np.append(Vec[:(1.0-Vec_l)/(Vec_s):][::-1],Vec[(1.0-Vec_l)/(Vec_s):]):
  if os.path.isfile('cords/p%s_%s.xyz'%(ply,j)) != True:
    os.system('cp cords/p%s_%s.xyz ./temp.xyz'%(ply,hold))
    Qs.iso_expand('temp.xyz',key,j,hold,nmol,ply)
    if j == Vec_l:
      hold == 1.0
    else:
      hold == j
  else:
    os.system('cp cords/p%s_%s.xyz ./'%(ply,j))
  if os.path.isfile('wave/wave_QHA_p%s_v%s.npy'%(ply,j)) != True:
    wvn = prop.eig_wvn('p%s_%s.xyz'%(ply,j),'%s'%(key))[1]
    np.save('wave_QHA_p%s_v%s'%(ply,j),wvn)
    os.system('mv wave_QHA_p%s_v%s.npy wave/'%(ply,j))
  os.system('mv p%s_%s.xyz cords/'%(ply,j))
os.system('rm temp.xyz')

###Calculating Thermal properties
properR = np.zeros((len(Vec),len(T),8,2))
  #Quantum in column 0
  #Classical in column 1
  #Rows: 0 U, 1 Fv, 2 PV, 3 G, 4 Sv, 5 H, 6 V, 7 Uv

print "Calculating thermal properties"
for j in np.arange(len(Vec)):
  wvn = np.load('wave/wave_QHA_p%s_v%s.npy'%(ply,Vec[j]))
  if any(wvn < -1.0) == True:
    properR[j,:,:,:] = np.nan
  else:
    properR[j,:,0,:] = prop.lat_ener('cords/p%s_%s.xyz'%(ply,Vec[j]),'%s'%(key))/nmol
    for i in np.arange(len(T)):
      if T[i] != 0:
        properR[j,i,1,:] = prop.vib_helm(T[i],wvn,nmol)
        properR[j,i,4,:] = prop.vib_ent(T[i],wvn,nmol)
        properR[j,i,7,:] = prop.vib_ener(T[i],wvn,nmol)
      pv = prop.PV('%s'%(fil),prs)
      properR[j,i,2,:] = pv[1]*Vec[j]**3/nmol
      properR[j,i,3,:] = np.sum(properR[j,i,:3,:],axis=0)
      properR[j,i,5,:] = properR[j,0,0,0] + properR[j,0,2,0] + properR[j,0,7,:]
      properR[j,i,6,:] = pv[0]*Vec[j]**3

###Saving raw data
np.save('Raw_QHA_%s_%s'%(ply,nmol),properR)

###Finding the minimum Gibbs free energy path
proper = np.zeros((len(T),8,2))

for i in np.arange(len(T)):
  for j in np.arange(len(Vec)):
    if properR[j,i,3,0] == np.nanmin(properR[:,i,3,0]):
      proper[i,:,0] = properR[j,i,:,0]
    if properR[j,i,3,1] == np.nanmin(properR[:,i,3,1]):
      proper[i,:,1] = properR[j,i,:,1]

end = time.time()

#np.save('UqQHA_%s_%s'%(ply,nmol),proper[:,0,0])
#np.save('UcQHA_%s_%s'%(ply,nmol),proper[:,0,1])
#np.save('FqQHA_%s_%s'%(ply,nmol),proper[:,1,0])
#np.save('FcQHA_%s_%s'%(ply,nmol),proper[:,1,1])
#np.save('GqQHA_%s_%s'%(ply,nmol),proper[:,3,0])
np.save('GcQHA_%s_%s'%(ply,nmol),proper[:,3,1])
#np.save('SqQHA_%s_%s'%(ply,nmol),proper[:,4,0])
#np.save('ScQHA_%s_%s'%(ply,nmol),proper[:,4,1])
#np.save('HqQHA_%s_%s'%(ply,nmol),proper[:,5,0])
#np.save('HcQHA_%s_%s'%(ply,nmol),proper[:,5,1])
#np.save('VqQHA_%s_%s'%(ply,nmol),proper[:,6,0])
np.save('VcQHA_%s_%s'%(ply,nmol),proper[:,6,1])
np.save('T_QHA',T)
#np.save('Time_iQHA_%s'%(ply),end-start)
