#!/usr/bin/env python
###Isotropic QHA based off of center of mass movement using the a gradient approach to determine expansion

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
from scipy.optimize import fsolve
from optparse import OptionParser

###Input options
parser = OptionParser()
parser.add_option('-f', dest = 'fil', help = '.xyz file to be converted')
parser.add_option('-k', dest = 'key', help = 'key file for parameters')
parser.add_option('-n', dest = 'nmol', help = 'number of molecules in .xyz file')
parser.add_option('-p', dest = 'ply', help = 'polymorph characteristic')
       #Optional inputs#
parser.add_option('-P', dest = 'prs', help = 'Pressure of the system in atm', default = 1)
parser.add_option('-T', dest = 'temp', help = '.npy file containing desired temperature points', default = False)
parser.add_option('-R', dest = 'RK', help ='Runge-Kuta step size in [K]',default =150.)
parser.add_option('-r', dest = 'rev', help='Input -1 to go in reverse from 300K', default=0)

(options, args) = parser.parse_args()
fil = options.fil
key = options.key
nmol = int(options.nmol)
ply = int(options.ply)
prs = float(options.prs)
temp = options.temp
RK = float(options.RK)
rev = int(options.rev)


###Temperature range
T_ex = np.arange(RK,301,RK)
T_ex = np.insert(T_ex,0,0.1)
if rev == -1:
  print "Performing gradient search in reverse from 300K point"
  T_ex = T_ex[::-1]
  print T_ex

print "Searching for expanded and compressed structure"
#Seeing if a directory is made to save expanded strucutres
if os.path.isdir('cords') != True:
  os.system('mkdir cords')
if os.path.isdir('wave') != True:
  os.system('mkdir wave')

os.system('cp %s temp.xyz'%(fil))
path = 'R%s-%s'%(RK,max(T_ex))

Dv = np.zeros(len(T_ex)-1)
if os.path.isfile('dV_p%s_%s.npy'%(ply,path)) == True:
  print "Using dV/dT data from previous run"
  Dv = np.load('dV_p%s_%s.npy'%(ply,path))

proper = np.zeros((len(T_ex),13,2))
  #Quantum in column 0
  #Classical in column 1
  #Rows: 0 U, 1 Fv, 2 PV, 3 G, 4 Sv, 5 H, 6 V,
  #      7 a-vec, 8 b-vec, 9 c-vec, 10 alpha angle
  #      11 beta angle, 12 gamma angle

for i in range(len(T_ex)):
  PV = prop.PV('temp.xyz',prs)
  if T_ex[i] != T_ex[-1]:
    if Dv[i] == 0.0:
      K = np.zeros(4)
      h = T_ex[i+1] - T_ex[i]
      t = np.array([0,h/2,h/2,h]).astype(float)
      RK_cord = 'RK_temp.xyz'
      for j in range(len(K)):
        if j != 0:
          RK_i = (((K[j-1]*t[j] + PV[0])/PV[0])**(1/3.) - 1.)*PV[2]
          Qs.expand('temp.xyz',key,RK_i[0],RK_i[1],RK_i[2],0.0,0.0,0.0,nmol,ply)
          os.system('mv p%s_temp.xyz %s'%(ply,RK_cord))
        else:
          os.system('cp temp.xyz %s'%(RK_cord))
        os.system('cp %s p%s_%s_T%s.xyz'%(RK_cord,ply,path+'_k%s-%s'%(j+1,T_ex[i]),T_ex[i]+t[j]))
        K[j] = Qs.iso_grad(RK_cord,key,prs,ply,nmol,T_ex[i]+t[j],path+'_k%s-%s'%(j+1,T_ex[i]))

        if j == 0:
          os.system('cp iso_p%s_%s.npy wave/iso_p%s_%s_T%s.npy'%(ply,path+'_k%s-%s_T%s'%(j+1,T_ex[i],T_ex[i]+t[j]),ply,path,T_ex[i]))        
        files = [f for f in os.listdir('./') if f.startswith("p%s_%s"%(ply,path+'_k%s-%s'%(j+1,T_ex[i])))]
        wiles = [f for f in os.listdir('./') if f.startswith("iso_p%s_%s"%(ply,path+'_k%s-%s'%(j+1,T_ex[i])))]
        for k in range(len(files)):
          os.system('rm %s'%(files[k]))
          os.system('rm %s' %(wiles[k]))
    
      Dv[i] = np.sum(np.multiply(K,np.array([1./6,1./3,1./3,1./6])))
      np.save('dV_p%s_%s'%(ply,path),Dv)
  else:
    np.save('wave/iso_p%s_%s_T%s.npy'%(ply,path,T_ex[i]),prop.eig_wvn('temp.xyz',key)[1])
  wvn = np.load('wave/iso_p%s_%s_T%s.npy'%(ply,path,T_ex[i]))

  if any(wvn < -1.0) == True:
    proper[:,:,:] = np.nan
  else:
    proper[:,0,:] = prop.lat_ener('temp.xyz',key)/nmol
    if T_ex[i] != 0:
      proper[i,1,:] = prop.vib_helm(T_ex[i],wvn,nmol)
      proper[i,4,:] = prop.vib_ent(T_ex[i],wvn,nmol)
      proper[i,7,:] = prop.vib_ener(T_ex[i],wvn,nmol)
    proper[i,2,:] = PV[1]/nmol
    proper[i,3,:] = np.sum(proper[i,:3,:],axis=0)
    proper[i,5,:] = proper[0,0,0] + proper[0,2,0] + proper[0,7,:]
    proper[i,6,:] = PV[0]
    proper[i,7:,0] = PV[2]
    proper[i,7:,1] = PV[2]

  os.system('cp temp.xyz cords/p%s_%s_T%s.xyz'%(ply,path,T_ex[i]))

  if T_ex[i] != T_ex[-1]:
    if os.path.isfile('cords/p%s_%s_T%s.xyz'%(ply,path,T_ex[i+1])) != True:
      i_T = (((Dv[i]*(T_ex[i+1]-T_ex[i]) + PV[0])/PV[0])**(1/3.) - 1.)*PV[2]
      Qs.expand('temp.xyz',key,i_T[0],i_T[1],i_T[2],0.0,0.0,0.0,nmol,ply)
      os.system('mv p%s_temp.xyz temp.xyz'%(ply))
    else:
      os.system('cp cords/p%s_%s_T%s.xyz ./temp.xyz'%(ply,path,T_ex[i+1]))

os.system('rm RK_temp.xyz temp.xyz')

#np.save('Uq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,0,0])
#np.save('Uc_Iso_Grd_%s_%s'%(ply,nmol),proper[:,0,1])
#np.save('Fq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,1,0])
#np.save('Fc_Iso_Grd_%s_%s'%(ply,nmol),proper[:,1,1])
#np.save('Gq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,3,0])
np.save('Gc_Iso_Grd_%s'%(path),proper[:,3,1])
#np.save('Sq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,4,0])
#np.save('Sc_Iso_Grd_%s_%s'%(ply,nmol),proper[:,4,1])
#np.save('Hq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,5,0])
#np.save('Hc_Iso_Grd_%s_%s'%(ply,nmol),proper[:,5,1])
#np.save('Vq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,6,0])
np.save('Vc_Iso_Grd_%s'%(path),proper[:,6,1])
#np.save('LATq_Iso_Grd_%s_%s'%(ply,nmol),proper[:,7:,0])
#np.save('LATc_Iso_Grd_%s_%s'%(ply,nmol),proper[:,7:,1])
np.save('T_Iso_Grd_%s'%(path),T_ex)

print "Done"
