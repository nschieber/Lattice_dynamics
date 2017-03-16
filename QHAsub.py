# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:35:38 2016

@author: nabraham
"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import itertools as it
import properties as prop
#######################################################################################################
#######################################################################################################
def iso_expand(fil,key,vec_n,vec_o,nmol,ply):
    #Expands or compresses the 
    #fil is the Tinker .xyz file for the desired output
    #key is the corresponding Tinker .key file for
    #frac is the fraction to change lattice vectors by
    #nmol is the number of molecules per cell specified
    with open('%s'%(fil)) as f:
        z = [lines.split() for lines in f]
        dummy = z[0:2]
        nstrc = z[0:2]
        vec = np.array(dummy[1]).astype(float)
        n = len(z[2:])/nmol       ###Number of molecules
        afile = np.array(list(it.izip_longest(*z[2:], fillvalue=' '))).T
        #This next part depends on atom connections, I couldn't think of a more efficient way to do this
        if len(afile[2]) == 10:
            l = np.column_stack((afile[0:,0],afile[0:,1],afile[0:,2].astype(float),afile[0:,3].astype(float),afile[0:,4].astype(float),afile[0:,5],afile[0:,6],afile[0:,7],afile[0:,8],afile[0:,9]))
        if len(afile[2]) == 9:
            l = np.column_stack((afile[0:,0],afile[0:,1],afile[0:,2].astype(float),afile[0:,3].astype(float),afile[0:,4].astype(float),afile[0:,5],afile[0:,6],afile[0:,7],afile[0:,8]))
        if len(afile[2]) == 8:
            l = np.column_stack((afile[0:,0],afile[0:,1],afile[0:,2].astype(float),afile[0:,3].astype(float),afile[0:,4].astype(float),afile[0:,5],afile[0:,6],afile[0:,7]))
    dummy[0] = str('    %s  DUMMY SOLID T=   0.00000\n' %(nmol))
    dummy[1] = str('    %s    %s    %s    %s    %s    %s\n'%(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5]))
    for j in np.arange(0,nmol,1):
        a = j + 1
        D = np.mean(l[j*n:(a)*n,2:5].astype(float),axis = 0)
        d = ('   %s    DM   %s    %s    %s    57\n' %(a,D[0],D[1],D[2]))
        dummy.append(d)
        l[j*n:(a)*n,2:5] = np.subtract(l[j*n:(a)*n,2:5].astype(float),D)
    with open('dummy.xyz', 'w') as file:
        file.writelines(dummy)
    ###Dummy file to fractional form
    os.system("crystal dummy.xyz -k %s 2" %(key))
    ###Adjusting box vectors of dummy file
    with open('dummy.frac','r') as frac:
            dummy_f = frac.readlines()
    vec[0:3] = vec[0:3]*vec_n/vec_o
    vecln = str('     %s    %s    %s   %s   %s   %s  \n' %(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5]))
    dummy_f[1] = vecln
    with open('dummy.frac_temp','w') as file:
        file.writelines(dummy_f)
    ###Dummy file back to cartesian coordinates
    os.system("crystal dummy.frac_temp -k %s 1" %(key))
    lfile = pd.read_table('dummy.xyz_2',skiprows=2,header = None, sep=r"\s*", engine='python')
    lfile = np.column_stack((lfile[2],lfile[3],lfile[4]))
    ###Adjustment of new molecular coordinates
    for j in np.arange(0,nmol,1):
        l[j*n:(j+1)*n,2:5] = np.subtract(l[j*n:(j+1)*n,2:5].astype(float),-lfile[j])
    ###Minimization of new structure
    nstrc[1] = dummy_f[1] +'\n'
    for j in l:
        ns = '   '
        for k in j:
            ns = ns + '%s    ' %(k)
        ns = ns + '\n'
        nstrc.append(ns)
    nst0 = '    '
    for j in nstrc[0]:
        nst0 = nst0 + '%s  '%(j)
    nstrc[0] = nst0 + '\n'
    with open('%s_2'%(fil),'w') as fil10:
        fil10.writelines(nstrc)
    os.system('minimize %s_2 -k %s 0.01' %(fil, key))
    os.system('mv %s_3 p%s_%s.xyz'%(fil,ply,vec_n))
    os.system('rm dummy.frac_temp dummy.xyz_2 dummy.frac dummy.xyz %s_2' %(fil))

#######################################################################################################
#######################################################################################################

def wvn_gru(wvn_ref,V_ref,V_new,Gru):
    #Calcuates the new wavenumbers from the Grunseian parameter
    fact = np.diag(np.power(V_new/V_ref,-1*Gru))
    wvn_new = np.dot(wvn_ref,fact)
    return wvn_new

#######################################################################################################
#######################################################################################################

def expand(fil,key,a,b,c,alp,bet,gam,nmol,ply):
    #Finds an expaned or compressed strucutre with specific changes to each lattice parameter
    #fil initial input strucutre
    #key keyfile to use

    ##All lengths are in Ang.
    #a change in first lattice vector
    #b change in seccon lattice vector
    #c change in third lattice vector

    ##All angles are in deg.
    #alp change in first lattice angle
    #bet change in second lattice angle
    #gam change in third lattice angle

    #nmol number of molecules in fil
    #ply is polymorph marker
    import itertools as it
    import numpy as np
    np.set_printoptions(suppress=True)

    delt = np.array([a,b,c,alp,bet,gam])

    with open(fil) as f:
      afile = np.array(list(it.izip_longest(*[lines.split() for lines in f], fillvalue=' '))).T
    cords = afile[2:,2:5].astype(float)
    vec = afile[1,:6].astype(float)
    atom = len(cords[:,0])/nmol
    
    dummy = ['%s     \n'%(nmol),'  %s  %s  %s  %s  %s  %s\n'%(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5])]
    for i in range(nmol):
        a = i + 1
        D = np.mean(cords[i*atom:(a)*atom],axis = 0)
        d = ('   %s    DM   %s    %s    %s    57\n' %(a,D[0],D[1],D[2]))
        dummy.append(d)
        cords[i*atom:(a)*atom] = np.subtract(cords[i*atom:(a)*atom],D)
    with open('dummy.xyz', 'w') as file:
        file.writelines(dummy)
    
    os.system("crystal dummy.xyz -k %s 2 &> /dev/null" %(key))

    with open('dummy.frac','r') as frac:
            dummy_f = frac.readlines()
    vec = vec + delt
    dummy_f[1] = str('     %s    %s    %s   %s   %s   %s  \n' %(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5]))
    with open('dummy.frac_temp','w') as file:
        file.writelines(dummy_f)

    os.system("crystal dummy.frac_temp -k %s 1 &> /dev/null" %(key))
    
    with open('dummy.xyz_2') as f:
      dummy = np.array(list(it.izip_longest(*[lines.split() for lines in f], fillvalue=' '))).T
      D = dummy[2:,2:5].astype(float)
    for i in range(nmol):
        a = i + 1
        cords[i*atom:(a)*atom] = np.subtract(cords[i*atom:(a)*atom],-D[i])
    afile[2:,2:5] = np.around(cords,decimals=8).astype(str)
    afile[1,:6] = vec
    
    outfile = [] 
    for i in afile:
      out_str = ''
      for j in i:
        out_str = out_str + '   %s'%(j)
      out_str = out_str + '\n'
      outfile.append(out_str)

    with open('%s_2'%(fil),'w') as fil_out:
        fil_out.writelines(outfile)
    os.system('minimize %s_2 -k %s 0.01 &> /dev/null' %(fil, key))
    os.system('mv %s_3 p%s_temp.xyz'%(fil,ply))
    os.system('rm dummy.frac_temp dummy.xyz_2 dummy.frac dummy.xyz %s_2' %(fil))

#######################################################################################################
#######################################################################################################

def vec_gradient (fil_0,wvn_0,fil_1,wvn_1,fil_2,wvn_2,key,T_0,dT,P,nmol,ind):
    PV_0 = prop.PV(fil_0,P)
    PV_1 = prop.PV(fil_1,P)
    PV_2 = prop.PV(fil_2,P)

    U_0 = prop.lat_ener(fil_0,key)/nmol
    U_1 = prop.lat_ener(fil_1,key)/nmol
    U_2 = prop.lat_ener(fil_2,key)/nmol

    G_T0_V0 = PV_0[1]/nmol + U_0 + prop.vib_helm(T_0,wvn_0,nmol)[1]
    G_T1_V0 = PV_0[1]/nmol + U_0 + prop.vib_helm(T_0 + dT,wvn_0,nmol)[1]
    G_T2_V0 = PV_0[1]/nmol + U_0 + prop.vib_helm(T_0 + 2*dT,wvn_0,nmol)[1]

    G_T1_V1 = PV_1[1]/nmol + U_1 + prop.vib_helm(T_0 + dT,wvn_1,nmol)[1]

    G_T0_V2 = PV_2[1]/nmol + U_2 + prop.vib_helm(T_0,wvn_2,nmol)[1]
    G_T1_V2 = PV_2[1]/nmol + U_2 + prop.vib_helm(T_0 + dT,wvn_2,nmol)[1]
    G_T2_V2 = PV_2[1]/nmol + U_2 + prop.vib_helm(T_0 + 2*dT,wvn_2,nmol)[1]

    ddG_dTda = (G_T2_V2 - G_T2_V0 - G_T0_V2 + G_T0_V0)/(4*dT*(PV_1[2][ind]-PV_0[2][ind]))
    ddG_dda = (G_T1_V2 - 2*G_T1_V1 + G_T1_V0)/(PV_1[2][ind]-PV_0[2][ind])**2
    da_dT = - ddG_dTda/ddG_dda
    return da_dT

#######################################################################################################
#######################################################################################################

def vol_gradient (fil_0,wvn_0,fil_1,wvn_1,fil_2,wvn_2,key,T_0,dT,P,nmol):
    PV_0 = prop.PV(fil_0,P)
    PV_1 = prop.PV(fil_1,P)
    PV_2 = prop.PV(fil_2,P)

    U_0 = prop.lat_ener(fil_0,key)/nmol
    U_1 = prop.lat_ener(fil_1,key)/nmol
    U_2 = prop.lat_ener(fil_2,key)/nmol

    G_T0_V0 = PV_0[1]/nmol + U_0 + prop.vib_helm(T_0,wvn_0,nmol)[1]
    G_T1_V0 = PV_0[1]/nmol + U_0 + prop.vib_helm(T_0 + dT,wvn_0,nmol)[1]
    G_T2_V0 = PV_0[1]/nmol + U_0 + prop.vib_helm(T_0 + 2*dT,wvn_0,nmol)[1]

    G_T1_V1 = PV_1[1]/nmol + U_1 + prop.vib_helm(T_0 + dT,wvn_1,nmol)[1]

    G_T0_V2 = PV_2[1]/nmol + U_2 + prop.vib_helm(T_0,wvn_2,nmol)[1]
    G_T1_V2 = PV_2[1]/nmol + U_2 + prop.vib_helm(T_0 + dT,wvn_2,nmol)[1]
    G_T2_V2 = PV_2[1]/nmol + U_2 + prop.vib_helm(T_0 + 2*dT,wvn_2,nmol)[1]

    ddG_dTdV = (G_T2_V2 - G_T2_V0 - G_T0_V2 + G_T0_V0)/(4*dT*(PV_1[0]-PV_0[0]))
    ddG_ddV = (G_T1_V2 - 2*G_T1_V1 + G_T1_V0)/(PV_1[0]-PV_0[0])**2
    dV_dT = - ddG_dTdV/ddG_ddV
    return dV_dT

#######################################################################################################
#######################################################################################################

def aniso_grad (fil,key,P,ply,nmol,T_0,path):
    dT = 0.1
    Dl = prop.PV(fil,P)[2]*0.00005
    Dl[3:] = 0.0001
    mark = ['a','b','c','alp','bet','gam']
    wvn_1 = prop.eig_wvn(fil,key)[1]

    da = np.zeros((1,6))[0]
    for i in range(6):
      dl = np.zeros((1,6))[0]
      dl[i] = Dl[i]

      expand(fil,key,dl[0],dl[1],dl[2],dl[3],dl[4],dl[5],nmol,ply)
      wvn_2 = prop.eig_wvn('p%s_temp.xyz'%(ply),key)[1]
      fil_2 = 'p%s_ex%s_%s.xyz'%(ply,path,mark[i])
      os.system('mv p%s_temp.xyz %s'%(ply,fil_2))

      expand(fil,key,-dl[0],-dl[1],-dl[2],-dl[3],-dl[4],-dl[5],nmol,ply)
      wvn_0 = prop.eig_wvn('p%s_temp.xyz'%(ply),key)[1]
      fil_0 = 'p%s_ex%s_-%s.xyz'%(ply,path,mark[i])
      os.system('mv p%s_temp.xyz %s'%(ply,fil_0))

      da[i] = vec_gradient(fil_0,wvn_0,fil,wvn_1,fil_2,wvn_2,key,T_0,dT,P,nmol,i)
    return da

#######################################################################################################
#######################################################################################################

def iso_grad (fil,key,P,ply,nmol,T_0,path):
    dT = 0.01
    dl = prop.PV(fil,P)[2]*0.01

    wvn_1 = prop.eig_wvn(fil,key)[1]
    np.save('iso_p%s_%s_T%s.npy'%(ply,path,T_0),wvn_1)

    expand(fil,key,dl[0],dl[1],dl[2],0.0,0.0,0.0,nmol,ply)
    wvn_2 = prop.eig_wvn('p%s_temp.xyz'%(ply),key)[1]
    fil_2 = 'p%s_%s_T%s+.xyz'%(ply,path,T_0)
    os.system('mv p%s_temp.xyz %s'%(ply,fil_2))
    np.save('iso_p%s_%s_T%s+.npy'%(ply,path,T_0),wvn_2)

    expand(fil,key,-dl[0],-dl[1],-dl[2],0.0,0.0,0.0,nmol,ply)
    wvn_0 = prop.eig_wvn('p%s_temp.xyz'%(ply),key)[1]
    fil_0 = 'p%s_%s_T%s-.xyz'%(ply,path,T_0)
    os.system('mv p%s_temp.xyz %s'%(ply,fil_0))
    np.save('iso_p%s_%s_T%s-.npy'%(ply,path,T_0),wvn_0)
    dV = vol_gradient(fil_0,wvn_0,fil,wvn_1,fil_2,wvn_2,key,T_0,dT,P,nmol)
    return dV

#######################################################################################################
#######################################################################################################

def iso_grad_gru (fil,key,P,ply,nmol,T_0,path,gru,wvn_ref,V_ref):
    dT = 0.01
    PV = prop.PV(fil,P)
    dl = PV[2]*0.01

    wvn_1 = wvn_gru(wvn_ref,V_ref,PV[0],gru)

    expand(fil,key,dl[0],dl[1],dl[2],0.0,0.0,0.0,nmol,ply)
    wvn_2 = wvn_gru(wvn_ref,V_ref,PV[0]*(1.01)**3,gru)
    fil_2 = 'p%s_%s_T%s+.xyz'%(ply,path,T_0)
    os.system('mv p%s_temp.xyz %s'%(ply,fil_2))

    expand(fil,key,-dl[0],-dl[1],-dl[2],0.0,0.0,0.0,nmol,ply)
    wvn_0 = wvn_gru(wvn_ref,V_ref,PV[0]*(0.99)**3,gru)
    fil_0 = 'p%s_%s_T%s-.xyz'%(ply,path,T_0)
    os.system('mv p%s_temp.xyz %s'%(ply,fil_0))
    dV = vol_gradient(fil_0,wvn_0,fil,wvn_1,fil_2,wvn_2,key,T_0,dT,P,nmol)
    return dV

