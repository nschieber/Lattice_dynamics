# -*- coding: utf-8 -*-
"""
Created on Mon Nov 7, 2016

@author: nabraham
"""

import os
import sys
import subprocess
import numpy as np
import itertools as it
import properties as prop

def Hes_analyt(fil,key):
    #Returns the Analytical matrix from Tinker testhess run routine
    fname = os.path.splitext(os.path.basename('%s'%(fil)))[0]
    if os.path.isfile('%s.hes'%(fname)) != True:
      os.system('testhess %s -k %s Y N'%(fil,key))
    with open(fil,'r') as l:
      line = [lines.split() for lines in l]
      line = np.array(list(it.izip_longest(*line[:2],fillvalue=' '))).T
      atm = int(line[0,0])
    hold = []
    with open(fname +'.hes','r') as l:
      line = [lines.split() for lines in l]
      line = np.array(list(it.izip_longest(*line,fillvalue=' '))).T
    for i in line:
      try:
        for j in i:
          hold.append(float(j))
      except ValueError:
        pass
    hess_a = np.diag(hold[0:3*atm])
    hess_a[np.triu_indices(3*atm,1)] = hold[3*atm:]
    for i in np.arange(0,3*atm,1):
      for j in np.arange(0,3*atm,1):
        if i > j:
          hess_a[i,j] = hess_a[j,i]
    os.system('rm %s.hes'%(fname))
    return hess_a

def Hes_numer(fil,ana):
    #Returns the matrix of the comparison of Analytical and Numberical Hessian from Tinker testhess
    hess_n = ana
    with open('%s'%(fil),'r') as l:
      line = [lines.split() for lines in l]
      start = int(subprocess.check_output("sed -n '/1st Atom/=' TH_G.txt",shell=True))+1
      line = np.array(list(it.izip_longest(*line[start:],fillvalue=' '))).T
    for i in line:
      if i[0] == ' ':
        break
      else:
        for j in np.arange(0,4,1):
          if i[j] == '(X)':
            i[j] = 0
          elif i[j] == '(Y)':
            i[j] = 1
          elif i[j] == '(Z)':
            i[j] = 2
        i = i.astype(float)
        hess_n[(i[0]-1)*3+i[1],(i[2]-1)*3+i[3]] = i[5]
        hess_n[(i[2]-1)*3+i[3],(i[0]-1)*3+i[1]] = i[5]
    return hess_n

def Hes_mass(mat,fil,nmol):
    m = prop.mass(fil,nmol)
    m2 = np.zeros(len(m)*3)
    for i in np.arange(len(m)):
      m2[i*3:(i+1)*3] = 1/np.sqrt(m[i])
    m3 = m2.copy()
    for i in np.arange(nmol-1):
      m3 = np.append(m3,m2)
    m3 = np.diag(m3)
    #Mmat = m4*mat*m4
    Mmat = np.dot(m3,np.dot(mat,m3))
#    Mmat = np.multiply(m3,np.multiply(m3,mat.T).T)
    return Mmat

def num_wvn(fil,key,nmol,step):
    ###Takes in a target file and outputs the numerical wavenumbers
    ###Calculates the numerical hessian
    conv = 418.4 #Converting kcal to g*Ang**2/ps**2
    lsp = 0.0299792458 #Speed of light in cm/ps
    #
    os.system('testhess %s -k %s Y Y G %s Ang > TH_G.txt'%(fil,key,step))
    hes_anl = Hes_analyt('%s'%(fil),'%s'%(key))
    hes_num = Hes_numer('TH_G.txt',hes_anl.copy())
    hes_num = Hes_mass(hes_num,fil,nmol)
    eig = np.linalg.eigvalsh(hes_num)
    Nwvn = (np.sqrt(conv)/(2*np.pi*lsp))*np.multiply(np.sign(eig),np.sqrt(np.absolute(eig)))
    os.system('rm TH_G.txt')
    return np.sort(Nwvn)



