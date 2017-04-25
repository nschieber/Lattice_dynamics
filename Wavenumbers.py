#!/usr/bin/env python

import os
import sys
import subprocess
import numpy as np
import Expand as Ex
import Properties as Pr


##########################################
#####             Input              #####
##########################################
def Call_Wavenumbers(Method,*positional_parameters,**keyword_parameters):
    """
    This function helps to direct how the wavenumbers will be calculated and calls other functions to calculate and return the wavenumbers

    **Required Inputs
    Method = Harmonic approximation ('HA');
             Stepwise Isotropic QHA ('SiQ');
             Stepwise Isotropic QHA w/ Gruneisen Parameter ('SiQg');
             Gradient Isotropic QHA ('GiQ');
             Gradient Isotropic QHA w/ Gruneisen Parameter ('GiQg');
             Gradient Anisotropic QHA ('GaQ');
             Gradient Anistoropic QHA w/ Gruneisen Parameter ('GaQg');

    **Optional Inputs
    Output = Name of file to put wavenumbers into, if it already exists it will be loaded
    Gruneisen = Gruneisen parameters found with Setup_Isotropic_Gruneisen
    Wavenumber_Reference = Reference wavenumbers for Gruneisen parameter, will be output from Setup_Isotropic_Gruneisen
    Volume_Reference = Reference volume of strucutre for Wavenumber_Reference, will be output from Setup_Isotropic_Gruneisen
    New_Volume = Volume of new structure to calculate wavenumbers for
    Gruneisen_Vol_FracStep = Volumetric stepsize to expand lattice minimum structure to numerically determine the Gruneisen parameter
    molecules_in_cord = number of molecules in Coordinate_file
    Coordinate_file = File containing lattice parameters and atom coordinates
    Parameter_file = Optional input for program
    Program = 'tink' for Tinker Molecular Modeling
              'test' for a test run
    """

    keyword_parameters['Output',
                       'Gruneisen','Wavenumber_Reference','Volume_Reference','New_Volume',
                       'Gruneisen_Vol_FracStep','molecules_in_coord',
                       'Coordinate_file','Parameter_file',
                       'Program']

    if os.path.isfile(Output) == True:
      #If the output file exists, it will be opened and returned
      Wavenumbers = np.load(Output)
      return Wavenumbers

    if ('Program' in keyword_parameters):
      Program = keyword_parameters['Program']

    else:
      if Method == ('SiQg' or 'GiQg'):
        # Methods that use the gruneisen parameter
        if ('Gruneisen' in keyword_parameters) and ('Wavenumber_Reference' in keyword_parameters) and ('Volume_Reference' in keyword_parameters) and ('New_Volume' in keyword_parameters):
          # Assigning reference parameters for the isotropic Gruneisen parameter
          Gruneisen = keyword_parameters['Gruneisen']
          Wavenumber_Reference = keyword_parameters['Wavenumber_Reference']
          Volume_Reference = keyword_parameters['Volume_Reference']
          New_Volume = keyword_parameters['New_Volume']

          # Calculating the wavenumbers of the new Isotropically expanded structure 
          Wavenumbers = Get_Iso_Gruneisen_Wavenumbers(Gruneisen,Wavenumber_Reference,Volume_Reference,New_Volume)
          return Wavenumbers

        else:
          # If the Gruneisen parameter has yet to be determined, here it will be calculated
          # It is assumed that the input Coordinate_file is the lattice minimum strucutre
          Gruneisen, Wavenumber_Reference, Volume_Reference = Setup_Isotropic_Gruneisen(Coordinate_file,Parameter_file,Program,Gruneisen_Vol_FracStep,molecules_in_coord)
          return Gruneisen, Wavenumber_Reverence, Volume_Reference

      elif Method == ('GaQg'):
        # Anisotropic Gruneisen parameter will be a little different than the other two
        print "Anisotropic Gruneisen parameter not yet implimented here"
        sys.exit()

      elif Method == ('SiQ' or 'GiQ' or 'GaQ'):
        # Directly computing the wavenumbers for a specific program, given a coordinate file
        Coordinate_file = keyword_parameters['Coordinate_file'] # Assigning the coordinate file
        
        if Program == 'tink':
          Parameter_file = keyword_parameters['Parameter_file'] # Assigning the parameter file
          Wavenumbers = Tinker_Wavenumber(Coordinate_file,Parameter_file )
        elif Program == 'test':
          Wavenumbers = Test_Wavenumber(Coordinate_file)
        return Wavenumbers


##########################################
#####   TINKER MOLECULAR MODELING    #####
##########################################
def Tinker_Wavenumber(Coordinate_file,Parameter_file): #Was Tink_WVN
    """
    Calls the vibrate executable of Tinker Molecular Modeling and extracts the wavenumbers

    **Required Inputs
    Coordinate_file = Tinker .xyz file for crystal structure
    Parameter_file = Tinker .key file specifying the force field parameter
    """

    # Calling Tinker's vibrate executable and extracting the eigenvalues and wavenumbers of the respective Hessian and mass-weighted Hessian
    eigenvalues_and_wavenumbers = subprocess.check_output("vibrate %s -k %s  CR |  grep -oP '[-+]*[0-9]*\.[0-9]{2,9}'"%(Coordinate_file,Parameter_file),shell=True)

    # Splitting the outputs into array form
    eigenvalues_and_wavenumbers = eigenvalues_and_wavenumbers.split('\n')
    eigenvalues_and_wavenumbers_hold = []
    for i in eigwvn:
      if i == '':
        pass
      else:
        eigenvalues_and_wavenumbers_hold.append(float(i))

    # Extracting the wavenumbers and assuring they're sorted from lowest to highest
    Wavenumbers = np.sort(np.array(eigwvn2[len(eigwvn2)/2:]))
    return Wavenumbers


##########################################
#####              Test              #####
##########################################
def Test_Wavenumber(Coordinate_file):
    """
    This function takes a set of lattice parameters in a .npy file and returns a set of wavenumbers
    Random funcitons can be input here to run different tests and implimented new methods efficiently

    **Required Inputs
    Coordinate_file = File containing lattice parameters and atom coordinates
    """
    Wavenumbers = 0.
    return Wavenumbers


##########################################
#####           Gruneisen            #####
##########################################
def Setup_Isotropic_Gruneisen(Coordinate_file,Program,Gruneisen_Vol_FracStep,molecules_in_coord,*positional_parameters,**keyword_parameters):
    """
    This function calculates the Isotropic Gruneisen parameters for a given coordinate file.
    Calculated numerically given a specified volume fraction stepsize
    ******Eventually! Impliment a second order Gruneisen parameter in here

    **Required Inputs
    Coordinate_file = File containing lattice parameters and atom coordinates
    Program = 'tink' for Tinker Molecular Modeling
              'test' for a test run
    Gruneisen_Vol_FracStep = Volumetric stepsize to expand lattice minimum structure to numerically determine the Gruneisen parameter
    molecules_in_cord = number of molecules in Coordinate_file

    **Optional inputs
    Parameter_file = Optional input for program
    """

    keyword_parameters['Parameter_file']

    # Change in lattice parameters for expanded structure
    dLattice_Parameters = Ex.dvec_Iso((1+Gruneisen_Vol_FracStep),Program,Coordinate_file)
    # Creating expanded structure and putting it in temp.*
    Ex.Expand(Coordinate_file,Parameter_file,Program,dLattice_Parameters,molecules_in_coord,'temp')

    if Program == 'tink':
      Parameter_file = keyword_parameters['Parameter_file'] # Keyfile for Tinker
      Wavenumber_Reference = Tink_WVN(Coordinate_file,Parameter_file) # Wavenumbers for lattice structure
      Wavenumber_expand = Tink_WVN('temp.xyz',Parameter_file) # Wavenumbers expanded strucutre
      lattice_parameters = Pr.Tink_LATVEC(Coordinate_file) # Extracting lattice parameters from lattice structure
      file_ending = '.xyz' # File ending for coordinate files
    elif Program == 'test':
      Wavenumber_Reference = Test_WVN(Coordinate_file,Parameter_file) # Wavenumbers for lattice structure
      Wavenumber_expand = Test_WVN('temp.npy',Parameter_file) # Wavenumbers expanded strucutre
      lattice_parameters = Pr.Test_LATVEC(Coordinate_file) # Extracting lattice parameters from lattice structure
      file_ending = '.npy' # File ending for 'coordinate' files

    Volume_Reference = Pr.Volume(lattice_parameters) # Reference volume
    Volume_expand = Volume_Reference + Gruneisen_Vol_FracStep*Volume_Reference # Volume of expanded strucutre

    Gruneisen = -(np.log(Wavenumber_Reference) - np.log(Wavenumber_expand))/(np.log(Volume_Reference) - np.log(Volume_expand))
    Gruneisen[:3] = 0
    os.system('rm temp'+file_ending)
    return Gruneisen, Wavenumber_Reference, Volume_Reference


def Get_Iso_Gruneisen_Wavenumbers(Gruneisen,Wavenumber_Reference,Volume_Reference,New_Volume): #Was Iso_GRU_New
    """
    This function calculates new wavenumber for an isotropically expanded strucutre using the gruneisen parameter
    ******Eventually! Impliment a second order Gruneisen parameter in here

    **Required Inputs
    Gruneisen = Gruneisen parameters found with Setup_Isotropic_Gruneisen
    Wavenumber_Reference = Reference wavenumbers for Gruneisen parameter, will be output from Setup_Isotropic_Gruneisen
    Volume_Reference = Reference volume of strucutre for Wavenumber_Reference, will be output from Setup_Isotropic_Gruneisen
    New_Volume = Volume of new structure to calculate wavenumbers for
    """

    Wavenumbers = np.dot(Wavenumber_Reference,np.diag(np.power(New_Volume/Volume_Reference,-1*Gruneisen)))
    return Wavenumbers



