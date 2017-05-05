#!/usr/bin/env python

import os
import sys
import subprocess
import numpy as np
import Expand as Ex
import ThermodynamicProperties as Pr


##########################################
#                 Input                  #
##########################################
def Call_Wavenumbers(Method, **keyword_parameters):
    """
    This function helps to direct how the wavenumbers will be calculated and calls other functions to calculate and 
    return the wavenumbers

    **Required Inputs
    Method = Harmonic approximation ('HA');
             Stepwise Isotropic QHA ('SiQ');
             Stepwise Isotropic QHA w/ Gruneisen Parameter ('SiQg');
             Gradient Isotropic QHA ('GiQ');
             Gradient Isotropic QHA w/ Gruneisen Parameter ('GiQg');
             Gradient Anisotropic QHA ('GaQ');
             Gradient Anisotropic QHA w/ Gruneisen Parameter ('GaQg');

    **Optional Inputs
    Output = Name of file to put wavenumbers into, if it already exists it will be loaded
    Gruneisen = Gruneisen parameters found with Setup_Isotropic_Gruneisen
    Wavenumber_Reference = Reference wavenumbers for Gruneisen parameter, will be output from Setup_Isotropic_Gruneisen
    Volume_Reference = Reference volume of structure for Wavenumber_Reference, will be output from 
    Setup_Isotropic_Gruneisen
    New_Volume = Volume of new structure to calculate wavenumbers for
    Gruneisen_Vol_FracStep = Volumetric stepsize to expand lattice minimum structure to numerically determine the
    Gruneisen parameter
    molecules_in_coord = number of molecules in Coordinate_file
    Coordinate_file = File containing lattice parameters and atom coordinates
    Parameter_file = Optional input for program
    Program = 'Tinker' for Tinker Molecular Modeling
              'Test' for a test run
    """
    # If the output file exists, it will be opened and returned

    if 'Output' is keyword_parameters:
        if os.path.isfile(keyword_parameters['Output']):
            wavenumbers = np.load(keyword_parameters['Output'])
            return wavenumbers

    if (Method == 'SiQg') or (Method == 'GiQg'):
        # Methods that use the Gruneisen parameter
        if ('Gruneisen' in keyword_parameters) and ('Wavenumber_Reference' in keyword_parameters) and \
                ('Volume_Reference' in keyword_parameters) and ('New_Volume' in keyword_parameters):
            # Calculating the wavenumbers of the new Isotropically expanded structure
            wavenumbers = Get_Iso_Gruneisen_Wavenumbers(keyword_parameters['Gruneisen'],
                                                                keyword_parameters['Wavenumber_Reference'],
                                                                keyword_parameters['Volume_Reference'],
                                                                keyword_parameters['New_Volume'])
            return wavenumbers

        else:
            # If the Gruneisen parameter has yet to be determined, here it will be calculated
            # It is assumed that the input Coordinate_file is the lattice minimum strucutre
            Gruneisen, Wavenumber_Reference, Volume_Reference = \
                        Setup_Isotropic_Gruneisen(keyword_parameters['Coordinate_file'],
                                                  keyword_parameters['Program'],
                                                  keyword_parameters['Gruneisen_Vol_FracStep'],
                                                  keyword_parameters['molecules_in_coord'], Parameter_file=
                                                  keyword_parameters['Parameter_file'])
            return Gruneisen, Wavenumber_Reference, Volume_Reference

    elif Method == 'GaQg':
        print "Anisotropic Gruneisen parameter not yet implemented here"
        sys.exit()

    elif (Method == 'SiQ') or (Method == 'GiQ') or (Method == 'GaQ') or (Method == 'HA'):
        # Directly computing the wavenumbers for a specific program, given a coordinate file
        if keyword_parameters['Program'] == 'Tinker':
            wavenumbers = Tinker_Wavenumber(keyword_parameters['Coordinate_file'], keyword_parameters['Parameter_file'])
        elif keyword_parameters['Program'] == 'Test':
            wavenumbers = Test_Wavenumber(keyword_parameters['Coordinate_file'])
        return wavenumbers


##########################################
#       TINKER MOLECULAR MODELING        #
##########################################
def Tinker_Wavenumber(Coordinate_file, Parameter_file): #Was Tink_WVN
    """
    Calls the vibrate executable of Tinker Molecular Modeling and extracts the wavenumbers

    **Required Inputs
    Coordinate_file = Tinker .xyz file for crystal structure
    Parameter_file = Tinker .key file specifying the force field parameter
    """
    # Calling Tinker's vibrate executable and extracting the eigenvalues and wavenumbers of the respective
    # Hessian and mass-weighted Hessian
    eigenvalues_and_wavenumbers = subprocess.check_output("vibrate %s -k %s  CR |  grep -oP '[-+]*[0-9]*\.[0-9]{2,9}'"
                                                          % (Coordinate_file, Parameter_file), shell=True)
    # Splitting the outputs into array form
    eigenvalues_and_wavenumbers = eigenvalues_and_wavenumbers.split('\n')
    eigenvalues_and_wavenumbers_hold = []
    for i in eigenvalues_and_wavenumbers:
        if i == '':
            pass
        else:
            eigenvalues_and_wavenumbers_hold.append(float(i))

    # Extracting the wavenumbers and assuring they're sorted from lowest to highest
    wavenumbers = np.sort(np.array(eigenvalues_and_wavenumbers_hold[len(eigenvalues_and_wavenumbers_hold)/2:]))
    return wavenumbers


##########################################
#                  Test                  #
##########################################
def Test_Wavenumber(Coordinate_file):
    """
    This function takes a set of lattice parameters in a .npy file and returns a set of wavenumbers
    Random funcitons can be input here to run different tests and implimented new methods efficiently

    **Required Inputs
    Coordinate_file = File containing lattice parameters and atom coordinates
    """
    wavenumbers = np.array([0.,0.,0.,52.,380.,1570.,3002.])
    lattice_parameters= np.load(Coordinate_file)
    for i in np.arange(3,len(wavenumbers[3:])+1):
        wavenumbers[i] = wavenumbers[i]*(1/3.)*(((lattice_parameters[0]-16)/6)**2+((lattice_parameters[1]-12)/5)**2+((lattice_parameters[2]-23)/11)**2)
    return wavenumbers


##########################################
#               Gruneisen                #
##########################################
def Setup_Isotropic_Gruneisen(Coordinate_file, Program, Gruneisen_Vol_FracStep, molecules_in_coord,
                              **keyword_parameters):
    """
    This function calculates the Isotropic Gruneisen parameters for a given coordinate file.
    Calculated numerically given a specified volume fraction stepsize
    ******Eventually! Impliment a second order Gruneisen parameter in here

    **Required Inputs
    Coordinate_file = File containing lattice parameters and atom coordinates
    Program = 'Tinker' for Tinker Molecular Modeling
              'Test' for a test run
    Gruneisen_Vol_FracStep = Volumetric stepsize to expand lattice minimum structure to numerically determine the 
    Gruneisen parameter
    molecules_in_coord = number of molecules in Coordinate_file

    **Optional inputs
    Parameter_file = Optional input for program
    """
    # Change in lattice parameters for expanded structure
    dLattice_Parameters = Ex.Isotropic_Change_Lattice_Parameters((1+Gruneisen_Vol_FracStep), Program, Coordinate_file)

    # Determining wavenumbers of lattice strucutre and expanded strucutre
    # Also, assigning a file ending name for the nex coordinate file (program dependent)
    if Program == 'Tinker':
        Ex.Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'temp',
                            dlattice_parameters=dLattice_Parameters, Parameter_file=keyword_parameters['Parameter_file'])
        Wavenumber_Reference = Tinker_Wavenumber(Coordinate_file, keyword_parameters['Parameter_file'])
        Wavenumber_expand = Tinker_Wavenumber('temp.xyz', keyword_parameters['Parameter_file'])
        lattice_parameters = Pr.Tinker_Lattice_Parameters(Coordinate_file)
        file_ending = '.xyz'
    elif Program == 'Test':
        Ex.Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'temp',
                            dlattice_parameters=dLattice_Parameters)
        Wavenumber_Reference = Test_Wavenumber(Coordinate_file)
        Wavenumber_expand = Test_Wavenumber('temp.npy')
        lattice_parameters = Pr.Test_Lattice_Parameters(Coordinate_file)
        file_ending = '.npy'

    # Calculating the volume of the lattice minimum and expanded structure
    Volume_Reference = Pr.Volume(lattice_parameters=lattice_parameters)
    Volume_expand = Volume_Reference + Gruneisen_Vol_FracStep*Volume_Reference

    # Calculating the Gruneisen parameter and zeroing out the parameters for the translational modes
    Gruneisen = np.zeros(len(Wavenumber_Reference))
    Gruneisen[3:] = -(np.log(Wavenumber_Reference[3:]) - np.log(Wavenumber_expand[3:]))/(np.log(Volume_Reference) -
                                                                                 np.log(Volume_expand))

    # Removing extra files created in process
    os.system('rm temp'+file_ending)
    return Gruneisen, Wavenumber_Reference, Volume_Reference


def Get_Iso_Gruneisen_Wavenumbers(Gruneisen, Wavenumber_Reference, Volume_Reference, New_Volume): #Was Iso_GRU_New
    """
    This function calculates new wavenumber for an isotropically expanded strucutre using the gruneisen parameter
    ******Eventually! Impliment a second order Gruneisen parameter in here

    **Required Inputs
    Gruneisen = Gruneisen parameters found with Setup_Isotropic_Gruneisen
    Wavenumber_Reference = Reference wavenumbers for Gruneisen parameter, will be output from Setup_Isotropic_Gruneisen
    Volume_Reference = Reference volume of strucutre for Wavenumber_Reference, will be output from Setup_Isotropic_Gruneisen
    New_Volume = Volume of new structure to calculate wavenumbers for
    """
    wavenumbers = np.diag(np.power(New_Volume/Volume_Reference, -1*Gruneisen))
    wavenumbers = np.dot(Wavenumber_Reference, wavenumbers)
    return wavenumbers



