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
def Call_Wavenumbers(Method, min_RMS_gradient, **keyword_parameters):
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
    Crystal_matrix_Reference
    New_Crystal_matrix
    Gruneisen_Lat_FracStep
    """
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
            # If there is a saved Gruneisen parameter and set of wavenumbers
            if os.path.isfile(keyword_parameters['Output'] + '_GRUwvn_' + Method + '.npy') and os.path.isfile(
                                            keyword_parameters['Output'] + '_GRU_' + Method + '.npy'):
                print "   ...Using Gruneisen parameters from: " + keyword_parameters['Output'] + '_GRU_' \
                      + Method + '.npy'
                Gruneisen = np.load(keyword_parameters['Output'] + '_GRU_' + Method + '.npy')
                Wavenumber_Reference = np.load(keyword_parameters['Output'] + '_GRUwvn_' + Method + '.npy')
                Volume_Reference = Pr.Volume(Coordinate_file=keyword_parameters['Coordinate_file'],
                                             Program=keyword_parameters['Program'],
                                             Parameter_file=keyword_parameters['Parameter_file'])
            # If the Gruneisen parameter has yet to be determined, here it will be calculated
            # It is assumed that the input Coordinate_file is the lattice minimum strucutre
            else:
                Gruneisen, Wavenumber_Reference, Volume_Reference = \
                    Setup_Isotropic_Gruneisen(keyword_parameters['Coordinate_file'],
                                              keyword_parameters['Program'],
                                              keyword_parameters['Gruneisen_Vol_FracStep'],
                                              keyword_parameters['molecules_in_coord'], min_RMS_gradient,
                                              Parameter_file=keyword_parameters['Parameter_file'])
                print "   ... Saving reference wavenumbers and Gruneisen parameters to: " + \
                      keyword_parameters['Output'] + '_GRU_' + Method + '.npy'
                np.save(keyword_parameters['Output'] + '_GRU_' + Method, Gruneisen)
                np.save(keyword_parameters['Output'] + '_GRUwvn_' + Method, Wavenumber_Reference)
            return Gruneisen, Wavenumber_Reference, Volume_Reference

    elif Method == 'GaQg':
        if ('Gruneisen' in keyword_parameters) and ('Wavenumber_Reference' in keyword_parameters) and \
                ('Crystal_matrix_Reference' in keyword_parameters) and ('New_Crystal_matrix' in keyword_parameters):
            # Calculating the wavenumbers of the new Isotropically expanded structure
            wavenumbers = Get_Aniso_Gruneisen_Wavenumbers(keyword_parameters['Gruneisen'],
                                                          keyword_parameters['Wavenumber_Reference'],
                                                          keyword_parameters['Crystal_matrix_Reference'],
                                                          keyword_parameters['New_Crystal_matrix'])
            return wavenumbers

        else:
            if os.path.isfile(keyword_parameters['Output'] + '_GRUwvn_' + Method + '.npy') and os.path.isfile(
                                            keyword_parameters['Output'] + '_GRU_' + Method + '.npy'):
                print "   ...Using Gruneisen parameters from: " + keyword_parameters['Output'] + '_GRU_' \
                      + Method + '.npy'
                Gruneisen = np.load(keyword_parameters['Output'] + '_GRU_' + Method + '.npy')
                Wavenumber_Reference = np.load(keyword_parameters['Output'] + '_GRUwvn_' + Method + '.npy')
                if keyword_parameters['Program'] == 'Tinker':
                    Crystal_matrix_Reference = Ex.Lattice_parameters_to_Crystal_matrix(
                        Pr.Tinker_Lattice_Parameters(keyword_parameters['Coordinate_file']))
                elif keyword_parameters['Program'] == 'Test':
                    Crystal_matrix_Reference = Ex.Lattice_parameters_to_Crystal_matrix(
                        Pr.Test_Lattice_Parameters(keyword_parameters['Coordinate_file']))
		elif keyword_parameters['Program'] == 'CP2K':
                    Crystal_matrix_Reference = Ex.Lattice_parameters_to_Crystal_matrix(
                        Pr.CP2K_Lattice_Parameters(keyword_parameters['Coordinate_file']))
            # If the Gruneisen parameter has yet to be determined, here it will be calculated
            # It is assumed that the input Coordinate_file is the lattice minimum strucutre
            else:
                Gruneisen, Wavenumber_Reference, Crystal_matrix_Reference = \
                    Setup_Anisotropic_Gruneisen(keyword_parameters['Coordinate_file'], keyword_parameters['Program'],
                                                keyword_parameters['Gruneisen_Lat_FracStep'],
                                                keyword_parameters['molecules_in_coord'],
                                                min_RMS_gradient,
                                                Parameter_file=keyword_parameters['Parameter_file'])
                print "   ... Saving reference wavenumbers and Gruneisen parameters to: " + \
                      keyword_parameters['Output'] + '_GRU_/_GRUwvn' + Method + '.npy'
                np.save(keyword_parameters['Output'] + '_GRU_' + Method, Gruneisen)
                np.save(keyword_parameters['Output'] + '_GRUwvn_' + Method, Wavenumber_Reference)
            return Gruneisen, Wavenumber_Reference, Crystal_matrix_Reference

    elif (Method == 'SiQ') or (Method == 'GiQ') or (Method == 'GaQ') or (Method == 'HA'):
        # Directly computing the wavenumbers for a specific program, given a coordinate file
        if keyword_parameters['Program'] == 'Tinker':
            wavenumbers = Tinker_Wavenumber(keyword_parameters['Coordinate_file'], keyword_parameters['Parameter_file'])
	elif keyword_parameters['Program'] == 'CP2K':
	    wavenumbers = CP2K_Wavenumber(keyword_parameters['Coordinate_file'], keyword_parameters['Parameter_file'],keyword_parameters['cp2kroot'])
        elif keyword_parameters['Program'] == 'Test':
            wavenumbers = Test_Wavenumber(keyword_parameters['Coordinate_file'])
        return wavenumbers


##########################################
#       TINKER MOLECULAR MODELING        #
##########################################
def Tinker_Wavenumber(Coordinate_file, Parameter_file):
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
#                  CP2K                  #
##########################################

def CP2K_Wavenumber(coordinatefile, parameter_file, cp2kroot):
    wavenumbers = np.zeros((3,))
    wavenumfile = open(cp2kroot+'-VIBRATIONS-1.mol','r')
    lines = wavenumfile.readlines()
    iter = 2
    while '[FR-COORD]' not in lines[iter]:
        wave = lines[iter].split()
        wavenumbers = np.append(wavenumbers, float(wave[0]))
	iter = iter+1

    return wavenumbers

	 
    

##########################################
#                  Test                  #
##########################################
def Test_Wavenumber(Coordinate_file, function='Test2'):
    """
    This function takes a set of lattice parameters in a .npy file and returns a set of wavenumbers
    Random functions can be input here to run different tests and implimented new methods efficiently

    **Required Inputs
    Coordinate_file = File containing lattice parameters and atom coordinates
    """

    if function == 'Test1':
        wavenumbers = np.array([0., 0., 0., 52., 380., 1570., 3002.])
        lattice_parameters = np.load(Coordinate_file)
        for i in np.arange(3, len(wavenumbers[3:])+3):  # probably should be 3?
            wavenumbers[i] = wavenumbers[i]*(1/3.)*(((lattice_parameters[0]-16)/6)**2 + ((lattice_parameters[1] -
                                                                                          12)/5)**2 +
                                                    ((lattice_parameters[2] - 23)/11)**2)
    elif function == 'Test2':
        wavenumbers = np.arange(1, 200)
        wavenumbers = wavenumbers**(5.0/3.0)  # get something in the right range = 400^(4/3) = 2941
        wavenumbers[0:3] = [0, 0, 0]  # zero translation
        lattice_parameters = np.load(Coordinate_file)
        [refx, refy, refz] = [lattice_parameters[0]/10, lattice_parameters[1]/7, lattice_parameters[2]/12]
#        [refx,refy,refz] = [lattice_parameters[0]/5,lattice_parameters[1]/6,lattice_parameters[2]/8]
        for i in range(3,len(wavenumbers[3:])+3):
            wavenumbers[i] = wavenumbers[i]*(1.0/15.0)*(2*refx**4.8 + 10*refy**4.2 + 4*np.sin(2*np.pi*refx) +
                                                        3*refz**4.8)
    return wavenumbers


##########################################
#               Gruneisen                #
##########################################
def Setup_Isotropic_Gruneisen(Coordinate_file, Program, Gruneisen_Vol_FracStep, molecules_in_coord, min_RMS_gradient,
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
        Ex.Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'temp', min_RMS_gradient,
                            dlattice_parameters=dLattice_Parameters,
                            Parameter_file=keyword_parameters['Parameter_file'])
        Organized_wavenumbers = Tinker_Gru_organized_wavenumbers('Isotropic', Coordinate_file, 'temp.xyz', keyword_parameters['Parameter_file'])
        Wavenumber_Reference = Organized_wavenumbers[0] 
        Wavenumber_expand = Organized_wavenumbers[1]
        lattice_parameters = Pr.Tinker_Lattice_Parameters(Coordinate_file)
        file_ending = '.xyz'
    if Program == 'CP2K':
        Ex.Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'temp', min_RMS_gradient,
                            dlattice_parameters=dLattice_Parameters,
                            Parameter_file=keyword_parameters['Parameter_file'], cp2kroot = keyword_parameters['cp2kroot'])
        Organized_wavenumbers = CP2K_Gru_organized_wavenumbers('Isotropic', Coordinate_file, 'temp.xyz', keyword_parameters['Parameter_file'])
        Wavenumber_Reference = Organized_wavenumbers[0] 
        Wavenumber_expand = Organized_wavenumbers[1]
        lattice_parameters = Pr.CP2K_Lattice_Parameters(Coordinate_file)
        file_ending = '.pdb'
    elif Program == 'Test':
        Ex.Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'temp', min_RMS_gradient,
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


def Get_Iso_Gruneisen_Wavenumbers(Gruneisen, Wavenumber_Reference, Volume_Reference, New_Volume): 
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


def Setup_Anisotropic_Gruneisen(Coordinate_file, Program, Lattice_FracStep, molecules_in_coord, min_RMS_gradient, **keyword_parameters):
    # Determining the changes in the crystal matrix
    dcrystal_matrix_hold = Ex.Change_Crystal_Matrix(Lattice_FracStep, Program, Coordinate_file)

    # Setting the order in which the crystal matrix will be pulled from
    matrix_order = np.matrix([[0, 0], [1, 1], [2, 2], [0, 1], [0, 2], [1, 2]])

    for i in range(6):
        # Creating 6 expanded strucutures, each in one of the 6 lattice directions
        dcrystal_matrix = np.zeros((3, 3))
        dcrystal_matrix[matrix_order[i, 0], matrix_order[i, 1]] = dcrystal_matrix_hold[matrix_order[i, 0],
                                                                                       matrix_order[i, 1]]

        Ex.Expand_Structure(Coordinate_file, Program, 'crystal_matrix', molecules_in_coord, 'temp_' + str(i), min_RMS_gradient,
                            dcrystal_matrix=dcrystal_matrix, Parameter_file=keyword_parameters['Parameter_file'], cp2kroot = keyword_parameters['cp2kroot'])
        
    if Program == 'Tinker':
        Wavenumber_Reference = Tinker_Wavenumber(Coordinate_file, keyword_parameters['Parameter_file'], keyword_parameters['cp2kroot'])
        expanded_coordinates = ['temp_0.xyz','temp_1.xyz','temp_2.xyz','temp_3.xyz','temp_4.xyz','temp_5.xyz']
        lattice_parameters = Pr.Tinker_Lattice_Parameters(Coordinate_file)
        Crystal_matrix_Reference = Ex.Lattice_parameters_to_Crystal_matrix(lattice_parameters) 
        Organized_wavenumbers = Tinker_Gru_organized_wavenumbers('Anisotropic', Coordinate_file, expanded_coordinates, keyword_parameters['Parameter_file'])
        Wavenumber_Reference = Organized_wavenumbers[0]
        Gruneisen = np.zeros((len(Wavenumber_Reference), 3, 3))
        Wavenumber_expand = Organized_wavenumbers[1:]
        for i in range(6):
            lattice_parameters_expand = Pr.Tinker_Lattice_Parameters(expanded_coordinates[i])
            Crystal_matrix_Reference_expand = Ex.Lattice_parameters_to_Crystal_matrix(lattice_parameters_expand)
            eta = (Crystal_matrix_Reference_expand - Crystal_matrix_Reference)[matrix_order[i, 0], matrix_order[i, 1]] \
                / Crystal_matrix_Reference[matrix_order[i, 0], matrix_order[i, 1]]
            Gruneisen[3:, matrix_order[i, 0], matrix_order[i, 1]] = -(np.log(Wavenumber_expand[i, 3:]) -
                                                                      np.log(Wavenumber_Reference[3:])) / eta
            os.system('rm ' + expanded_coordinates[i])

    elif Program == 'CP2K':
        Wavenumber_Reference = CP2K_Wavenumber(Coordinate_file, keyword_parameters['Parameter_file'], keyword_parameters['cp2kroot'])
        expanded_coordinates = ['temp_0.pdb','temp_1.pdb','temp_2.pdb','temp_3.pdb','temp_4.pdb','temp_5.pdb']
        lattice_parameters = Pr.CP2K_Lattice_Parameters(Coordinate_file)
        Crystal_matrix_Reference = Ex.Lattice_parameters_to_Crystal_matrix(lattice_parameters) 
        Organized_wavenumbers = CP2K_Gru_organized_wavenumbers('Anisotropic', Coordinate_file, expanded_coordinates, keyword_parameters['Parameter_file'])
        Wavenumber_Reference = Organized_wavenumbers[0]
        Gruneisen = np.zeros((len(Wavenumber_Reference), 3, 3))
        Wavenumber_expand = Organized_wavenumbers[1:]
        for i in range(6):
            lattice_parameters_expand = Pr.CP2K_Lattice_Parameters(expanded_coordinates[i])
            Crystal_matrix_Reference_expand = Ex.Lattice_parameters_to_Crystal_matrix(lattice_parameters_expand)
            eta = (Crystal_matrix_Reference_expand - Crystal_matrix_Reference)[matrix_order[i, 0], matrix_order[i, 1]] \
                / Crystal_matrix_Reference[matrix_order[i, 0], matrix_order[i, 1]]
            Gruneisen[3:, matrix_order[i, 0], matrix_order[i, 1]] = -(np.log(Wavenumber_expand[i, 3:]) -
                                                                      np.log(Wavenumber_Reference[3:])) / eta
            os.system('rm ' + expanded_coordinates[i])

    elif Program == 'Test':
        expanded_coordinates = ['temp_0.npy','temp_1.npy','temp_2.npy','temp_3.npy','temp_4.npy','temp_5.npy']
        Wavenumber_Reference = Test_Wavenumber(Coordinate_file)
        Gruneisen = np.zeros((len(Wavenumber_Reference), 3, 3))
        lattice_parameters = Pr.Test_Lattice_Parameters(Coordinate_file)
        keyword_parameters['Parameter_file'] = ''
        Crystal_matrix_Reference = Ex.Lattice_parameters_to_Crystal_matrix(lattice_parameters)
        for i in range(6):
            Wavenumber_expand = Test_Wavenumber(expanded_coordinates[i])
            lattice_parameters_expand = Pr.Test_Lattice_Parameters(expanded_coordinates[i])
            Crystal_matrix_Reference_expand = Ex.Lattice_parameters_to_Crystal_matrix(lattice_parameters_expand)
            eta = (Crystal_matrix_Reference_expand - Crystal_matrix_Reference)[matrix_order[i, 0], matrix_order[i, 1]] \
                / Crystal_matrix_Reference[matrix_order[i, 0], matrix_order[i, 1]]
            Gruneisen[3:, matrix_order[i, 0], matrix_order[i, 1]] = -(np.log(Wavenumber_expand[3:]) -
                                                                      np.log(Wavenumber_Reference[3:])) / eta
            os.system('rm ' + expanded_coordinates[i])
    return Gruneisen, Wavenumber_Reference, Crystal_matrix_Reference


def Get_Aniso_Gruneisen_Wavenumbers(Gruneisen, Wavenumber_Reference, Crystal_matrix_Reference, New_Crystal_matrix):
    matrix_order = np.matrix([[0, 0], [1, 1], [2, 2], [0, 1], [0, 2], [1, 2]])

    wavenumbers = np.zeros(len(Wavenumber_Reference))

    for i in np.arange(3, len(wavenumbers), 1):
        hold = 0
        for j in range(6):
            hold = hold + -1 * (New_Crystal_matrix[matrix_order[j, 0], matrix_order[j, 1]] -
                                Crystal_matrix_Reference[matrix_order[j, 0], matrix_order[j, 1]]) / \
                          Crystal_matrix_Reference[matrix_order[j, 0], matrix_order[j, 1]] * \
                          Gruneisen[i, matrix_order[j, 0], matrix_order[j, 1]]
        wavenumbers[i] = Wavenumber_Reference[i] * np.exp(hold)
    return wavenumbers

def Tinker_Gru_organized_wavenumbers(Expansion_type, Coordinate_file, Expanded_Coordinate_file, Parameter_file):
    from munkres import Munkres, print_matrix
    m = Munkres()

    number_of_modes = 3*Pr.Tinker_atoms_per_molecule(Coordinate_file, 1)

    if Expansion_type == 'Isotropic':
        wavenumbers = np.zeros((2, number_of_modes))
        eigenvectors = np.zeros((2, number_of_modes, number_of_modes))
        wavenumbers[0], eigenvectors[0] = Tinker_Wavenumber_and_Vectors(Coordinate_file, Parameter_file)
        wavenumbers[1], eigenvectors[1] = Tinker_Wavenumber_and_Vectors(Expanded_Coordinate_file, Parameter_file)
    elif Expansion_type == 'Anisotropic':
        wavenumbers = np.zeros((7, number_of_modes))
        eigenvectors = np.zeros((7, number_of_modes, number_of_modes))
        wavenumbers[0], eigenvectors[0] = Tinker_Wavenumber_and_Vectors(Coordinate_file, Parameter_file)
        for i in xrange(1,7):
            wavenumbers[i], eigenvectors[i] = Tinker_Wavenumber_and_Vectors(Expanded_Coordinate_file[i-1], Parameter_file)


    # Weighting the modes matched together
    wavenumbers_out = np.zeros((len(wavenumbers[:, 0]), number_of_modes))
    wavenumbers_out[0] = wavenumbers[0]
    for k in xrange(1, len(wavenumbers[:, 0])):
        weight = np.zeros((number_of_modes - 3, number_of_modes - 3))
        for i in xrange(3, number_of_modes):
            diff = np.linalg.norm(np.dot(eigenvectors[0, i], eigenvectors[k, i]))/(np.linalg.norm(eigenvectors[0, i])*np.linalg.norm(eigenvectors[k, i]))
            if diff > 0.95:
                weight[i - 3] = 10000000.
                weight[i - 3, i - 3] = 1. - diff
            else:
                for j in xrange(3, number_of_modes):
                    weight[i - 3, j - 3] = 1 - np.linalg.norm(np.dot(eigenvectors[0, i], eigenvectors[k, j]))/(np.linalg.norm(eigenvectors[0, i])*np.linalg.norm(eigenvectors[k, j]))

        # Using the Hungarian algorithm to match wavenumbers
        Wgt = m.compute(weight)
        x,y = zip(*Wgt)
        z = np.column_stack((x,y))
        z = z +3

    # Re-organizing the expanded wavenumbers
        for i in z:
            wavenumbers_out[k, i[0]] = wavenumbers[k, i[1]]
    return wavenumbers_out

def CP2K_Gru_organized_wavenumbers(Expansion_type, Coordinate_file, Expanded_Coordinate_file, Parameter_file):
    from munkres import Munkres, print_matrix
    m = Munkres()

    number_of_modes = 3*Pr.CP2K_atoms_per_molecule(Coordinate_file, 1)

    if Expansion_type == 'Isotropic':
        wavenumbers = np.zeros((2, number_of_modes))
        eigenvectors = np.zeros((2, number_of_modes, number_of_modes))
        wavenumbers[0], eigenvectors[0] = CP2K_Wavenumber_and_Vectors(Coordinate_file, Parameter_file)
        wavenumbers[1], eigenvectors[1] = CP2K_Wavenumber_and_Vectors(Expanded_Coordinate_file, Parameter_file)
    elif Expansion_type == 'Anisotropic':
        wavenumbers = np.zeros((7, number_of_modes))
        eigenvectors = np.zeros((7, number_of_modes, number_of_modes))
        wavenumbers[0], eigenvectors[0] = CP2K_Wavenumber_and_Vectors(Coordinate_file, Parameter_file)
        for i in xrange(1,7):
            wavenumbers[i], eigenvectors[i] = CP2K_Wavenumber_and_Vectors(Expanded_Coordinate_file[i-1], Parameter_file)


    # Weighting the modes matched together
    wavenumbers_out = np.zeros((len(wavenumbers[:, 0]), number_of_modes))
    wavenumbers_out[0] = wavenumbers[0]
    for k in xrange(1, len(wavenumbers[:, 0])):
        weight = np.zeros((number_of_modes - 3, number_of_modes - 3))
        for i in xrange(3, number_of_modes):
            diff = np.linalg.norm(np.dot(eigenvectors[0, i], eigenvectors[k, i]))/(np.linalg.norm(eigenvectors[0, i])*np.linalg.norm(eigenvectors[k, i]))
            if diff > 0.95:
                weight[i - 3] = 10000000.
                weight[i - 3, i - 3] = 1. - diff
            else:
                for j in xrange(3, number_of_modes):
                    weight[i - 3, j - 3] = 1 - np.linalg.norm(np.dot(eigenvectors[0, i], eigenvectors[k, j]))/(np.linalg.norm(eigenvectors[0, i])*np.linalg.norm(eigenvectors[k, j]))

        # Using the Hungarian algorithm to match wavenumbers
        Wgt = m.compute(weight)
        x,y = zip(*Wgt)
        z = np.column_stack((x,y))
        z = z +3

    # Re-organizing the expanded wavenumbers
        for i in z:
            wavenumbers_out[k, i[0]] = wavenumbers[k, i[1]]
    return wavenumbers_out


def Tinker_Wavenumber_and_Vectors(Coordinate_file, Parameter_file):
    # Calling Tinker's vibrate executable and extracting the eigenvectors and wavenumbers of the respective
    # Hessian and mass-weighted Hessian
    os.system('cp ' + Coordinate_file + ' vector_temp.xyz')
    output = subprocess.check_output("vibrate vector_temp.xyz -k %s  A |  grep -oP '[-+]*[0-9]*\.[0-9]{2,9}'"
                                                          % (Parameter_file), shell=True)

    os.system('rm vector_temp.*')

    # Finding the number modes in the system
    number_of_modes = 3*Pr.Tinker_atoms_per_molecule(Coordinate_file, 1)

    # Splitting the outputs into array form
    output = output.split('\n')
    output.remove('')
    output = np.array(output).astype(float)

    # Grabbing the wavenumbers
    wavenumbers = np.array(output[number_of_modes: number_of_modes*2]).astype(float)

    # Grabbing the eigenvectors
    eigenvectors = np.zeros((number_of_modes, number_of_modes))
    for i in range(number_of_modes):
        start = number_of_modes*(i + 2) + i + 1
        finish = start + number_of_modes
        eigenvectors[i] = output[start: finish] /np.sqrt(np.sum(output[start: finish]**2))
    return wavenumbers, eigenvectors


def CP2K_Wavenumber_and_Vectors(Coordinate_file, Parameter_file):
    # Calling Tinker's vibrate executable and extracting the eigenvectors and wavenumbers of the respective
    # Hessian and mass-weighted Hessian
    os.system('cp ' + Coordinate_file + ' vector_temp.xyz')
    output = subprocess.check_output("vibrate vector_temp.xyz -k %s  A |  grep -oP '[-+]*[0-9]*\.[0-9]{2,9}'"
                                                          % (Parameter_file), shell=True)

    os.system('rm vector_temp.*')

    # Finding the number modes in the system
    number_of_modes = 3*Pr.Tinker_atoms_per_molecule(Coordinate_file, 1)

    # Splitting the outputs into array form
    output = output.split('\n')
    output.remove('')
    output = np.array(output).astype(float)

    # Grabbing the wavenumbers
    wavenumbers = np.array(output[number_of_modes: number_of_modes*2]).astype(float)

    # Grabbing the eigenvectors
    eigenvectors = np.zeros((number_of_modes, number_of_modes))
    for i in range(number_of_modes):
        start = number_of_modes*(i + 2) + i + 1
        finish = start + number_of_modes
        eigenvectors[i] = output[start: finish] /np.sqrt(np.sum(output[start: finish]**2))
    return wavenumbers, eigenvectors


