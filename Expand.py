#!/usr/bin/env python
import os
import sys
import itertools as it
import numpy as np
import ThermodynamicProperties as Pr
import Wavenumbers as Wvn
import subprocess
import fileinput

##########################################
#                 Input                  #
##########################################
def Call_Expansion(Method, Purpose, Program, Coordinate_file, molecules_in_coord, min_RMS_gradient, **keyword_parameters):
    """
    :param Method: Harmonic approximation ('HA');
                   Stepwise Isotropic QHA ('SiQ');
                   Stepwise Isotropic QHA w/ Gruneisen Parameter ('SiQg');
                   Gradient Isotropic QHA ('GiQ');
                   Gradient Isotropic QHA w/ Gruneisen Parameter ('GiQg');
                   Gradient Anisotropic QHA ('GaQ');
                   Gradient Anisotropic QHA w/ Gruneisen Parameter ('GaQg');
    :param Purpose: 'expand' expanding a strucutre
                    'local_gradient' finding the local gradient of expansion for a strucutre
    :param Program: 'Tinker' Tinker molecular modeling
                    'Test' Test case
    :param Coordinate_file: file containing lattice parameters (and coordinates)
    :param molecules_in_coord: number of molecules in Coordinate_file
    :param keyword_parameters: Parameter_file, volume_fraction_chage, matrix_parameters_fraction_change, Temperature,
    Pressure, matrix_parameters_fraction_change, Statistical_mechanics, Gruneisen, 
    Wavenumber_Reference, Volume_Reference, Output
    
    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    volume_fraction_change: fractional volumetric step size (either for numerical gradient or amount to expand by).
    matrix_parameters_fraction_change: fraction to change crystal matrix parameter by
    Temperature: in Kelvin
    Pressure: in atm
    Statistical_mechanics: 'Classical' Classical mechanics
                           'Quantum' Quantum mechanics
    Gruneisen: isotropic Gruneisen parameter
    Wavenumber_reference: reference wavenumbers for the Gruneisen parameter
    Volume_reference: reference volume for the Gruneisen parameter
    Output: string to name expanded coordinate file
    dcrystal_matrix: changes in the crystal matrix
    crystal_matrix_reference
    """
    # If there is no parameter_file, one is just inputted so there is no errors
    if 'Parmaeter_file' in keyword_parameters:
        keyword_parameters['Parameter_file'] == ''

    # Expanded a strucutre being input
    if Purpose == 'expand':
        if (Method == 'GiQ') or (Method == 'GiQg') or (Method == 'SiQ') or (Method == 'SiQg'):
            dlattice_parameters = Isotropic_Change_Lattice_Parameters(keyword_parameters['volume_fraction_change'],
                                                                      Program, Coordinate_file)
            Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord,
                             keyword_parameters['Output'], min_RMS_gradient, 
                             Parameter_file=keyword_parameters['Parameter_file'],
                             dlattice_parameters=dlattice_parameters)
        elif (Method == 'GaQ') or (Method == 'GaQg'):
            Expand_Structure(Coordinate_file, Program, 'crystal_matrix', molecules_in_coord,
                             keyword_parameters['Output'], min_RMS_gradient,
                             Parameter_file=keyword_parameters['Parameter_file'],
                             dcrystal_matrix=keyword_parameters['dcrystal_matrix'])

    # Fining the local gradient of expansion for inputted strucutre
    elif Purpose == 'local_gradient':
        if Method == 'GiQ':
            isotropic_local_gradient, wavenumbers, volume = \
                Isotropic_Local_Gradient(Coordinate_file, Program, keyword_parameters['Temperature'],
                                         keyword_parameters['Pressure'], keyword_parameters['volume_fraction_change'],
                                         molecules_in_coord,
                                         keyword_parameters['Statistical_mechanics'], Method, min_RMS_gradient,
                                         Parameter_file=keyword_parameters['Parameter_file'])
            return isotropic_local_gradient, wavenumbers, volume
        elif Method == 'GiQg':
            isotropic_local_gradient, wavenumbers, volume = \
                Isotropic_Local_Gradient(Coordinate_file, Program, keyword_parameters['Temperature'],
                                         keyword_parameters['Pressure'], keyword_parameters['volume_fraction_change'],
                                         molecules_in_coord,
                                         keyword_parameters['Statistical_mechanics'], Method, min_RMS_gradient,
                                         Parameter_file=keyword_parameters['Parameter_file'],
                                         Gruneisen=keyword_parameters['Gruneisen'],
                                         Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                         Volume_Reference=keyword_parameters['Volume_Reference'])
            return isotropic_local_gradient, wavenumbers, volume
        elif Method == 'GaQ':
            anisotropic_local_gradient, wavenumbers = \
                Anisotropic_Local_Gradient(Coordinate_file, Program, keyword_parameters['Temperature'],
                                           keyword_parameters['Pressure'],
                                           keyword_parameters['matrix_parameters_fraction_change'],
                                           molecules_in_coord,
                                           keyword_parameters['Statistical_mechanics'], Method,
                                           keyword_parameters['Aniso_LocGrad_Type'], min_RMS_gradient,
                                           Parameter_file=keyword_parameters['Parameter_file'])
            return anisotropic_local_gradient, wavenumbers
        elif Method == 'GaQg':
            anisotropic_local_gradient, wavenumbers = \
                Anisotropic_Local_Gradient(Coordinate_file, Program, keyword_parameters['Temperature'],
                                           keyword_parameters['Pressure'],
                                           keyword_parameters['matrix_parameters_fraction_change'],
                                           molecules_in_coord,
                                           keyword_parameters['Statistical_mechanics'], Method,
                                           keyword_parameters['Aniso_LocGrad_Type'], min_RMS_gradient,
                                           Parameter_file=keyword_parameters['Parameter_file'],
                                           Gruneisen=keyword_parameters['Gruneisen'],
                                           Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                           Crystal_matrix_Reference=keyword_parameters['Crystal_matrix_Reference'])
            return anisotropic_local_gradient, wavenumbers


##########################################
#       TINKER MOLECULAR MODELING        #
##########################################
def Return_Tinker_Coordinates(Coordinate_file):
    """
    This function opens a Tinker .xyz for a crystal and returns the 3x(number of atoms) matrix

    **Required Inputs
    Coordinate_file = Tinker .xyz file for a crystal
    """
    with open(Coordinate_file) as f:
        # Opening xyz coordinate file to expand
        coordinates = np.array(list(it.izip_longest(*[lines.split() for lines in f], fillvalue=' '))).T
    coordinates = coordinates[2:, 2:5].astype(float)
    return coordinates


def Output_Tinker_New_Coordinate_File(Coordinate_file, Parameter_file, coordinates, lattice_parameters, Output, min_RMS_gradient):
    """
    This function takes a new set of coordinates and utilizes a previous coordinate file as a template to produce a new
    Tinker .xyz crystal file
    The structure is also minimized

    **Required Inputs
    Coordinate_file = Tinker .xyz file for a crystal
    Parameter_file = Tinker .key file with force field parameters
    coordinates = New coordinates in a 3x(number of atoms) matrix
    lattice_parameters = lattice parameters as an array ([a,b,c,alpha,beta,gamma])
    Output = file name of new .xyz file
    """   
    Ouput_Tinker_Coordinate_File(Coordinate_file, Parameter_file, coordinates, lattice_parameters, Output)
    Tinker_minimization(Parameter_file, Output + '.xyz', Output, min_RMS_gradient)

def Ouput_Tinker_Coordinate_File(Coordinate_file, Parameter_file, coordinates, lattice_parameters, Output):
    """
    This function takes a new set of coordinates and utilizes a previous coordinate file as a template to produce a new
    Tinker .xyz crystal file

    **Required Inputs
    Coordinate_file = Tinker .xyz file for a crystal
    Parameter_file = Tinker .key file with force field parameters
    coordinates = New coordinates in a 3x(number of atoms) matrix
    lattice_parameters = lattice parameters as an array ([a,b,c,alpha,beta,gamma])
    Output = file name of new .xyz file
    """
    with open(Coordinate_file) as f:
        # Opening xyz coordinate file to expand
        coordinates_template = np.array(list(it.izip_longest(*[lines.split() for lines in f], fillvalue=' '))).T

    coordinates_template[2:, 2:5] = np.around(coordinates, decimals=8).astype('str')
    coordinates_template[1, :6] = lattice_parameters.astype(str)
    string_coordinates = ''
    for i in range(len(coordinates_template[:, 0])):
        for j in range(len(coordinates_template[i, :])):
            string_coordinates = string_coordinates + '    ' + coordinates_template[i, j]
        string_coordinates = string_coordinates + '\n'

    with open(Output + '.xyz', 'w') as file_out:
        file_out.write(string_coordinates)


def Tinker_minimization(Parameter_file, Coordinate_file, Output, min_RMS_gradient):
    with open('minimization.out', 'a') as myfile:
        myfile.write("======================== New Minimization ========================\n")
    run_min = True
    count = 0
    subprocess.call(['cp', Coordinate_file, 'Temp_min_0.xyz'])
    while run_min == True:
        output = subprocess.check_output(['minimize', 'Temp_min_' + str(count) + '.xyz', '-k', Parameter_file, str(min_RMS_gradient)])
        count = count + 1
        with open('minimization.out', 'a') as myfile:
            myfile.write(output)
        output = output.split('\n')
        subprocess.call(['mv', 'Temp_min_' + str(count - 1) + '.xyz_2', 'Temp_min_' + str(count) + '.xyz'])
        if output[-6] == ' LBFGS  --  Normal Termination due to SmallGrad':
            run_min = False
            subprocess.call(['mv', 'Temp_min_' + str(count) + '.xyz', Output + '.xyz'])
        elif count == 10:
            run_min = False
            subprocess.call(['mv', 'Temp_min_1.xyz', Output + '.xyz'])
            print "      Could not minimize strucutre to tolerance after 10 shake cycles"
            print "      Using strucutre from first minimzation"
        else:
            if count == 1:
                print "   ... Structure did not minimze to tolerance, shaking molecule and re-minimizing"
            coordinates = Return_Tinker_Coordinates('Temp_min_' + str(count) + '.xyz')
            coordinates = coordinates + np.random.randint(0, 10, size=(len(coordinates), 3))*min_RMS_gradient/100.
            lattice_parameters = Pr.Tinker_Lattice_Parameters('Temp_min_' + str(count) + '.xyz')
            Ouput_Tinker_Coordinate_File('Temp_min_' + str(count) + '.xyz', Parameter_file, coordinates, lattice_parameters, 'Temp_min_' + str(count))
    for i in xrange(11):
        if os.path.isfile('Temp_min_' + str(i) + '.xyz'):
            subprocess.call(['rm', 'Temp_min_' + str(i) + '.xyz'])


##########################################
#       CP2K MOLECULAR MODELING        #
##########################################
def Return_CP2K_Coordinates(Coordinate_file):
    """
    This function opens a Tinker .xyz for a crystal and returns the 3x(number of atoms) matrix

    **Required Inputs
    Coordinate_file = Tinker .xyz file for a crystal
    """
    with open(Coordinate_file) as f:
        # Opening xyz coordinate file to expand
        coordinates = np.array(list(it.izip_longest(*[lines.split() for lines in f], fillvalue=' '))).T
    coords = np.zeros((len(coordinates)-3,3))
    coords[:,:] = coordinates[2:-1, 3:6].astype(float)
    return coords


def Output_CP2K_New_Coordinate_File(Coordinate_file, Parameter_file, coordinates, lattice_parameters, Output, min_RMS_gradient):
    """
    This function takes a new set of coordinates and utilizes a previous coordinate file as a template to produce a new
    Tinker .xyz crystal file
    The structure is also minimized

    **Required Inputs
    Coordinate_file = Tinker .xyz file for a crystal
    Parameter_file = Tinker .key file with force field parameters
    coordinates = New coordinates in a 3x(number of atoms) matrix
    lattice_parameters = lattice parameters as an array ([a,b,c,alpha,beta,gamma])
    Output = file name of new .xyz file
    """   
    Ouput_Tinker_Coordinate_File(Coordinate_file, Parameter_file, coordinates, lattice_parameters, Output)
    CP2K_minimization(Parameter_file, Output + '.xyz', Output, min_RMS_gradient)

def Ouput_CP2K_Coordinate_File(Coordinate_file, Parameter_file, coordinates, lattice_parameters, Output):
    """
    This function takes a new set of coordinates and utilizes a previous coordinate file as a template to produce a new
    .pdb crystal file

    **Required Inputs
    Coordinate_file = Tinker .xyz file for a crystal
    Parameter_file = Tinker .key file with force field parameters
    coordinates = New coordinates in a 3x(number of atoms) matrix
    lattice_parameters = lattice parameters as an array ([a,b,c,alpha,beta,gamma])
    Output = file name of new .xyz file
    """
    xstr = []
    numatoms = np.shape(coordinates)[0]
    for d in range(3):                  
        xstr.append("%9.3f" % (lattice_param[d]))
    for d in range(3,6):                  
        xstr.append("%7.2f" % (lattice_param[d]))
    with open(Output + '.pdb', 'w') as file_out:
        file_out.write('REMARK'+ '\n')
        file_out.write('CRYST1    '+str(lattice_param[0])+xstr[1]+xstr[2]+xstr[3]+xstr[4]+xstr[5]+'\n')
        for x in range(numatoms):

            if (x+1) % 2 == 1:
                ty = 'C'
            else:
                ty = 'H'
            line =  '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4s}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',x+1, ty,'','','','','',coordinates[x,0],coordinates[x,0],coordinates[x,0], 0.0,0.0,ty,'')
            file_out.write(line+'\n')
        file_out.write('END')


def CP2K_minimization(Parameter_file, Coordinate_file, Output, min_RMS_gradient):
    new_coord_file = 'Temp_min_0.pdb'
    subprocess.call(['cp', Coordinate_file, new_coord_file])
    new_param_file = 'Temp_min_0.inp'
    subprocess.call(['cp', parameter_file, new_param_file])
    #lattice_param = Pr.CP2K_Lattice_Parameters(Coordinate_file)
    volume = '{:4.8s}'.format(str(lattice_param[0]*lattice_param[2]*lattice_param[1]))
    
    tempnamestr = 's/TEMPORARY/BNZ_VOL_'+volume+'/g'
    subprocess.call(["sed", "-i", "-e",  tempnamestr, new_param_file])
    
    tempcoordstr = 's/TEMPCOORDFILE/'+new_coord_file+'/g'
    subprocess.call(["sed", "-i", "-e",  tempcoordstr, new_param_file])
    
    tempalphastr = 's/ALPHA/{:4.5s}/g'.format(str(lattice_param[3]))
    subprocess.call(["sed", "-i", "-e",  tempalphastr, new_param_file])
    
    tempbetastr = 's/BETA/{:4.5s}/g'.format(str(lattice_param[4]))
    subprocess.call(["sed", "-i", "-e",  tempbetastr, new_param_file])

    tempgammastr = 's/GAMMA/{:4.5s}/g'.format(str(lattice_param[5]))
    subprocess.call(["sed", "-i", "-e",  tempgammastr, new_param_file])
    
    tempaaastr = 's/AAA/{:4.5s}/g'.format(str(lattice_param[0]))
    subprocess.call(["sed", "-i", "-e",  tempaaastr, new_param_file])
    
    tempbbbstr = 's/BBB/{:4.5s}/g'.format(str(lattice_param[1]))
    subprocess.call(["sed", "-i", "-e",  tempbbbstr, new_param_file])

    tempcccstr = 's/CCC/{:4.5s}/g'.format(str(lattice_param[2]))
    subprocess.call(["sed", "-i", "-e",  tempcccstr, new_param_file])

    temprunstr = 's/TEMPRUNTYPE/GEO_OPT/g'
    subprocess.call(["sed", "-i", "-e",  temprunstr, new_param_file])


    subprocess.call(['mpirun', '-np', '112', 'cp2k.popt', '-i', INPUTFILE, '>', OUTPUTNAME])

    



##########################################
#                  TEST                  #
##########################################
def Return_Test_Coordinates():
    """
    This funciton returns coordinates for the test system
    Because there are no real coordiantes, it just returns a matrix of ones as a holder for the coordinates
    """
    coordinates = np.ones((1, 3))
    return coordinates


def Output_Test_New_Coordinate_File(lattice_parameters, Output):
    """
    This function saves the lattice parameters in a .npy file for the test Program
    """
    np.save(Output, lattice_parameters)


##########################################
#          Assistant Functions           #
##########################################
def Lattice_parameters_to_Crystal_matrix(lattice_parameters):
    """
    This function takes the lattice parameters and returns the crystal matrix

    **Required Inputs
    lattice_parameters = crystal lattice parameters as an array ([a,b,c,alpha,beta,gamma])
    """
    # Computing pieces of the crystal lattice matrix
    Vxx = lattice_parameters[0]
    Vxy = lattice_parameters[1]*np.cos(np.radians(lattice_parameters[5]))
    Vxz = lattice_parameters[2]*np.cos(np.radians(lattice_parameters[4]))

    Vyy = lattice_parameters[1]*np.sin(np.radians(lattice_parameters[5]))
    Vyz = lattice_parameters[2]*(np.cos(np.radians(lattice_parameters[3])) - np.cos(np.radians(lattice_parameters[4]))*np.cos(np.radians(lattice_parameters[5])))/np.sin(np.radians(lattice_parameters[5]))

    Vzz = np.sqrt(lattice_parameters[2]**2 - Vxz**2 - Vyz**2)
    # Combining the pieces of the matrix together
    crystal_matrix = np.matrix([[Vxx, Vxy, Vxz], [0., Vyy, Vyz], [0., 0., Vzz]])
    return crystal_matrix


def Crystal_matrix_to_Lattice_parameters(crystal_matrix):
    """
    This function takes a crystal lattice matrix and return the lattice parameters

    **Required Inputs
    crystal_matrix = crystal lattice matrix ([[Vxx,Vxy,Vxz],
                                              [0. ,Vyy,Vyz],
                                              [0. ,0. ,Vzz]])
    """
    # Computing lattice parameters
    a = np.sqrt(crystal_matrix[0, 0]**2)
    b = np.sqrt(crystal_matrix[0, 1]**2 + crystal_matrix[1, 1]**2)
    c = np.sqrt(crystal_matrix[0, 2]**2 + crystal_matrix[1, 2]**2 + crystal_matrix[2, 2]**2)
    gamma = np.degrees(np.arccos(crystal_matrix[0, 1]/b))
    beta = np.degrees(np.arccos(crystal_matrix[0, 2]/c))
    alpha = np.degrees(np.arccos(crystal_matrix[1,2]*np.sin(np.radians(gamma))/c + np.cos(np.radians(beta))*np.cos(np.radians(gamma))))

    # Assuring that the parameters returned are all positive
    if alpha < 0.0:
        alpha = 180. + alpha
    if beta < 0.0:
        beta = 180. + beta
    if gamma < 0.0:
        gamma = 180. + gamma

    # Creating an array of lattice parameters
    lattice_parameters = np.array([a, b, c, alpha, beta, gamma])
    return lattice_parameters


def Isotropic_Change_Lattice_Parameters(volume_fraction_change, Program, Coordinate_file):
    """
    This function returns the change in lattice parameters for isotropic expansion/compression based off of a given 
    change in volume fraction

    **Required Inputs
    volume_fraction_change = Volume of the new desired strucutre over the volume of the previous structure
    Program = 'Tinker' for Tinker Molecular Modeling
              'Test' for a test run
    Coordinate_file = file containing lattice parameters of the previous strucuture
    """
    # Calling the lattice parameters
    if Program == 'Tinker':
        lattice_parameters = Pr.Tinker_Lattice_Parameters(Coordinate_file)
    elif Program == 'Test':
        lattice_parameters = Pr.Test_Lattice_Parameters(Coordinate_file)
    elif Program == 'CP2K':
        lattice_parameters = Pr.CP2K_Lattice_Parameters(Coordinate_file)

    # Calculating the new isotropic lattice parameters
    dlattice_parameters = lattice_parameters*volume_fraction_change**(1/3.) - lattice_parameters

    # Setting changes in angles to zero because they do not change in isotropic expansion
    dlattice_parameters[3:] = 0.
    return dlattice_parameters


def Change_Crystal_Matrix(matrix_parameters_fraction_change, Program, Coordinate_file):
    """
    This function returns a matrix of changes in crystal matrix parameters based off of a general fraction
    This function is used for local gradient anisotropic expansion, each change in parameter is tested separately in 
    that case. We are combining them in one matrix here though

    **Required Inputs
    matrix_parameters_fraction_change = New parameter over old parameter
    Program = 'Tinker' for Tinker Molecular Modeling
              'Test' for a test run
    Coordinate_file = file containing lattice parameters of the previous structure
    """
    # Calling the lattice parameters
    if Program == 'Tinker':
        lattice_parameters = Pr.Tinker_Lattice_Parameters(Coordinate_file)
    elif Program == 'Test':
        lattice_parameters = Pr.Test_Lattice_Parameters(Coordinate_file)
    elif Program == 'CP2K':
        lattice_parameters = Pr.CP2K_Lattice_Parameters(Coordinate_file)

    # Computing the crystal matrix of the coordinate file and determining the change in parameters based off of the
    # fractional input
    dcrystal_matrix = Lattice_parameters_to_Crystal_matrix(lattice_parameters)*matrix_parameters_fraction_change

    # If lattice parameter is smaller than 5 in the coordinate file than the fractional change is increased as if it
    # was 5.
    # This is for cases where the angles are ~90 degrees, it lets us fully test the gradient due to temperature.
    for i in range(3):
        for j in np.arange(i+1, 3):
            if dcrystal_matrix[i, j] < matrix_parameters_fraction_change:
                dcrystal_matrix[i, j] = matrix_parameters_fraction_change
    return dcrystal_matrix
    

##########################################
#            General Expansion           #
##########################################
def Expand_Structure(Coordinate_file, Program, Expansion_type, molecules_in_coord, Output, min_RMS_gradient, **keyword_parameters):
    """
    This function expands a coordinate file either based off of an inputted change in lattice vectors or crystal 
        lattice matrix

    **Required Inputs
    Coordinate_file = file containing lattice parameters (and coordinates)
    Program = 'Tinker' for Tinker Molecular Modeling
              'Test' for a test run
    Expansion_type = 'lattice parameters' expanding the structure by lattice parameters ([a,b,c,alpha,beta,gamma])
                   = 'crystal_matrix' expanding the strucutre by changes in the crystal matrix
    molecules_in_coord = number of molecules in the coordinate file
    Output = string to name expanded coordinate file
    
    **Optional Inputs
    dlattice_vectors = Changes in lattice parameters
    dcrystal_matrix = changes in crystal matrix
    Parameter_file = program specific file containingforce field parameters
    """
    if Program == 'Test':
        lattice_parameters = Pr.Test_Lattice_Parameters(Coordinate_file)
        if Expansion_type == 'lattice_parameters':
            lattice_parameters = lattice_parameters + keyword_parameters['dlattice_parameters']
        elif Expansion_type == 'crystal_matrix':
            crystal_matrix = Lattice_parameters_to_Crystal_matrix(lattice_parameters)
            crystal_matrix = crystal_matrix + keyword_parameters['dcrystal_matrix']
            lattice_parameters = Crystal_matrix_to_Lattice_parameters(crystal_matrix)
        Output_Test_New_Coordinate_File(lattice_parameters, Output)

    else:
        if Program == 'Tinker':
            coordinates = Return_Tinker_Coordinates(Coordinate_file)
            lattice_parameters = Pr.Tinker_Lattice_Parameters(Coordinate_file)
	elif Program == 'CP2K':
	    coordinates = Return_CP2K_Coordinates(Coordinate_file)
	    lattice_parameters = Pr.CP2K_Lattice_Parameters(Coordinate_file)

        crystal_matrix = Lattice_parameters_to_Crystal_matrix(lattice_parameters)

        Coordinate_center_of_mass = np.zeros((molecules_in_coord, 3))
        atoms_per_molecule = len(coordinates[:, 0])/molecules_in_coord

        for i in range(int(molecules_in_coord)):
            Coordinate_center_of_mass[i, :] = np.mean(coordinates[i*atoms_per_molecule:(i+1)*atoms_per_molecule],
                                                      axis=0)
            coordinates[i*atoms_per_molecule:(i+1)*atoms_per_molecule] = \
                np.subtract(coordinates[i*atoms_per_molecule:(i+1)*atoms_per_molecule], Coordinate_center_of_mass[i, :])

        # Center of mass coordinates converted to fractional coordinates
        Coordinate_center_of_mass = np.dot(np.linalg.inv(crystal_matrix), Coordinate_center_of_mass.T).T
        if Expansion_type == 'lattice_parameters':
            lattice_parameters = lattice_parameters + keyword_parameters['dlattice_parameters']
            crystal_matrix = Lattice_parameters_to_Crystal_matrix(lattice_parameters)
        elif Expansion_type == 'crystal_matrix':
            crystal_matrix = crystal_matrix + keyword_parameters['dcrystal_matrix']
            lattice_parameters = Crystal_matrix_to_Lattice_parameters(crystal_matrix)
        Coordinate_center_of_mass = np.dot(crystal_matrix, Coordinate_center_of_mass.T).T
        for i in range(int(molecules_in_coord)):
            coordinates[i*atoms_per_molecule:(i+1)*atoms_per_molecule] = \
                np.subtract(coordinates[i*atoms_per_molecule:(i+1)*atoms_per_molecule],
                            -1*Coordinate_center_of_mass[i, :])

        if Program == 'Tinker':
            Output_Tinker_New_Coordinate_File(Coordinate_file, keyword_parameters['Parameter_file'], coordinates,
                                              lattice_parameters, Output, min_RMS_gradient)
        elif Program == 'cp2k':
            Output_CP2K_New_Coordinate_File(Coordinate_file, keyword_parameters['Parameter_file'], coordinates,
                                              lattice_parameters, Output, min_RMS_gradient)

##########################################
#       Local Gradient of Expansion      #
##########################################
def Isotropic_Local_Gradient(Coordinate_file, Program, Temperature, Pressure, LocGrd_Vol_FracStep,
                             molecules_in_coord, Statistical_mechanics, Method, min_RMS_gradient, **keyword_parameters):
    """
    This function calculates the local gradient of isotropic expansion for a given coordinate file
    
    :param Coordinate_file: file containing lattice parameters (and coordinates)
    :param Program: 'Tinker' Tinker molecular modeling
                    'Test' Test case
    :param Temperature: in Kelvin
    :param Pressure: in atm
    :param LocGrd_Vol_FracStep: fractional volumetric step size for numerical gradient 
    :param molecules_in_coord: number of molecules in Coordinate_file
    :param Statistical_mechanics: 'Classical' Classical mechanics
                                  'Quantum' Quantum mechanics
    :param Method: 'GiQ' Gradient isotropic QHA
                   'GiQg' Gradient isotropic QHA with Gruneisen Parameter
    :param keyword_parameters: Parameter_file, Gruneisen, Wavenumber_reference, Volume_reference
    
    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    Gruneisen: isotropic Gruneisen parameter
    Wavenumber_reference: reference wavenumbers for the Gruneisen parameter
    Volume_reference: reference volume for the Gruneisen parameter
    """
    # Assigning general names for expanded and compressed structures
    if Program == 'Tinker':
        coordinate_plus = 'plus.xyz'
        coordinate_minus = 'minus.xyz'
    elif Program == 'Test':
        coordinate_plus = 'plus.npy'
        coordinate_minus = 'minus.npy'
        keyword_parameters['Parameter_file'] = ''
    elif Program == 'CP2K':
        coordinate_plus = 'plus.pdb'
        coordinate_minus = 'minus.pdb'
        keyword_parameters['Parameter_file'] = ''

    # Determining the change in lattice parameter for isotropic expansion
    dlattice_parameters = Isotropic_Change_Lattice_Parameters(1 + LocGrd_Vol_FracStep, Program, Coordinate_file)

    # Building the isotropically expanded and compressed strucutres
    Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'plus', min_RMS_gradient,
                     dlattice_parameters=dlattice_parameters, Parameter_file=keyword_parameters['Parameter_file'])
    Expand_Structure(Coordinate_file, Program, 'lattice_parameters', molecules_in_coord, 'minus', min_RMS_gradient,
                     dlattice_parameters=-1*dlattice_parameters, Parameter_file=keyword_parameters['Parameter_file'])

    # Determining the volume of Coordinate_file
    volume = Pr.Volume(Program=Program, Coordinate_file=Coordinate_file)
    # Calculating wavenumbers coordinate_file, plus.*, and minus.*
    if Method == 'GiQ':
        wavenumbers = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file=Coordinate_file,
                                           Parameter_file=keyword_parameters['Parameter_file'],
                                           Program=Program)
        wavenumbers_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file=coordinate_plus,
                                                Parameter_file=keyword_parameters['Parameter_file'],
                                                Program=Program)
        wavenumbers_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file=coordinate_minus,
                                                 Parameter_file=keyword_parameters['Parameter_file'],
                                                 Program=Program)
    else:
        wavenumbers = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                           Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                           Volume_Reference=keyword_parameters['Volume_Reference'],
                                           New_Volume=volume)
        wavenumbers_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                                Volume_Reference=keyword_parameters['Volume_Reference'],
                                                New_Volume=volume + volume*LocGrd_Vol_FracStep)
        wavenumbers_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                 Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                                 Volume_Reference=keyword_parameters['Volume_Reference'],
                                                 New_Volume=volume - volume*LocGrd_Vol_FracStep)

    # If temperature is zero, we assume that the local gradient is the same at 0.1K
    if Temperature == 0.:
        Temperature = 0.1

    # Checking if wavenumbers are largely negative
    if any(wavenumbers < -1.) or any(wavenumbers_plus < -1.) or any (wavenumbers_minus < -1.):
        print "   ... Check wavenumbers:"
        print "   ...... Wavenumbers:" + str(wavenumbers[:3])
        print "   ...... Expanded:" + str(wavenumbers_plus[:3])
        print "   ...... Compressed:" + str(wavenumbers_minus[:3])

    # Calculating the numerator of the local gradient -dS/dV
    numerator = -(Pr.Vibrational_Entropy(Temperature, wavenumbers_plus, Statistical_mechanics)/molecules_in_coord -
                  Pr.Vibrational_Entropy(Temperature, wavenumbers_minus, Statistical_mechanics)/molecules_in_coord) / \
                 (2*LocGrd_Vol_FracStep*volume)

    # Calculating the denominator of the local gradient d**2G/dV**2
    denominator = (Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_plus, coordinate_plus,
                                        Statistical_mechanics, molecules_in_coord,
                                        Parameter_file=keyword_parameters['Parameter_file']) -
                   2*Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers, Coordinate_file,
                                          Statistical_mechanics, molecules_in_coord,
                                          Parameter_file=keyword_parameters['Parameter_file']) +
                   Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_minus, coordinate_minus,
                                        Statistical_mechanics, molecules_in_coord,
                                        Parameter_file=keyword_parameters['Parameter_file'])) / \
                  ((volume*LocGrd_Vol_FracStep)**2)

    # Removing excess files
    os.system('rm '+coordinate_plus+' '+coordinate_minus)
    return -numerator/denominator, wavenumbers, volume



def Anisotropic_Local_Gradient(Coordinate_file, Program, Temperature, Pressure, LocGrd_LatParam_FracStep,
                               molecules_in_coord, Statistical_mechanics, Method,
                               Hessian_number, min_RMS_gradient, **keyword_parameters):
    """
    This function calculates the local gradient of anisotropic expansion for a given coordinate file
    
    :param Coordinate_file: file containing lattice parameters (and coordinates)
    :param Program: 'Tinker' Tinker molecular modeling
                    'Test' Test case
    :param Temperature: in Kelvin
    :param Pressure: in atm
    :param LocGrd_LatParam_FracStep: 
    :param molecules_in_coord: number of molecules in Coordinate_file
    :param Statistical_mechanics: 'Classical' Classical mechanics
                                  'Quantum' Quantum mechanics
    :param Method: 'GaQ' Gradient anisotropic QHA
                   'GaQg' Gradient anisotropic QHA with Gruneisen Parameter
    :param Hessian_number: 73 Hessians to calculate the complete anistropic gradient
                           25 for d**2G_dUdU only calculating the diagonals and off-diags. of the upper left 3x3 matrix
                           19 for d**2G_dUdU only calculating the upper left 3x3 matrix
                           13 for d**2G_dUdU only calculating the diagonals
                           7  for d**2G_dUdU only calculating the upper left 3x3 matrix daigonals
    :param keyword_parameters: Parameter_file, Gruneisen, Wavenumber_reference, Volume_reference

    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    Gruneisen: isotropic Gruneisen parameter
    Wavenumber_reference: reference wavenumbers for the Gruneisen parameter
    Crystal_matrix_Reference
    """
    # Determining the file ending of the coordinate files
    if Program == 'Tinker':
        file_ending = '.xyz'
    elif Program == 'CP2K':
        file_ending = '.pdb'
    elif Program == 'Test':
        file_ending = '.npy'
        keyword_parameters['Parameter_file'] = ''

    # Preparing the matrix with each entry as d**2G/(du*du')
    dG_U = np.zeros((6, 6))
    # Preparing the vector with each entry as d*G/du*dT
    dS_dU = np.zeros(6)
   
    # Modified anisotropic Local Gradient
    if Hessian_number == 73:
        diag_limit = 6
        off_diag_limit = 6
    elif Hessian_number == 25:
        diag_limit = 6
        off_diag_limit = 3
    elif Hessian_number == 19:
        diag_limit = 3
        off_diag_limit = 3
    elif Hessian_number == 13:
        diag_limit = 6
        off_diag_limit = 0
    elif Hessian_number == 7:
        diag_limit = 3
        off_diag_limit = 0

    # Determining the changes in crystal lattice parameters
    dcrystal_matrix_hold = Change_Crystal_Matrix(LocGrd_LatParam_FracStep, Program, Coordinate_file)

    # Setting the order of lattice parameter to perform (diagonals first followed by off diagonals)
    matrix_order = np.matrix([[0, 0], [1, 1], [2, 2], [0, 1], [0, 2], [1, 2]])

    # Retrieving the wavenumbers and crystal matrix of the initial structure
    if Method == 'GaQ':
        if Program == 'Tinker':
            wavenumbers = Wvn.Tinker_Wavenumber(Coordinate_file, Parameter_file=keyword_parameters['Parameter_file'])
        if Program == 'CP2K':
            wavenumbers = Wvn.CP2K_Wavenumber(Coordinate_file, Parameter_file=keyword_parameters['Parameter_file'], cp2kroot = keyword_parameters['Parameter_file'])
        elif Program == 'Test':
            wavenumbers = Wvn.Test_Wavenumber(Coordinate_file)
    elif Method == 'GaQg':
        if Program == 'Tinker':
            crystal_matrix = Lattice_parameters_to_Crystal_matrix(Pr.Tinker_Lattice_Parameters(Coordinate_file))
        if Program == 'CP2K':
            crystal_matrix = Lattice_parameters_to_Crystal_matrix(Pr.CP2K_Lattice_Parameters(Coordinate_file))
        elif Program == 'Test':
            crystal_matrix = Lattice_parameters_to_Crystal_matrix(Pr.Test_Lattice_Parameters(Coordinate_file))
        wavenumbers = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                           Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                           Crystal_matrix_Reference=keyword_parameters['Crystal_matrix_Reference'],
                                           New_Crystal_matrix=crystal_matrix)

    # If temperature is zero, we assume that the local gradient is the same at 0.1K
    if Temperature == 0.:
        Temperature = 0.1

    # Calculating the diagonals of dG_U and the vector dS_dU
    for i in range(diag_limit):
        dcrystal_matrix = np.zeros((3, 3))
        lattice_stepsize = dcrystal_matrix_hold[matrix_order[i, 0], matrix_order[i, 1]]
        # Building the change in the crystal matrix for the lattice parameter at hand
        dcrystal_matrix[matrix_order[i, 0], matrix_order[i, 1]] = lattice_stepsize

        # Finding the structures by increasing and decreasing the lattice parameter
        Expand_Structure(Coordinate_file, Program, 'crystal_matrix', molecules_in_coord, 'p', min_RMS_gradient,
                         dcrystal_matrix=dcrystal_matrix, Parameter_file=keyword_parameters['Parameter_file'])
        Expand_Structure(Coordinate_file, Program, 'crystal_matrix', molecules_in_coord, 'm', min_RMS_gradient,
                         dcrystal_matrix=-1*dcrystal_matrix, Parameter_file=keyword_parameters['Parameter_file'])

        # Calculating the wavenumbers for the additional strucutre
        if Method == 'GaQ':
            wavenumbers_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file='p'+file_ending, Program=Program,
                                                    Parameter_file=keyword_parameters['Parameter_file'])
            wavenumbers_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file='m'+file_ending, Program=Program,
                                                     Parameter_file=keyword_parameters['Parameter_file'])
        elif Method == 'GaQg':
            wavenumbers_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                    Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                                    Crystal_matrix_Reference=
                                                    keyword_parameters['Crystal_matrix_Reference'],
                                                    New_Crystal_matrix=crystal_matrix + dcrystal_matrix)
            wavenumbers_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                     Wavenumber_Reference=keyword_parameters['Wavenumber_Reference'],
                                                     Crystal_matrix_Reference=
                                                     keyword_parameters['Crystal_matrix_Reference'],
                                                     New_Crystal_matrix=crystal_matrix - dcrystal_matrix)

        dS_dU[i] = (Pr.Vibrational_Entropy(Temperature, wavenumbers_plus, Statistical_mechanics) / molecules_in_coord -
                     Pr.Vibrational_Entropy(Temperature, wavenumbers_minus, Statistical_mechanics) /
                     molecules_in_coord) / (2 * lattice_stepsize)

#        if np.absolute(dS_dU[i]) < 9.e-06:
#            dS_dU[i] = 0.

#        dG = (Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_plus, 'p'+file_ending,
#                                           Statistical_mechanics, molecules_in_coord,
#                                           Parameter_file=keyword_parameters['Parameter_file']) -
#                      Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_minus, 'm'+file_ending,
#                                           Statistical_mechanics, molecules_in_coord,
#                                           Parameter_file=keyword_parameters['Parameter_file']))/(2*lattice_stepsize)

#        if np.absolute(dG) > 9.e-06:
        dG_U[i, i] = (Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_plus, 'p'+file_ending,
                                           Statistical_mechanics, molecules_in_coord,
                                           Parameter_file=keyword_parameters['Parameter_file']) -
                      2*Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers, Coordinate_file,
                                             Statistical_mechanics, molecules_in_coord,
                                             Parameter_file=keyword_parameters['Parameter_file']) +
                      Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_minus, 'm'+file_ending,
                                           Statistical_mechanics, molecules_in_coord,
                                           Parameter_file=keyword_parameters['Parameter_file']))/(lattice_stepsize**2)

        # Calculating the off-diagonals
        if i < off_diag_limit:
            for j in np.arange(i + 1, off_diag_limit):
                dcrystal_matrix_2 = np.zeros((3, 3))
                lattice_stepsize_2 = dcrystal_matrix_hold[matrix_order[j, 0], matrix_order[j, 1]]
                dcrystal_matrix_2[matrix_order[j, 0], matrix_order[j, 1]] = lattice_stepsize_2

                # Expanding the structures along a different crystal matrix parameter
                Expand_Structure('p'+file_ending, Program, 'crystal_matrix', molecules_in_coord, 'pp', min_RMS_gradient,
                                 dcrystal_matrix=dcrystal_matrix_2, Parameter_file=keyword_parameters['Parameter_file'])
                Expand_Structure('p'+file_ending, Program, 'crystal_matrix', molecules_in_coord, 'pm', min_RMS_gradient,
                                 dcrystal_matrix=-1*dcrystal_matrix_2,
                                 Parameter_file=keyword_parameters['Parameter_file'])
                Expand_Structure('m'+file_ending, Program, 'crystal_matrix', molecules_in_coord, 'mm', min_RMS_gradient,
                                 dcrystal_matrix=-1*dcrystal_matrix_2,
                                 Parameter_file=keyword_parameters['Parameter_file'])
                Expand_Structure('m'+file_ending, Program, 'crystal_matrix', molecules_in_coord, 'mp', min_RMS_gradient,
                                 dcrystal_matrix=dcrystal_matrix_2, Parameter_file=keyword_parameters['Parameter_file'])

                # Calculating the wavenumbers of the new strucuture
                if Method == 'GaQ':
                    wavenumbers_plus_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file='pp' + file_ending, 
                                                                 Program=Program, 
                                                                 Parameter_file=keyword_parameters['Parameter_file'])
                    wavenumbers_plus_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file='pm' + file_ending,
                                                                  Program=Program,
                                                                  Parameter_file=keyword_parameters['Parameter_file'])
                    wavenumbers_minus_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file='mm' + file_ending,
                                                                   Program=Program,
                                                                   Parameter_file=keyword_parameters['Parameter_file'])
                    wavenumbers_minus_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Coordinate_file='mp' + file_ending,
                                                                  Program=Program,
                                                                  Parameter_file=keyword_parameters['Parameter_file'])
                elif Method == 'GaQg':
                    wavenumbers_plus_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                                 Wavenumber_Reference=
                                                                 keyword_parameters['Wavenumber_Reference'],
                                                                 Crystal_matrix_Reference=
                                                                 keyword_parameters['Crystal_matrix_Reference'],
                                                                 New_Crystal_matrix=
                                                                 crystal_matrix + dcrystal_matrix + dcrystal_matrix_2)
                    wavenumbers_plus_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                                  Wavenumber_Reference=
                                                                  keyword_parameters['Wavenumber_Reference'],
                                                                  Crystal_matrix_Reference=
                                                                  keyword_parameters['Crystal_matrix_Reference'],
                                                                  New_Crystal_matrix=
                                                                  crystal_matrix + dcrystal_matrix - dcrystal_matrix_2)
                    wavenumbers_minus_minus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                                   Wavenumber_Reference=
                                                                   keyword_parameters['Wavenumber_Reference'],
                                                                   Crystal_matrix_Reference=
                                                                   keyword_parameters['Crystal_matrix_Reference'],
                                                                   New_Crystal_matrix=
                                                                   crystal_matrix - dcrystal_matrix - dcrystal_matrix_2)
                    wavenumbers_minus_plus = Wvn.Call_Wavenumbers(Method, min_RMS_gradient, Gruneisen=keyword_parameters['Gruneisen'],
                                                                  Wavenumber_Reference=
                                                                  keyword_parameters['Wavenumber_Reference'],
                                                                  Crystal_matrix_Reference=
                                                                  keyword_parameters['Crystal_matrix_Reference'],
                                                                  New_Crystal_matrix=
                                                                  crystal_matrix - dcrystal_matrix + dcrystal_matrix_2)

                # Calculating the diagonal elements of dG_U
                dG_U[i, j] = (Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_plus_plus,
                                                   'pp'+file_ending, Statistical_mechanics, molecules_in_coord,
                                                   Parameter_file=keyword_parameters['Parameter_file']) -
                              Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_plus_minus, 'pm' +
                                                   file_ending, Statistical_mechanics, molecules_in_coord,
                                                   Parameter_file=keyword_parameters['Parameter_file']) -
                              Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_minus_plus, 'mp' +
                                                   file_ending, Statistical_mechanics, molecules_in_coord,
                                                   Parameter_file=keyword_parameters['Parameter_file']) +
                              Pr.Gibbs_Free_Energy(Temperature, Pressure, Program, wavenumbers_minus_minus, 'mm' +
                                                   file_ending, Statistical_mechanics, molecules_in_coord,
                                                   Parameter_file=keyword_parameters['Parameter_file']
                                                   ))/(4*lattice_stepsize*lattice_stepsize_2)
#                if dG_U[i, j] < 9.e-06:
#                    dG_U[i, j] = 0.

                dG_U[j, i] = dG_U[i, j]
                # Removing excess files
                os.system('rm pm' + file_ending + ' mp' + file_ending + ' pp' + file_ending + ' mm' + file_ending)
        # Removing excess files
        os.system('rm p' + file_ending + ' m' + file_ending)

    # Putting together the total matrix for the changes in lattice parameters with temperatures
#    dCrystal_Matrix = np.linalg.solve(dG_U, dS_dU)
    dCrystal_Matrix = np.linalg.lstsq(dG_U, dS_dU)[0]
    for i in range(len(dCrystal_Matrix)):
        if np.absolute(dCrystal_Matrix[i]) < 9e-05:
            dCrystal_Matrix[i] = 0.

    dCrystal_Matrix = np.matrix([[dCrystal_Matrix[0], dCrystal_Matrix[3], dCrystal_Matrix[4]],
                                 [0., dCrystal_Matrix[1], dCrystal_Matrix[5]], 
                                 [0., 0., dCrystal_Matrix[2]]])


    return dCrystal_Matrix, wavenumbers
