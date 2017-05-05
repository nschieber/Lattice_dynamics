#!/usr/bin/env python
import os
import sys
import numpy as np
import Expand as Ex
import Properties as Pr
import Wavenumbers as Wvn

##########################################
#           Numerical Methods            #
##########################################
def Runge_Kutta_Fourth_Order(Method, Coordinate_file, Program, Temperature, Pressure, LocGrd_Temp_step,
                             molecules_in_coord, Statistical_mechanics, RK4_stepsize, **keyword_parameters):
    """
    This function determines the gradient of thermal expansion of a strucutre between two temperatures using
    a forth order Runge-Kutta numerical analysis
    :param Method: Gradient Isotropic QHA ('GiQ');
                   Gradient Isotropic QHA w/ Gruneisen Parameter ('GiQg');
                   Gradient Anisotropic QHA ('GaQ');
    :param Coordinate_file: file containing the lattice parameters (and coordinates)
    :param Program: 'tink' for Tinker Molecular Modeling
                    'test' for a test run
    :param Temperature: in Kelvin
    :param Pressure: in atm
    :param LocGrd_Temp_step: temperature step size to use in local gradient 
    :param molecules_in_coord: number of molecules in coordinate file
    :param Statistical_mechanics: 'C' Classical mechanics
                                  'Q' Quantum mechanics
    :param RK4_stepsize: stepsize for runge-kutta 4th order
    :param keyword_parameters: Parameter_file, LocGrd_Vol_FracStep, LocGrd_LatParam_FracStep, Gruneisen, 
    Wavenumber_Reference, Volume_Reference, Aniso_LocGrad_Type
    
    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    LocGrd_Vol_FracStep: isotropic volume fractional stepsize for local gradient
    LocGrd_LatParam_FracStep: anisotropic crystal matrix fractional stepsize for local gradient
    Gruneisen: Gruneisen parameters found with Setup_Isotropic_Gruneisen
    Wavenumber_Reference: Reference wavenumbers for Gruneisen parameter
    Volume_Reference: Reference volume of structure for Wavenumber_Reference
    Aniso_LocGrad_Type: 73 Hessians to calculate the complete anistropic gradient
                        25 for d**2G_dUdU only calculating the diagonals and off-diags. of the upper left 3x3 matrix
                        19 for d**2G_dUdU only calculating the uppder left 3x3 matrix
                        13 for d**2G_dUdU only calculating the diagonals
                        7  for d**2G_dUdU only calculating the upper left 3x3 matrix daigonals
    """
    # Setting up program specific file endings and giving parameter files blank names to avoid errors
    if Program == 'tink':
        file_ending = '.xyz'
    elif Program == 'test':
        file_ending = '.npy'
        keyword_parameters['Parameter_file'] = ''

    RK_multiply = np.array([1./6.,1./3.,1./3.,1./6.])

    # Copying the coordinate file to a seperate file to work with
    os.system('cp ' + Coordinate_file + ' RK4' + file_ending)

    # Setting the different temperature stepsizes
    temperature_steps = np.array([0., RK4_stepsize/2., RK4_stepsize/2., RK4_stepsize])

    # Setting RK_4 array/matix and general parameters that aren't required for specific methods
    if (Method == 'GiQ') or (Method == 'GiQg'):
        RK_grad = np.zeros(4)
        if Method == 'GiQ':
            keyword_parameters['Gruneisen'] = 0.
            keyword_parameters['Wavenumber_Reference'] = 0.
            keyword_parameters['Volume_Reference'] = 0.
    elif (Method == 'GaQ') or (Method == 'GaQg'):
        RK_grad = np.zeros((4, 3, 3))
        if Method == 'GaQ':
            keyword_parameters['Gruneisen'] = 0.
            keyword_parameters['Wavenumber_Reference'] = 0.
            keyword_parameters['Volume_Reference'] = 0.

    # Calculating the RK gradients for the overall numerical gradient
    for i in range(4):
        if (Method == 'GiQ') or (Method == 'GiQg'):
            RK_grad[i], wavenumbers_hold, volume_hold = Ex.Call_Expansion(Method, 'local_gradient', Program,
                                                                         'RK4' + file_ending, molecules_in_coord,
                                                                         Temperature=Temperature, Pressure=Pressure,
                                                                         volume_fraction_change=keyword_parameters[
                                                                              'LocGrd_Vol_FracStep'],
                                                                         LocGrd_Temp_step=LocGrd_Temp_step,
                                                                         Statistical_mechanics=Statistical_mechanics,
                                                                         Parameter_file=
                                                                         keyword_parameters['Parameter_file'],
                                                                         Gruneisen=keyword_parameters['Gruneisen'],
                                                                         Wavenumber_Reference=
                                                                         keyword_parameters['Wavenumber_Reference'],
                                                                         Volume_Reference=
                                                                         keyword_parameters['Volume_Reference'])
        elif (Method == 'GaQ') or (Method == 'GaQg'):
            RK_grad[i], wavenumbers_hold = Ex.Call_Expansion(Method, 'local_gradient', Program, 'RK4' + file_ending,
                                                            molecules_in_coord, Temperature=Temperature,
                                                            Pressure=Pressure, matrix_parameters_fraction_change=
                                                            keyword_parameters['LocGrd_LatParam_FracStep'],
                                                            LocGrd_Temp_step= LocGrd_Temp_step,
                                                            Statistical_mechanics=Statistical_mechanics,
                                                            Parameter_file=keyword_parameters['Parameter_file'],
                                                            Gruneisen=keyword_parameters['Gruneisen'],
                                                            Wavenumber_Reference=
                                                            keyword_parameters['Wavenumber_Reference'],
                                                            Volume_Reference=keyword_parameters['Volume_Reference'],
                                                            Aniso_LocGrad_Type=keyword_parameters['Aniso_LocGrad_Type'])
            volume_hold = 0.
        if i == 0:
            wavenumbers = wavenumbers_hold
            volume = volume_hold
        if i != 3:
            if (Method == 'GiQ') or (Method == 'GiQg'):
                dcrystal_matrix = 0.
                volume_fraction_change = (volume + RK_grad[i]*temperature_steps[i+1])/volume
            elif (Method == 'GaQ') or (Method == 'GaQg'):
                dcrystal_matrix = RK_grad[i]*temperature_steps[i+1]
                volume_fraction_change = 0.
            dvolume_or_matrix = RK_grad[i]*temperature_steps[i+1]
            # Expanding the crystal to the next stepsize
            Ex.Call_Expansion(Method, 'expand', Program, Coordinate_file, molecules_in_coord,
                              Parameter_file=keyword_parameters['Parameter_file'], dcrystal_matrix=dcrystal_matrix,
                              volume_fraction_change=volume_fraction_change, Output='RK4')
        # Multiplying the found gradient by the fraction it will contribute to the overall gradient
        RK_grad[i] = RK_grad[i]*RK_multiply[i]

    # Summing all RK gradients for the overall numerical gradient
    numerical_gradient = np.sum(RK_grad, axis=0)

    # Removign excess files
    os.system('rm RK4'+file_ending)
    return numerical_gradient, wavenumbers, volume


##########################################
#      Stepwise Isotropic Expansion      #
##########################################
def Isotropic_Stepwise_Expansion(StepWise_Vol_StepFrac, StepWise_Vol_LowerFrac, StepWise_Vol_UpperFrac, Coordinate_file,
                                 Program, Temperature, Pressure, Output, Method, molecules_in_coord, Wavenum_Tol,
                                 Statistical_mechanics, **keyword_parameters):
    """
    This function performs stepwise isotropic QHA either with or without the gruneisen parameter
    :param StepWise_Vol_StepFrac: volumetric fraction step
    :param StepWise_Vol_LowerFrac: lower bound on the fraction to compress to
    :param StepWise_Vol_UpperFrac: uppder gound on the fraction to expand to
    :param Coordinate_file: file containing the lattice parameters (and coordinates)
    :param Program: 'tink' for Tinker Molecular Modeling
                    'test' for a test run
    :param Temperature: array of temperatures in Kelvin
    :param Pressure: in atm
    :param Output: string for outputted files
    :param Method: Stepwise Isotropic QHA ('SiQ');
                   Stepwise Isotropic QHA w/ Gruneisen Parameter ('SiQg');
    :param molecules_in_coord: number of molecules in coordinate file
    :param Wavenum_Tol: lowest tollerable wavenumbers (some small negative wavenumbers are okay)
    :param Statistical_mechanics: 'C' Classical mechanics
                                  'Q' Quantum mechanics
    :param keyword_parameters: Parameter_file, Gruneisen_Vol_FracStep
    
    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    Gruneisen_Vol_FracStep: volume fraction step used to determine the Gruneisen parameter 
    """
    # Setting file endings and determining how many wavenumbers there will be
    if Program == 'tink':
        file_ending = '.xyz'
        number_of_wavenumbers = Pr.Tinker_atoms_per_molecule(Coordinate_file, 1)*3
    elif Program == 'test':
        file_ending = '.npy'
        number_of_wavenumbers = len(Wvn.Test_Wavenumber(Coordinate_file))
        keyword_parameters['Parameter_file'] = ''

    # Setting up array of volume fractions from the lattice structure
    volume_fraction = np.arange(StepWise_Vol_LowerFrac, StepWise_Vol_UpperFrac, StepWise_Vol_StepFrac)
    volume_fraction = np.append(volume_fraction[:(1.0 - StepWise_Vol_UpperFrac)/StepWise_Vol_StepFrac:][::-1],
                                volume_fraction[(1.0 - StepWise_Vol_LowerFrac)/StepWise_Vol_StepFrac:])

    # Setting up a matrix to store the wavenumbers in
    wavenumbers = np.zeros((len(volume_fraction), number_of_wavenumbers+1))
    wavenumbers[:, 0] = volume_fraction

    # Setting parameters for the Gruneisen parameter and loading in previously found wavenumbers for SiQ
    if Method == 'SiQg':
        Gruneisen, Wavenumber_Reference, Volume_Reference = Wvn.Call_Wavenumbers(Method,Coordinate_file=Coordinate_file,
                                                                                 Program=Program,
                                                                                 Gruneisen_Vol_FracStep=
                                                                                 keyword_parameters[
                                                                                     'Gruneisen_Vol_FracStep'],
                                                                                 molecules_in_coord=molecules_in_coord,
                                                                                 Parameter_file=
                                                                                 keyword_parameters['Parameter_file'])
    elif Method == 'SiQ':
        Gruneisen = 0.
        Wavenumber_Reference = 0.
        Volume_Reference = 0.
        if os.path.isfile(Output + '_' + Method + '_WVN.npy'):
            wavenumbers = np.append(wavenumbers, np.load(Output + '_' + Method + '_WVN.npy'), axis=0)

    # setting a matrix for properties versus temperature and pressure
    properties = np.zeros((len(volume_fraction), len(Temperature), 14))

    # Finding all expanded structures
    previous_volume = 1.0
    lattice_volume = Pr.Volume(Program=Program, Coordinate_file=Coordinate_file)
    os.system('cp ' + Coordinate_file + ' ' + Output + '_' + Method + str(previous_volume) + file_ending)
    for i in range(len(volume_fraction)):
        if os.path.isfile('Cords/' + Output + '_' + Method + str(volume_fraction[i]) + file_ending):
            # Skipping structures if they've already been constructed
            os.system('cp Cords/' + Output + '_' + Method + str(volume_fraction[i]) + file_ending + ' ./')
        else:
            Ex.Call_Expansion(Method, 'expand', Program, Output + '_' + Method + str(previous_volume) + file_ending,
                              molecules_in_coord, Parameter_file=keyword_parameters['Parameter_file'],
                              volume_fraction_change=(volume_fraction[i]/previous_volume),
                              Output=Output + '_' + Method + str(volume_fraction[i]))
        # Calculating wavenumbers of new expanded strucutre
        wavenumbers[i, 1:] = Wvn.Call_Wavenumbers(Method, Program=Program, Gruneisen=Gruneisen,
                                                  Wavenumber_Reference=Wavenumber_Reference,
                                                  Volume_Reference=Volume_Reference,
                                                  New_Volume=volume_fraction[i]*lattice_volume,
                                                  Coordinate_file=(Output + '_' + Method +
                                                                   str(volume_fraction[i]) + file_ending),
                                                  Parameter_file=keyword_parameters['Parameter_file'])

        # Saving the wavenumbers if the Gruneisen parameter is not being used
        if Method == 'SiQ':
            np.save(Output + '_' + Method + '_WVN', wavenumbers[~np.all(wavenumbers == 0, axis=1)])

        # Calculating properties of systems with wavenumbers above user specified tollerance
        if all(wavenumbers[i, 1:] > Wavenum_Tol):
            properties[i, :, :] = Pr.Properties_with_Temperature(Output + '_' + Method + str(volume_fraction[i]) +
                                                                 file_ending, wavenumbers[i, 1:], Temperature, Pressure,
                                                                 Program, Statistical_mechanics, molecules_in_coord,
                                                                 Parameter_file=keyword_parameters['Parameter_file'])
        else:
            properties[i, :, :] = np.nan

        # Moving old strucutres around and adjusting parameters to find the next strucutre
        os.system('mv ' + Output + '_' + Method + str(previous_volume) + file_ending + ' Cords/')
        previous_volume = volume_fraction[i]
        if volume_fraction[i] == StepWise_Vol_LowerFrac:
            os.system('cp ' + Output + '_' + Method + str(previous_volume) + file_ending + ' Cords/')
            previous_volume = 1.0
            os.system('cp Cords/' + Output + '_' + Method + str(previous_volume) + file_ending + ' ./')
        if volume_fraction[i] == max(volume_fraction):
            os.system('mv ' + Output + '_' + Method + str(previous_volume) + file_ending + ' Cords/')

    # Saving the raw data before minimizing
    np.save(Output+'_raw', properties)

    # Building matrix for minimum Gibbs Free energy data across temperature range
    minimum_gibbs_properties = np.zeros((len(Temperature), 14))
    for i in range(len(properties[0, :, 0])):
        for j in range(len(properties[:, 0, 0])):
            if properties[j, i, 2] == np.nanmin(properties[:, i, 2]):
                minimum_gibbs_properties[i, :] = properties[j, i, :]
    return minimum_gibbs_properties

##########################################
#      Gradient Isotropic Expansion      #
##########################################
def Isotropic_Gradient_Expansion(Coordinate_file, Program, molecules_in_coord, Output, Method, Gradient_MaxTemp,
                                 Pressure, LocGrd_Vol_FracStep, LocGrd_Temp_step, Statistical_mechanics,
                                 NumAnalysis_step, NumAnalysis_method, **keyword_parameters):
    """
    This function calculated the isotropic gradient for thermal expansion and returns the properties along that path
    :param Coordinate_file: file containing the lattice parameters (and coordinates)
    :param Program: 'tink' for Tinker Molecular Modeling
                    'test' for a test run
    :param molecules_in_coord: number of molecules in coordinate file
    :param Output: string for outputted files
    :param Method: Gradient Isotropic QHA ('GiQ');
                   Gradient Isotropic QHA w/ Gruneisen Parameter ('GiQg');
    :param Gradient_MaxTemp: Maximum temperature in gradient method
    :param Pressure: in atm
    :param LocGrd_Vol_FracStep:  isotropic volume fractional stepsize for local gradient
    :param LocGrd_Temp_step: temperature step size to use in local gradient 
    :param Statistical_mechanics: 'C' Classical mechanics
                                  'Q' Quantum mechanics
    :param NumAnalysis_step: stepsize for numerical method
    :param NumAnalysis_method: 'RK4' Runge-Kutta 4th order
    :param keyword_parameters: Parameter_file, Gruneisen_Vol_FracStep
    
    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    Gruneisen_Vol_FracStep: volume fraction step used to determine the Gruneisen parameter 
    """
    # Setting file endings and determining how many wavenumbers there will be
    if Program == 'tink':
        file_ending = '.xyz'
        number_of_wavenumbers = Pr.Tinker_atoms_per_molecule(Coordinate_file, 1)*3
    elif Program == 'test':
        file_ending = '.npy'
        number_of_wavenumbers = len(Wvn.Test_Wavenumber(Coordinate_file))
        keyword_parameters['Parameter_file'] = ''

    # Setting the temperature array
    temperature = np.arange(0, Gradient_MaxTemp + 1, NumAnalysis_step)

    # Setting the volume gradient array to be filled
    volume_gradient = np.zeros((len(temperature) - 1, 2))
    volume_gradient[:, 0] = temperature[:len(temperature)-1]
    if os.path.isfile(Output + '_' + Method + '_dV.npy'):
        volume_gradient_hold = np.load(Output + '_' + Method + '_dV.npy')
        # If the temperatures line up, then the previous local gradients will be used
        if len(volume_gradient_hold[:, 0]) <= len(volume_gradient[:, 0]):
            if all(volume_gradient_hold[:, 0] == volume_gradient[:len(volume_gradient_hold[:, 0]), 0]):
                volume_gradient[:len(volume_gradient_hold[:, 0]), :] = volume_gradient_hold

    # Setting up a matrix to store the wavenumbers in
    wavenumbers = np.zeros((len(temperature), number_of_wavenumbers+1))
    wavenumbers[:, 0] = temperature

    # Setting parameters for the Gruneisen parameter and loading in previously found wavenumbers for SiQ
    if Method == 'GiQg':
        Gruneisen, Wavenumber_Reference, Volume_Reference = Wvn.Call_Wavenumbers(Method,Coordinate_file=Coordinate_file,
                                                                                 Program=Program,
                                                                                 Gruneisen_Vol_FracStep=
                                                                                 keyword_parameters[
                                                                                     'Gruneisen_Vol_FracStep'],
                                                                                 molecules_in_coord=molecules_in_coord,
                                                                                 Parameter_file=
                                                                                 keyword_parameters['Parameter_file'])
    elif Method == 'GiQ':
        Gruneisen = 0.
        Wavenumber_Reference = 0.
        Volume_Reference = 0.
        if os.path.isfile(Output + '_' + Method + '_WVN.npy'):
            wavenumbers_hold = np.append(wavenumbers, np.load(Output + '_' + Method + '_WVN.npy'), axis=0)
            # If the temperatures line up in the previous wavenumber matrix, it will be used in the current run
            if len(wavenumbers_hold[:, 0]) <= len(wavenumbers[:, 0]):
                if all(wavenumbers_hold[:, 0] == wavenumbers[:len(wavenumbers_hold[:, 0]), 0]):
                    wavenumbers[:len(wavenumbers_hold[:, 0]), :] = wavenumbers_hold

    # Setting up an array to store the properties
    properties = np.zeros((len(temperature), 14))

    # Holding lattice structures as the structure at 0K
    os.system('cp ' + Coordinate_file + ' ' + Output + '_' + Method + 'T' + str(temperature[0]) + file_ending)

    # Finding structures at higher temperatures
    for i in range(len(temperature) - 1):
        if any(wavenumbers[i, 1:] != 0.) and (volume_gradient[i, 1] != 0.):
            pass
        else:
            if NumAnalysis_method == 'RK4':
                volume_gradient[i, 1], wavenumbers[i, 1:], volume = \
                    Runge_Kutta_Fourth_Order(Method, Output + '_' + Method + 'T' + str(temperature[i]) + file_ending,
                                             Program, temperature[i], Pressure, LocGrd_Temp_step,
                                             molecules_in_coord, Statistical_mechanics, NumAnalysis_step,
                                             Parameter_file=keyword_parameters['Parameter_file'],
                                             Gruneisen=Gruneisen, Wavenumber_Reference=Wavenumber_Reference,
                                             Volume_Reference=Volume_Reference, LocGrd_Vol_FracStep=LocGrd_Vol_FracStep)

        # Saving wavenumbers and local gradient information
        if Method == 'GiQ':
            np.save(Output + '_' + Method + '_WVN', wavenumbers[~np.all(wavenumbers == 0, axis=1)])
        np.save(Output + '_' + Method + '_dV', volume_gradient[~np.all(volume_gradient == 0, axis=1)])

        # Populating the properties for the current temperature
        properties[i, :] = Pr.Properties(Output + '_' + Method + 'T' + str(temperature[i]) + file_ending,
                                         wavenumbers[i, 1:], temperature[i], Pressure, Program, Statistical_mechanics,
                                         molecules_in_coord, Parameter_file=keyword_parameters['Parameter_file'])

        # Expanding to the next strucutre
        Ex.Call_Expansion(Method, 'expand', Program, Output + '_' + Method + 'T' + str(temperature[i]) + file_ending,
                          molecules_in_coord, Parameter_file=keyword_parameters['Parameter_file'],
                          volume_fraction_change=(volume + volume_gradient[i, 1]*NumAnalysis_step)/volume,
                          Output=Output + '_' + Method + 'T' + str(temperature[i + 1]))
        os.system('mv ' + Output + '_' + Method + 'T' + str(temperature[i]) + file_ending + ' Cords/')
        if temperature[i + 1] == temperature[-1]:
            volume = Pr.Volume(Program=Program,Coordinate_file=Output + '_' + Method + 'T' + str(temperature[i + 1]) +
                               file_ending)
            wavenumbers[i+1, 1:] = Wvn.Call_Wavenumbers(Method, Gruneisen=Gruneisen,
                                                  Wavenumber_Reference=Wavenumber_Reference,
                                                  Volume_Reference=Volume_Reference, New_Volume=volume,
                                                  Parameter_file=keyword_parameters['Parameter_file'],
                                                  Program=Program,
                                                  Coordinate_file=Output + '_' + Method + 'T' + str(temperature[i + 1])
                                                                  + file_ending)
            properties[i+1, :] = Pr.Properties(Output + '_' + Method + 'T' + str(temperature[i + 1]) + file_ending,
                                               wavenumbers[i + 1, 1:], temperature[i + 1], Pressure, Program,
                                               Statistical_mechanics, molecules_in_coord,
                                               Parameter_file=keyword_parameters['Parameter_file'])
            os.system('mv ' + Output + '_' + Method + 'T' + str(temperature[i + 1]) + file_ending + ' Cords/')

    # Saving the raw data before minimizing
    np.save(Output+'_raw', properties)
    return properties
        
##########################################
##### Gradient Anisotropic Expansion #####
##########################################
def Ansotropic_Gradient_Expansion(Coordinate_file, Program, molecules_in_coord, Output, Method, Gradient_MaxTemp,
                                  Pressure, LocGrd_LatParam_FracStep, LocGrd_Temp_step, Statistical_mechanics,
                                  NumAnalysis_step, NumAnalysis_method, Aniso_LocGrad_Type, **keyword_parameters):
    """
    This function calculated the anisotropic gradient for thermal expansion and returns the properties along that path
    :param Coordinate_file: file containing the lattice parameters (and coordinates)
    :param Program: 'tink' for Tinker Molecular Modeling
                    'test' for a test run
    :param molecules_in_coord: number of molecules in coordinate file
    :param Output: string for outputted files
    :param Method: Gradient Anisotropic QHA ('GaQ');
                   Gradient Anisotropic QHA w/ Gruneisen Parameter ('GaQg');
    :param Gradient_MaxTemp: Maximum temperature in gradient method
    :param Pressure: in atm
    :param LocGrd_LatParam_FracStep: anisotropic crystal matrix fractional stepsize for local gradient
    :param LocGrd_Temp_step: temperature step size to use in local gradient 
    :param Statistical_mechanics: 'C' Classical mechanics
                                  'Q' Quantum mechanics
    :param NumAnalysis_step: stepsize for numerical method
    :param NumAnalysis_method: 'RK4' Runge-Kutta 4th order
    :param Aniso_LocGrad_Type: 73 Hessians to calculate the complete anistropic gradient
                               25 for d**2G_dUdU only calculating the diagonals and off-diags. of the upper left 3x3 matrix
                               19 for d**2G_dUdU only calculating the uppder left 3x3 matrix
                               13 for d**2G_dUdU only calculating the diagonals
                               7  for d**2G_dUdU only calculating the upper left 3x3 matrix daigonals
    :param keyword_parameters: Parameter_file
    
    Optional Parameters
    Parameter_file: program specific file containing force field parameters
    """
    #Optional: Parameter_file,

    # Setting file endings and determining how many wavenumbers there will be
    if Program == 'tink':
        file_ending = '.xyz'
        number_of_wavenumbers = Pr.Tinker_atoms_per_molecule(Coordinate_file, 1)*3
    elif Program == 'test':
        file_ending = '.npy'
        number_of_wavenumbers = len(Wvn.Test_Wavenumber(Coordinate_file))
        keyword_parameters['Parameter_file'] = ''

    # Setting the temperature array
    temperature = np.arange(0, Gradient_MaxTemp + 1, NumAnalysis_step)

    # Setting the volume gradient array to be filled
    lattice_gradient = np.zeros((len(temperature) - 1, 3, 4))
    lattice_gradient[:, 0, 0] = temperature[:len(temperature)-1]
    if os.path.isfile(Output + '_' + Method + '_du.npy'):
        lattice_gradient_hold = np.load(Output + '_' + Method + '_du.npy')
        # If the temperatures line up, then the previous local gradients will be used
        if len(lattice_gradient_hold[:, 0, 0]) <= lattice_gradient[:, 0, 0]:
            if all(lattice_gradient_hold[:, 0, 0] == lattice_gradient[:len(lattice_gradient_hold[:, 0, 0]), 0, 0]):
                lattice_gradient[:len(lattice_gradient_hold[:, 0, 0]), :, :] = lattice_gradient_hold

    # Setting up a matrix to store the wavenumbers in
    wavenumbers = np.zeros((len(temperature), number_of_wavenumbers+1))
    wavenumbers[:, 0] = temperature

    # Setting parameters for the Gruneisen parameter and loading in previously found wavenumbers for SiQ
    if Method == 'GaQg':
        print "Anisotropic Gruneisen not yet made"
        sys.exit()
    elif Method == 'GiQ':
        if os.path.isfile(Output + '_' + Method + '_WVN.npy'):
            wavenumbers_hold = np.append(wavenumbers, np.load(Output + '_' + Method + '_WVN.npy'), axis=0)
            # If the temperatures line up in the previous wavenumber matrix, it will be used in the current run
            if len(wavenumbers_hold[:, 0]) <= wavenumbers[:, 0]:
                if all(wavenumbers_hold[:, 0] == wavenumbers[:len(wavenumbers_hold[:, 0]), 0]):
                    wavenumbers[:len(wavenumbers_hold[:, 0]), :] = wavenumbers_hold

    # Setting up an array to store the properties
    properties = np.zeros((len(temperature), 14))

    # Holding lattice structures as the structure at 0K
    os.system('cp ' + Coordinate_file + ' ' + Output + '_' + Method + 'T' + str(temperature[0]) + file_ending)

    # Finding structures at higher temperatures
    for i in range(len(temperature) - 1):
        if any(wavenumbers[i, 1:] != 0.) and (lattice_gradient[i, 0, 1] != 0.):
            pass
        else:
            if NumAnalysis_method == 'RK4':
                lattice_gradient[i, :, 1:], wavenumbers[i, 1:], ignore = \
                    Runge_Kutta_Fourth_Order(Method, Output + '_' + Method + 'T' + str(temperature[i]) + file_ending,
                                             Program, temperature[i], Pressure, LocGrd_Temp_step,
                                             molecules_in_coord, Statistical_mechanics, NumAnalysis_step,
                                             Parameter_file=keyword_parameters['Parameter_file'],
                                             LocGrd_LatParam_FracStep=LocGrd_LatParam_FracStep,
                                             Aniso_LocGrad_Type=Aniso_LocGrad_Type)

        # Saving wavenumbers and local gradient information
        if Method == 'GaQ':
            np.save(Output + '_' + Method + '_WVN', wavenumbers[~np.all(wavenumbers == 0, axis=1)])
        np.save(Output + '_' + Method + '_du', lattice_gradient[~np.all(lattice_gradient[:, :, 1:] == 0, axis=1)])

        # Populating the properties for the current temperature
        properties[i, :] = Pr.Properties(Output + '_' + Method + 'T' + str(temperature[i]) + file_ending,
                                         wavenumbers[i, 1:], temperature[i], Pressure, Program, Statistical_mechanics,
                                         molecules_in_coord, Parameter_file=keyword_parameters['Parameter_file'])

        # Expanding to the next strucutre
        Ex.Call_Expansion(Method, 'expand', Program, Output + '_' + Method + 'T' + str(temperature[i]) + file_ending,
                          molecules_in_coord, Parameter_file=keyword_parameters['Parameter_file'],
                          dcrystal_matrix=lattice_gradient[i, :, 1:]*NumAnalysis_step,
                          Output=Output + '_' + Method + 'T' + str(temperature[i + 1]))

        os.system('mv ' + Output + '_' + Method + 'T' + str(temperature[i]) + file_ending + ' Cords/')
        if temperature[i + 1] == temperature[-1]:
            wavenumbers[i+1, 1:] = Wvn.Call_Wavenumbers(Method, Parameter_file=keyword_parameters['Parameter_file'],
                                                        Program=Program,
                                                        Coordinate_file=Output + '_' + Method + 'T' +
                                                                        str(temperature[i + 1]) + file_ending)
            properties[i+1, :] = Pr.Properties(Output + '_' + Method + 'T' + str(temperature[i + 1]) + file_ending,
                                               wavenumbers[i + 1, 1:], temperature[i + 1], Pressure, Program,
                                               Statistical_mechanics, molecules_in_coord,
                                               Parameter_file=keyword_parameters['Parameter_file'])
            os.system('mv ' + Output + '_' + Method + 'T' + str(temperature[i + 1]) + file_ending + ' Cords/')

    # Saving the raw data before minimizing
    np.save(Output+'_raw', properties)
    return properties
