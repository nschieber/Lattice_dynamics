#!/usr/bin/env python

import os
import numpy as np
import Properties as Pr
import Wavenumbers as Wvn
import Thermal_NumericalAnalysis as TNA
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-a', dest = 'Temperature',default='0,25,50,75,100')
parser.add_option('-b', dest = 'Pressure',default='1')
parser.add_option('-c', dest = 'Method',default='HA')
parser.add_option('-d', dest = 'Output',default='out')
parser.add_option('-e', dest = 'Program',default='tink')

parser.add_option('-g', dest = 'Coordinate_file',default='molecule.xyz')
parser.add_option('-i', dest = 'Parameter_file',default='keyfile.key')
parser.add_option('-j', dest = 'molecules_in_coord',default='4')
parser.add_option('-k', dest = 'properties_to_save',default='G,T')
parser.add_option('-l', dest = 'NumAnalysis_method',default='RK4')
parser.add_option('-m', dest = 'NumAnalysis_step',default='25.')
parser.add_option('-n', dest = 'LocGrd_Temp_step',default='0.01')
parser.add_option('-o', dest = 'LocGrd_Vol_FracStep',default='3e-02')
parser.add_option('-p', dest = 'LocGrd_LatParam_FracStep',default='5e-05')
parser.add_option('-q', dest = 'StepWise_Vol_StepFrac',default='1.5e-3')
parser.add_option('-r', dest = 'StepWise_Vol_LowerFrac',default='0.97')
parser.add_option('-s', dest = 'StepWise_Vol_UpperFrac',default='1.16')
parser.add_option('-t', dest = 'Statistical_mechanics',default='C')
parser.add_option('-u', dest = 'Gruneisen_Vol_FracStep',default='1.5e-3')
parser.add_option('-v', dest = 'Wavenum_Tol',default='-1.0')
parser.add_option('-w', dest = 'Gradient_MaxTemp',default='300.')
parser.add_option('-x', dest = 'Aniso_LocGrad_Type',default='73') 

(options, args) = parser.parse_args()
# Temperature array [K]
Temperature = np.array(options.Temperature.split(',')).astype(float)
# Pressura [atm]
Pressure = float(options.Pressure)
# Lattice Dynamic method
Method = options.Method
# Output string
Output = options.Output
# Program to run with
Program = options.Program

# Initial coordinate file
Coordinate_file = options.Coordinate_file
# Parameter or input file
Parameter_file = options.Parameter_file
# Numeber of molecules in coordinate file
molecules_in_coord = float(options.molecules_in_coord)
# Properties to output
properties_to_save = np.array(options.properties_to_save.split(','))
# Numerical analysis
NumAnalysis_method = options.NumAnalysis_method
# Numerical analysis step size
NumAnalysis_step  = float(options.NumAnalysis_step)
# Local gradient temperature step size
LocGrd_Temp_step = float(options.LocGrd_Temp_step)
# Local gradient volume fraction stepsize
LocGrd_Vol_FracStep = float(options.LocGrd_Vol_FracStep)
# Local gradient lattice parameter fraction step size
LocGrd_LatParam_FracStep = float(options.LocGrd_LatParam_FracStep)
# Stepwise volume fraction stepsize
StepWise_Vol_StepFrac = float(options.StepWise_Vol_StepFrac)
# Stepwise lowerbound on volume fraction
StepWise_Vol_LowerFrac = float(options.StepWise_Vol_LowerFrac)
# Stepwise upperbound on volume fraction
StepWise_Vol_UpperFrac = float(options.StepWise_Vol_UpperFrac)
# Statistical Mechanics
Statistical_mechanics  = options.Statistical_mechanics
# Gruneisen volume fraction change
Gruneisen_Vol_FracStep = float(options.Gruneisen_Vol_FracStep)
# Gruneisen volume fraction change
Wavenum_Tol  = float(options.Wavenum_Tol)
# Maximum temperature
Gradient_MaxTemp  = float(options.Gradient_MaxTemp)
# Numerical shortcuts for anistropic gradient
Aniso_LocGrad_Type = int(options.Aniso_LocGrad_Type)


if Method == 'HA':
    # Running the Harmonic Approximation
    if os.path.isfile(Output + '_' + Method + '_WVN.npy'):
        wavenumbers = np.load(Output + '_' + Method + '_WVN.npy')
    else:
        wavenumbers = Wvn.Call_Wavenumbers(Method, Program=Program, Coordinate_file=Coordinate_file,
                                           Parameter_file=Parameter_file)
        #np.save(Output + '_' + Method + '_WVN', wavenumbers)
    if all(wavenumbers > Wavenum_Tol):
        properties = Pr.Properties_with_Temperature(Coordinate_file, wavenumbers, Temperature, Pressure, Program,
                                                    Statistical_mechanics, molecules_in_coord,
                                                    Parameter_file=Parameter_file)
        np.save(Output+'_raw', properties)
        Pr.Save_Properties(properties, properties_to_save, Output, Method, Statistical_mechanics)
else:
    if os.path.isdir('Cords') != True:
        os.system('mkdir Cords')

if (Method == 'SiQ') or (Method == 'SiQg'):
    # Stepwise Isotropic QHA
    properties = TNA.Isotropic_Stepwise_Expansion(StepWise_Vol_StepFrac, StepWise_Vol_LowerFrac, StepWise_Vol_UpperFrac,
                                                  Coordinate_file, Program, Temperature, Pressure, Output, Method,
                                                  molecules_in_coord, Wavenum_Tol, Statistical_mechanics,
                                                  Parameter_file=Parameter_file,
                                                  Gruneisen_Vol_FracStep=Gruneisen_Vol_FracStep)
    Pr.Save_Properties(properties, properties_to_save, Output, Method, Statistical_mechanics)

if (Method == 'GiQ') or (Method == 'GiQg'):
    # Gradient Isotropic QHA
    properties = TNA.Isotropic_Gradient_Expansion(Coordinate_file, Program, molecules_in_coord, Output, Method,
                                                  Gradient_MaxTemp, Pressure, LocGrd_Vol_FracStep, LocGrd_Temp_step,
                                                  Statistical_mechanics, NumAnalysis_step, NumAnalysis_method,
                                                  Parameter_file=Parameter_file,
                                                  Gruneisen_Vol_FracStep=Gruneisen_Vol_FracStep)
    Pr.Save_Properties(properties, properties_to_save, Output, Method, Statistical_mechanics)

if (Method == 'GaQ') or (Method == 'GaQg'):
    properties = TNA.Ansotropic_Gradient_Expansion(Coordinate_file, Program, molecules_in_coord, Output, Method,
                                                   Gradient_MaxTemp, Pressure, LocGrd_LatParam_FracStep, LocGrd_Temp_step,
                                                   Statistical_mechanics, NumAnalysis_step, NumAnalysis_method,
                                                   Aniso_LocGrad_Type, Parameter_file=Parameter_file)
    Pr.Save_Properties(properties, properties_to_save, Output, Method, Statistical_mechanics)
