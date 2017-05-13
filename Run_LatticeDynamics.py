#!/usr/bin/env python

import os
import numpy as np
import ThermodynamicProperties as Pr
import Wavenumbers as Wvn
import Thermal_NumericalAnalysis as TNA

def Lattice_Dynamics(Temperature=[0.0, 25.0, 50.0, 75.0, 100.0], Pressure=1, Method='HA', Program='Test',
                     Output='out', Coordinate_file='molecule.xyz', Parameter_file='keyfile.key',
                     molecules_in_coord=1, properties_to_save=['G', 'T'], NumAnalysis_method='RK4',
                     NumAnalysis_step=25.0,
                     LocGrd_Temp_step=0.01, LocGrd_Vol_FracStep=3e-02,
                     LocGrd_LatParam_FracStep=5e-05, StepWise_Vol_StepFrac=1.5e-3,
                     StepWise_Vol_LowerFrac=0.97, StepWise_Vol_UpperFrac=1.16,
                     Statistical_mechanics='Classical', Gruneisen_Vol_FracStep=1.5e-3,
                     Wavenum_Tol=-1, Gradient_MaxTemp=300.0, Aniso_LocGrad_Type='73', Gruneisen_order=1):

    if Method == 'HA':
        # Running the Harmonic Approximation
        if os.path.isfile(Output + '_' + Method + '_WVN.npy'):
            wavenumbers = np.load(Output + '_' + Method + '_WVN.npy')
        else:
            wavenumbers = Wvn.Call_Wavenumbers(Method, Program=Program, Coordinate_file=Coordinate_file,
                                               Parameter_file=Parameter_file)
            np.save(Output + '_' + Method + '_WVN', wavenumbers)
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
                                                       Gradient_MaxTemp, Pressure, LocGrd_LatParam_FracStep,
                                                       LocGrd_Temp_step, Statistical_mechanics, NumAnalysis_step,
                                                       NumAnalysis_method, Aniso_LocGrad_Type,
                                                       Parameter_file=Parameter_file)
        Pr.Save_Properties(properties, properties_to_save, Output, Method, Statistical_mechanics)

         
if __name__ == '__main__':
    
    import argparse 
    
    parser = argparse.ArgumentParser(description='Calculate free energies as a function of T using lattice dynamics')
    
    parser.add_argument('-a', '--temperature', dest='Temperature', nargs='+', type=float,
                        default=[0.0, 25.0, 50.0, 75.0, 100.0], help='Array of temperatures [K]')
    parser.add_argument('-b', '--pressure', dest='Pressure', default='1', help='Pressure (atm)')
    parser.add_argument('-c', '--method', dest='Method', default='HA', help='type of lattice dynamics to perform',
                        choices=['HA', 'SiQ', 'SiQg', 'GiQ', 'GiQg', 'GaQ', 'GaQg'])
    parser.add_argument('-d', '--output', dest='Output', default='out', help='Output string')
    parser.add_argument('-e', '--program', dest='Program', default='Tinker',
                        help='method to use to calculate wavenumbers', choices=['Tinker', 'Test'])
    parser.add_argument('-g', '--coordfile', dest='Coordinate_file', default='molecule.xyz', type=file,
                        help='coordinate file for crystal (in the format the program accepts)')
    parser.add_argument('-i', '--parameter file', dest='Parameter_file', default='keyfile.key', type=file,
                        help='parameter file for crystal (in the format the program accepts)')
    parser.add_argument('-j', '--molecules_in_coord', dest='molecules_in_coord', default=4, type=int,
                        help='number of molecules in the coordinate file')
    parser.add_argument('-k', '--properties', dest='properties_to_save', nargs='+', default=['G', 'T'],
                        help='calculated properties to store')
    parser.add_argument('-l', '--integration_method', dest='NumAnalysis_method', default='RK4', choices=['RK4'],
                        help='choice of method for temperature gradient integration')
    parser.add_argument('-m', '--integration_step', dest='NumAnalysis_step', default=25.0, type=float,
                        help='step size for temperature gradient integration')
    parser.add_argument('-n', '--gradient_tdelta', dest='LocGrd_Temp_step', default=0.01, type=float,
                        help='steps size for numerical calculation of gradient in T (in K)')
    #MRS: seems like volume and box lattice fraction parameters should be combined?
    parser.add_argument('-o', '--gradient_vdelta', dest='LocGrd_Vol_FracStep', default=3e-02, type=float,
                        help='step size fnumerical calculation of gradient in volume fraction (in fraction)')
    parser.add_argument('-p', '--gradient_boxdelta', dest='LocGrd_LatParam_FracStep', default=5e-05, type=float,
                        help='steps size for box lattice parameter fraction (in fraction)')
    parser.add_argument('-q', '--step_volfrac', dest='StepWise_Vol_StepFrac', default=1.5e-3, type=float,
                        help='step size for stepwise volume change')
    parser.add_argument('-r', '--step_volmin', dest='StepWise_Vol_LowerFrac', default=0.97, type=float,
                        help='minimum relative change for stepwise volume change')
    parser.add_argument('-s', '--step_volmax', dest='StepWise_Vol_UpperFrac', default=1.16, type=float,
                        help='maximum relative change for stepwise volume change' )
    parser.add_argument('-t', '--statistics', dest='Statistical_Mechanics', default='Classical',
                        choices=['Quantum', 'Classical'],
                        help='type of statistics to use for calculating the lattice dynamics')
    parser.add_argument('-u', '--gruneisen_vol_step', dest='Gruneisen_Vol_FracStep', type=float, default=1.5e-3,
                        help='Gruneisen volume fraction change')
    parser.add_argument('-v', '--wavenum_tol', dest='Wavenum_Tol', default=-1.0, type=float,
                        help='minimum negative wavenumber accepted (negative wavenumbers ignored)')
    parser.add_argument('-w', '--maxtemp', dest='Gradient_MaxTemp', default=300.0, type=float,
                        help='maximum temperature change for gradient approach')
    parser.add_argument('-x', '--anisotropy', dest='Aniso_LocGrad_Type', choices=['73'], default='73',
                        help='anisotropy approximation used')
    parser.add_argument('-y', '--order', dest='Gruneisen_order', choices=['First','Second'], default='First',
                        help='Order of Gruneisen parameter')

    args = parser.parse_args()

    Lattice_Dynamics(Temperature=args.Temperature,
                     Pressure=args.Pressure,
                     Method=args.Method,
                     Program=args.Program,
                     Output=args.Output,
                     Coordinate_file=args.Coordinate_file,
                     Parameter_file=args.Parameter_file,
                     molecules_in_coord=args.molecules_in_coord,
                     properties_to_save=args.properties_to_save,
                     NumAnalysis_method=args.NumAnalysis_method,
                     NumAnalysis_step=args.NumAnalysis_step,
                     LocGrd_Temp_step=args.LocGrd_Temp_step,
                     LocGrd_Vol_FracStep=args.LocGrd_Vol_FracStep,
                     LocGrd_LatParam_FracStep=args.LocGrd_LatParam_FracStep,
                     StepWise_Vol_StepFrac=args.StepWise_Vol_StepFrac,
                     StepWise_Vol_LowerFrac=args.StepWise_Vol_LowerFrac,
                     StepWise_Vol_UpperFrac=args.StepWise_Vol_UpperFrac,
                     Statistical_mechanics=args.Statistical_mechanics,
                     Gruneisen_Vol_FracStep=args.options.Gruneisen_Vol_FracStep,
                     Wavenum_Tol=args.Wavenum_Tol,
                     Gradient_MaxTemp=args.Gradient_MaxTemp,
                     Aniso_LocGrad_Type=args.Aniso_LocGrad_Type,
                     Gruneisen_order=args.Gruneisen_order)
