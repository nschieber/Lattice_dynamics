#!/usr/bin/env python
import subprocess

##Input file

### General Inputs ###
## Temperature
# Import numpy array of temperatures
Temperature = '0,25,50,75,100,125,150,175,200,225,250,275,300'

## Pressure
# Only set up for single pressure values
Pressure = '1'

## Method
# Only one method is accepted currently
# Option - Description:
# 'HA'   - Harmonic approximation
# 'SiQ'  - Stepwise isotropic quasi-harmonic approximation
# 'SiQg' - Stepwise isotropic quasi-harmonic approximation with Gruneisen parameter
# 'GiQ'  - Gradient isotropic quasi-harmonic approximation
# 'GiQg' - Gradient isotropic quasi-harmonic approximation with Gruneisen parameter
# 'GaQ'  - Gradient anisotropic quasi-harmonic approximation
# 'GaQg' - Gradient anisotropic quasi-harmonic approximation with Gruneisen parameter
Method = 'GiQ'

## Output
# Output name for files
# Example: 'out' would output Gibbs energy as 'out_G.npy'
Output = 'out'

## Program
# Program for minimization and vibrational spectra calculation
# Option - Description:
# 'tink' - Tinker Molecular Modeling Package
# 'test' - Input functions in *.py to determine thermal expansion landscape and *.py to determine changes in wavenumber
Program = 'tink'

## Statistical mechanics
# Option - Description
# 'C' - Classical
# 'Q' - Quantum
Statistical_mechanics = 'C'

## Structure
# Coordinate file
Coordinate_file = 'molecule.xyz'

## Molecular parameters
# Input file for particular program
Parameter_file = 'keyfile.key'

## Nuber of molecules
#****I want to eventually replace this with a program specific code to look at connectivity
molecules_in_coord = '8'

## Properties
# Properties to output in individual files
# An <out>_raw.npy file will be output for all properties
#****Eventually make a program to take the <out>_raw.npy file and re-output other properties not taken from initial run
# Options - Description
# 'G' - Gibbs free energy [kcal/mol]
# 'S' - Entropy [kcal/(mol*T)]
# 'T' - Temperature [K]
# 'P' - Pressure [atm]
# 'Av' - Helmholtz vibrational energy [kcal/mol]
# 'V' - Volume [Ang.^3]
# 'u' - Lattice parameters [Ang., Ang., Ang., Deg., Deg., Deg.]
# 'U' - Potential energy [kcal/mol]
properties = 'G,u,V,T'


### Gradient options ###
## Numerical analysis for thermal expansion
# Option - Description
# 'none' - Uses the local gradient to predict next step
# 'RK4'  - Runge-kutta 4th order
NumAnalysis_method = 'RK4'

## Stepsize numerical analysis
# Right now only temperature is the only input option
NumAnalysis_step = '25.'

## Local gradient presets
# These options are tuned best for a wide array of options
# Temperature stepsize
LocGrd_Temp_step = '0.01'
# Isotropic volume fraction change
LocGrd_Vol_FracStep = '3e-02'
# Anisotropic lattice fraction change
LocGrd_LatParam_FracStep = '5e-02'

### Stepwise options ###
# Stepwise thermal expansion is only set up for isotropic expansion
# It is far too expansive to run anistropic expansion stepwise
# Volume fraction stepsize
StepWise_Vol_StepFrac = '1.5e-3'
# Volume fraction lowerbound
StepWise_Vol_LowerFrac = '0.99'
# Volume fraction upperbound
StepWise_Vol_UpperFrac = '1.02'

### Gruneisen options ###
# Volume fraction gruneisen change
Gruneisen_Vol_FracStep = '1.5e-3'
# 

### Wavenumber tolerance
Wavenum_Tol = '-1.0'

### Maximum temperature for gradient method
Gradient_MaxTemp = '100.'

### Number of Hessians for anistropic local gradient
# Option - Description
# '73'   - Compute every vibrational spectra for local gradient
# '25'   - For d**2G/dUdU compute vibrational spectra for the diagonal and the off diagonals for the upper left 3x3
# matrix
Aniso_LocGrad_Type = '73'

# Gruneisen Order
# Option - Description
# '1'    - gamma = dln(omega)/dln(V)
# '2'    - gamma = V/omega*[domega/dV + d^2omega/dV^2
Gruneisen_order = '1'


### Run Program ###
subprocess.call(['Run_LatticeDynamics.py',
                 '-a', Temperature,
                 '-b', Pressure,
                 '-c', Method,
                 '-d', Output,
                 '-e', Program,
                 '-g', Coordinate_file,
                 '-i', Parameter_file,
                 '-j', molecules_in_coord,
                 '-k', properties,
                 '-l', NumAnalysis_method,
                 '-m', NumAnalysis_step,
                 '-n', LocGrd_Temp_step,
                 '-o', LocGrd_Vol_FracStep,
                 '-p', LocGrd_LatParam_FracStep,
                 '-q', StepWise_Vol_StepFrac,
                 '-r', StepWise_Vol_LowerFrac,
                 '-s', StepWise_Vol_UpperFrac,
                 '-t', Statistical_mechanics,
                 '-u', Gruneisen_Vol_FracStep,
                 '-v', Wavenum_Tol,
                 '-w', Gradient_MaxTemp,
                 '-x', Aniso_LocGrad_Type])
