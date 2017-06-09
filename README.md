Codes created and maintained by: Nate Abraham and Michael Shirts; 
                                 Shirts Group; 
                                 Department of Chemical and Biological Engineering; 
                                 University of Colorado Boulder

Currently, all lattice dynamic codes are designed to run for a desired temperature range and a single pressure.

Supported Programs:
- Tinker Molecular Modeling: Our code assumes that the Tinker binaries are callable without a direct path given
- Test systems: We have provided two funcitons to compute the  lattice energy and wavenumbers. More information can
be found in Test_Systems/Test/0README.md.

Running the scripts:
1) Contained in this directory, is an input.py file containing all of the user specified options. Move a copy of this
input file to another directory containing the coordinate file (and parameter file).
2) Make necessary changes to the input file. Most of the settings have already been tuned for the best results, for a
large subset of molecules and shouldn't need to be changed.
3) To run the lattice dynamic calculation, in the command line of the working directory type:
Run_LatticeDynamics.py -i input.inp

Scripts:
-  	ThermodynamicProperties.py: Contains functions to calculate individual properties for lattice dynamics. 
  Cleaned up from previous version and seperates a lot of functions.
- Wavenumbers.py: Contains all wavenumber functions.
- Expand.py: Contains functions to expand unit cells or determine the local gradient of thermal expansion.
- input.py: General input file with flags for user inputs and runs the general lattice dynamic script at the end.
  'python input.py' will run the main code as long as Run_LatticeDynamics.py is callable from the current directory.
- Run_LatticeDynamics.py: main code to call subroutines and functions and output properties.

Old Scripts:
* These files are to be phased out and are used for testing new code
- Anisotropic_grad.py: Anistropic QHA with a RK4 gradient numerial analysis for thermal expansion. Only looks at 
      the diagonal of the 6x6 gradient matrix.
- H_aprox.py: Harmonic approximation
- Isotropic_grad.py: Isotropic QHA with RK4 gradient numerical analysis for thermal expansion.
- Isotropic_grad_gry.py: Isotropic QHA + Gruneisen parameter with RK4 gradient numerical analysis for thermal 
      expansion.
- NumHes.py: Call Tinker's numerical Hessian analysis script and outputs the numerical vibrational spectra.
      ** Currently, the edited code of Tinker (in echarge2.f) produces segfaults when using the testhess binary.
- properties.py: Contains functions to calculate individual properties for lattice dynamics.
- QHAsub.py: Contains general scripts for thermal expansion and determining local gradients of thermal expansion.
- QHA_iso_step.py: Isotropic QHA using a stepwise approach for thermal expansion.
- QHA_iso_step_gru.py: Isotropic QHA + Gruneisen parameter using a stepwise approach for thermal expansion.


