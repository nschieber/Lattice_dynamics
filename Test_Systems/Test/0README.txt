This system is to provide fast and efficient results for testing new code.

input.py: already modified for run just type 'python input.py' to run.
test.npy: contains an array of "lattice minimum" parameters for the test run.
    * Parameters are [a,b,c,alpha,beta,gamma] in 0-2 Angstroms and 3-5 in Degrees
    ** Lattice minimum parameters are those that minimize Test_U in Properties.py

If a different test is desired, two codes need to be edited:
    1) In Properties.py the Test_U function defines the potential energy of the system
       dependent on the lattice paraemters in the given coordinate file.
       * If the function is changed, make sure the input coordinate file is at a minimum
         of the function.
    2) In Wavenumbers.py the Test_Wavenumber outputs wavenumbers dependent on the coordinate
       file that is input. As of now the code has a standard set of wavenumbers for the
       lattice minimum and changes each dependent on the input coordinate file.

NSA: The wavenumbers need to be better designed here. I think you would have a better idea of how to do this, I want to keep
     it consistent with the way the Tinker wavenumbers are called. We could utalize the Parameter_file variablet to input
     wavenumbers or a function.

