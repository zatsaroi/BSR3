# BSR3: B-spline atomic R-matrix codes, third version

BSR is a general program to calculate atomic continuum processes using the B-spline R-matrix method, including 
electron-atom and electron-ion scattering, and radiative processes such as bound-bound transitions, photoionization and polarizabilities. The calculations can be performed in LS-coupling or in an intermediate-coupling scheme by including terms of the Breit-Pauli Hamiltonian. 

The present version is the deep recomposition of the original version published in

      >  Computer Physics Communications 174 (2006) 273–356 

Numerous new features and extansions are added, see doc folder in this repository and the references:

* Atomic structure calculations using MCHF and BSR
  >Oleg Zatsarinny and Charlotte Froese Fischer  
  >Computer Physics Communications 180 (2009) 2041–2065  

* The B-spline R-matrix method for atomic processes: 
  application to atomic structure, electron collisions and photoionization
  >Oleg Zatsarinny and Klaus Bartschat  
  >J. Phys. B: At. Mol. Opt. Phys. 46 (2013) 112001

# Installation

put DEF_03 (see Libraries_03) in your home directory, make corrections in the file regarding your compiler,
then compile the libraries first, then all programs in raw, using the Makefiles for each program
