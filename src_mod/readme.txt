Brief description of the code flow and the modification. 

1. Ground state calculation is in gndstate.f90. When dft+u is turned on, it will first generate the density matrix then the potential matrix, with the two lines: 
call gendmatmt
call genvmatmt

