Brief description of the code flow and the modification. 

1. Ground state calculation is in gndstate.f90. When dft+u is turned on, it will first generate a density matrix then the potential matrix, with the two lines: 
call gendmatmt
call genvmatmt

2. gendmatmt (gendmatmt.f90) is not touched. 

3. genvmatmt (genvmatmt.f90) will generate the potential matrix: 
if (dftu /= 0) call vmatmtdu
This is modified by adding a few lines in the end, to make a copy of the potential matrix of the desired atom in "vmatmt_const". 

4. 
