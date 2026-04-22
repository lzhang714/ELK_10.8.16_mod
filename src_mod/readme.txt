
Brief description of the code flow and the modification: 

All added/modified lines are ended with "!LZ". 
So, for example, do a 
grep "!LZ" ./genvmatmt.f90
will show the modifications in it. 

1. 
Ground state calculation is carried out in gndstate.f90. 
When dft+u is turned on, it will first generate a density matrix then the potential matrix, with this two lines: 
call gendmatmt
call genvmatmt

2. 
subroutine gendmatmt (gendmatmt.f90) is not touched. 

3. 
subroutine genvmatmt (genvmatmt.f90) will generate the potential matrix, with this line: 
if (dftu /= 0) call vmatmtdu
This subroutine is modified by adding a few lines in the end, to make a copy of the potential matrix of the desired atom (when in the 1st scf iteration) into the newly defined quantity "vmatmt_const". 

4. 
subroutine vmatmtdu (vmatmtdu.f90) will do the real job to calculate the potential matrix. 
The most outer loop is over atoms. 
Within this loop, a few lines are added at the end to copy "vmatmt_const" back into the variable "vmatmt" to keep it constant. 

5. 
Necessary declaration and allocation of the new variables are added in moddftu.f90 and init0.f90.
