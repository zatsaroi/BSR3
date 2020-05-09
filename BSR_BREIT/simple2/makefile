FC = gfortran 
link = -O2 -I. 
comp = -c -O2 -I. 

#  make  PROGRAM with name:
 
N  = test_2conf
 
S  =	mod_boef_list.f90      coef_2conf.f90      test_coef_2conf.f90
	
O = det_exp.o $(S:.f90=.o)   
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) 
	rm *.o *.mod
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) $<
 
det_exp.o: det_exp.f
	$(FC) $(comp) det_exp.f

 
