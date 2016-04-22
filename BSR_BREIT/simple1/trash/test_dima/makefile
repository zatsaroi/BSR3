include ${HOME}/DEFS
 
#  make  PROGRAM with name:
 
N  = test_1conf
 
# Libraries to be searched:
 
L  =  $(ZCONFLS)  $(SHELLS)  $(ZCOM) 
 
S  =	mod_boef_list.f90        coef_1conf.f90           ee_coef.f90    \
	test_coef_1conf.f90
	
O = $(S:.f90=.o)  
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) $<
 
.f.o:
	$(FC) $(comp) $<
 
