include ${HOME}/DEFS_03
 
#  make  PROGRAM with name:
 
N  = mult3_tab
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(ZCOM) 
 
S  = mod_param.f90  Idet_simple.f90  mult_tab.f90  pri_mult.f90 r_conf.f90
	
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
 
