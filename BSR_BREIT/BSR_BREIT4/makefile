include ${HOME}/DEFS_04

#  make  PROGRAM with name:
 
N  = bsr_breit4
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(SHELLS) $(ZCOM) 
 
S  =	MOD_bsr_breit.f90                                              \
	BSR_BREIT.f90            breit_me.f90             add_res.f90              \
	conf_loop.f90            open_files.f90           read_arg.f90             \
	read_conf.f90            term_loop.f90            prepare_it_oper.f90    

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
	$(FC)  $(comp) $<

# gfortran:   -fbacktrace -fcheck=bounds 
