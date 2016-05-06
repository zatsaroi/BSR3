include ${HOME}/DEFS_03
 
N  = bsr_pol3
 
# Libraries to be searched:
 
L  =  $(BS) $(ZCONFLS) $(ZCOM) $(LAPACK)

S  =	mod_pol.f90              bsr_pol.f90              read_bsw.f90            \
	read_data.f90            r_arg.f90                r_bsrmat.f90            \
	r_dipmat.f90             solv_mat.f90             r_nortb.f90             \
	inf_bsr_pol.f90

 
O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod *~
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) $<
 
