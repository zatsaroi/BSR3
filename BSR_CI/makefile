include ${HOME}/DEFS_03
 
#  make  PROGRAM with name:
 
N  = bsr_ci
 
# Libraries to be searched:
 
L  =  $(BS) $(ZCONFLS) $(ZCOM) $(LAPACK)

S  = MOD_bsr_ci.f90  MOD_c_data.f90    MOD_l_core.f90  MOD_z_core.f90  \
     bsr_ci.f90      read_arg.f90      read_conf.f90   read_bsw.f90    \
     matrix_LS.f90   diag.f90          matrix_rel.f90  matrix_J.f90    \
     inform.f90      read_int_bnk.f90  I_data.f90      g_factor.f90    \
     read_data.f90
 
O = $(S:.f90=.o)  
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
#	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) -fcheck=all -Wall $<

