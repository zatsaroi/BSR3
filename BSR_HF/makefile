include ${HOME}/DEFS_03
 
#  make  PROGRAM with name:
 
N  = bsr_hf
 
# Libraries to be searched:
 
L  =  $(BS) $(ZCONFLS) $(SHELLS) $(ZCOM) $(LAPACK)

S  = MOD_bsr_hf.f90       MOD_hf_orbitals.f90    MOD_hf_energy.f90     MOD_boef_list.f90  \
     bsr_hf.f90           coef_ee_1conf.f90      I_int.f90             rk.f90             \
     get_case.f90         get_spline_param.f90   get_estimates.f90     check_tails.f90    \
     quadr_hf.f90         solve_hf.f90           boundary.f90                             \
     orthogonality.f90    bhl_hf.f90             rotate.f90            hf_matrix.f90      \
     update_hs.f90        write_cfg.f90          write_bsw.f90         summry.f90
 
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

