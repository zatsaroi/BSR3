include ${HOME}/DEFS_04

#  make  PROGRAM with name:
 
N  = mult4_mpi
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(SHELLS) $(ZCOM) 
 
S  =	mod_param.f90            mod_spinorb.f90          mod_termexp.f90          \
	mod_coef.f90             mod_zoef.f90             add_res.f90              \
	conf_loop.f90            def_ic.f90               det_mult.f90             \
	idet_fact.f90            mult.f90                 prepare.f90              \
	r_conf.f90               term_loop.f90                                     \
	stop_mpi.f90             read_arg.f90                                      \
	br_arg.f90               get_det_exp.f90          send_det_exp.f90         \
     conf_calc.f90            get_res_mpi.f90          send_res_mpi.f90
 
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

# -CB -traceback
 

 
