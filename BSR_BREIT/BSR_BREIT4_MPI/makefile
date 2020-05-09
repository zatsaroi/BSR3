include ${HOME}/DEFS_04

#  make  PROGRAM with name:   
 
N  = bsr_breit4_mpi
 
# Libraries to be searched:
 
L  =   $(ZCONFLS) $(SHELLS) $(ZCOM)
 
S  =	MOD_bsr_breit.f90        bsr_breit_mpi.f90                                 \
	add_it_oper.f90          add_res.f90              br_arg.f90               \
	breit_me.f90             conf_calc_mpi.f90        conf_loop_mpi.f90        \
	get_det_mpi.f90          get_res_mpi.f90                                   \
	open_files.f90           read_arg.f90                                      \
	read_conf.f90            send_det_mpi.f90         send_res_mpi.f90         \
	term_loop.f90            prepare_it_oper.f90

# mod_termexp.f90  mod_zoef.f90   

O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC_mpi) -o $(N) $(O) $(link) $(L) 
	mv $(N) $(BINDIR)
#	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC_mpi) $(comp)  $<
 
# -traceback -CB
