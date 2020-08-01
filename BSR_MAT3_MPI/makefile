include ${HOME}/DEFS_03

#  make  PROGRAM with name:  
 
N  = bsr_mat3_mpi
 
# Libraries to be searched:
 
L  =   $(BS) $(ZCONFLS) $(ZCOM) # $(LAPACK)
 
S  =	spline_moments.f90       spline_integrals.f90     rk_moments.f90  mrk_cell.f90 \
	mod_param_mpi.f90        mod_dets.f90             mod_ndef.f90             \
	mod_ndet.f90             mod_cmdata.f90           mod_intmat_mpi.f90       \
	l_core.f90               z_core.f90                                        \
	add_cdata.f90            add_coef.f90             br_arg.f90               \
	br_bsorb.f90             br_channel.f90           br_conf_LS.f90           \
	br_core.f90              br_dets.f90              br_grid.f90              \
	br_orb_LS.f90            br_phys_orb.f90          br_symc_LS.f90           \
	br_symt_LS.f90           br_target.f90            bsr_mat_mpi.f90          \
	bs_orth_mpi.f90          check_det.f90            check_mat_mpi.f90        \
	det_fact.f90             f_values.f90             gen_matrix.f90           \
	get_V.f90                idef_type.f90                                     \
	inf_bsr_mat.f90          i_data.f90               jsym_int.f90             \
	l_data.f90               merge_cdata.f90                                   \
	o_data.f90               pri_int.f90              read_arg.f90             \
	read_bsw.f90             read_conf.f90            read_data.f90            \
	record_matrix_mpi.f90    state_res_mpi.f90        stop_mpi.f90             \
	sub1_mpi.f90             vdet.f90                 z_data.f90

O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC_mpi) -o $(N) $(O) $(link) $(L)  
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC_mpi) -traceback -CB $(comp) $<
 
