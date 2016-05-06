include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = bsr_mat3
 
# Libraries to be searched:
 
L  =   $(BS) $(ZCONFLS) $(ZCOM) $(LAPACK)
 
S  =	mod_param.f90            mod_cmdata.f90           mod_dets.f90             \
	mod_intmat.f90           mod_ndef.f90             mod_ndet.f90             \
	l_core.f90               z_core.f90                                        \
	add_coef.f90             bsr_mat.f90              bs_orth.f90              \
	check_det.f90            det_fact.f90             f_values.f90             \
	gen_matrix.f90           get_V.f90                idef_type.f90            \
	ifind_channel.f90        I_data.f90               jsym_int.f90             \
	L_data.f90               merge_cdata.f90          O_data.f90               \
	pri_int.f90              read_arg.f90             read_bsw.f90             \
	read_conf.f90            read_data.f90            record_matrix.f90        \
	state_res.f90            sub1.f90                 vdet.f90                 \
	z_data.f90               check_mat.f90            add_cdata.f90            \
	bsr_mat_inf.f90
 
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
 
