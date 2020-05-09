include ${HOME}/DEFS_04

#  make  PROGRAM with name:
 
N  = bsr_mat4
M  = bsr_mat4_mpi
 
# Libraries to be searched:
 
L  =   $(BS) $(ZCONFLS) $(ZCOM)  
 
MO =	MOD_bsr_mat.f90          MOD_c_data.f90           l_core.f90               z_core.f90                                        

S  =	bsr_mat.f90              read_matrix.f90          record_matrix.f90        \
	state_res.f90            sub1.f90                 check_mat.f90            
	            
SS =	bs_orth.f90              check_det.f90            f_values.f90             \
	get_v.f90                gen_matrix.f90           get_orth_chan.f90        \
	idef_type.f90            jsym_int.f90                                      \
	I_data.f90               L_data.f90               O_data.f90               \
	read_arg.f90             read_bsw.f90             matrix_updates.f90       \
	read_conf.f90            read_data.f90            rw_mask.f90              \
	Z_data.f90               check_idef.f90           bsr_mat_inf.f90                       

SM = bsr_mat_mpi.f90          br_arg.f90               check_mat_mpi.f90        \
	collect_mpi.f90          read_matrix_mpi.f90      record_matrix_mpi.f90    \
	state_res_mpi.f90        stop_mpi.f90             sub1_mpi.f90                          
                      

OS = $(MO:.f90=.o) $(S:.f90=.o) $(SS:.f90=.o)

OM = $(MO:.f90=.o) $(SM:.f90=.o) $(SS:.f90=.o) 
 
$(N): $(OS)
	$(FC) -o $(N) $(OS) $(link) $(L) $(LAPACK)  
	mv $(N) $(BINDIR)
#	rm *.o *.mod 

$(M): $(OM)
	$(FC) -o $(M) $(OM) $(link) $(BSR_MPI_LIB) $(L) 
	mv $(M) $(BINDIR)
	rm *.o *.mod 

clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) -fbacktrace    -fcheck=bounds  $<


# -fbacktrace    -fcheck=bounds 
# -CB -traceback
#-g -fcheck=all -Wall -fbacktrace    -fcheck=bounds

