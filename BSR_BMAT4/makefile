include ${HOME}/DEFS_04

#  make  PROGRAM with name:   

NS = bsr_bmat4
NP = bsr_bmat4_mpi
 
# Libraries to be searched:
 
L  = $(BS) $(ZCONFLS) $(SHELLS) $(ZCOM)  $(BS)
 
MO =	MOD_bsr_breit.f90        MOD_bsr_mat.f90                                    

SG =	breit_me.f90             prepare_det_expn.f90     read_arg.f90             \
	rw_is_need.f90           read_conf.f90            term_loop.f90            \
	check_idef.f90           do_coef.f90              idef_type.f90            \
	read_data.f90            open_br.f90                      

BM = bsr_bmat.f90             add_res.f90              conf_loop.f90
          
BR = br_arg.f90  

BP =	bsr_bmat_mpi.f90         add_res_mpi.f90          conf_loop_mpi.f90        \
	get_res_mpi.f90          send_res_mpi.f90         conf_calc.f90            \
	get_det_exp.f90          send_det_exp.f90         
        

OS = $(MO:.f90=.o) $(SG:.f90=.o)  $(BM:.f90=.o)

OP = $(MO:.f90=.o) $(SG:.f90=.o)  $(BR:.f90=.o)  $(BP:.f90=.o)

$(NS): $(OS)
	$(FC) -o $(NS) $(OS) $(link)  $(L)  $(LAPACK)
	mv $(NS) $(BINDIR)
	rm *.o *.mod 


$(NP): $(OP)
	$(FC) -o $(NP) $(OP) $(link) $(BSR_MPI) $(L) 
	mv $(NP) $(BINDIR)
	rm *.o *.mod 

clean:
	rm -f *.o *.mod 
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp)  $<
 
# -CB -traceback
