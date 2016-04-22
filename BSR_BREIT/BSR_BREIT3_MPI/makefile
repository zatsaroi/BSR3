#  make  PROGRAM with name:   
 
N  = bsr_breit_mpi
 
FC = mpif90
 
LIBDIR   = ${HOME}/T/LIBS
MODDIR   = ${HOME}/T/LIBS
BINDIR   = ${HOME}/BIN
 
link  = -O2 -I. -I$(MODDIR) -L$(LIBDIR) 
comp  = -O2 -c -I. -I$(MODDIR) 


# -traceback -CB

# Libraries to be searched:
 
L  =   -lzconfLS -lshells -lzcom 
 
S  =	mod_param.f90            mod_spinorb.f90          mod_termexp.f90          \
	mod_boef.f90             mod_coef.f90             mod_zoef.f90             \
	add_it_oper.f90          add_res.f90              br_arg.f90               \
	bsr_breit_mpi.f90        conf_calc.f90            conf_loop.f90            \
	det_btreit.f90           get_det_exp.f90          get_res_mpi.f90          \
	idef_cme.f90             idet_fact.f90            incode_int.f90           \
	open_br.f90              prepare.f90              read_arg.f90             \
	r_conf.f90               send_det_exp.f90         send_res_mpi.f90         \
	term_loop.f90            zno_0ee.f90              zno_1ee.f90              \
	zno_2ee.f90              zno_breit.f90

O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L) $(LIBS) 
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
 
