SHELL = /bin/sh

#  make  PROGRAM with name:
 
N  = bsr_hdb3
 
FC = mpif90
 
LIBDIR   = ${HOME}/T/LIBS_03
MODDIR   = ${HOME}/T/LIBS_03
BINDIR   = ${HOME}/BIN
 
link  =  -O2 -I. -I$(MODDIR)
comp  =  -c -O2 -I. -I$(MODDIR) -I${MKLROOT}/include
 
# Libraries to be searched:
 
L  =  -L$(LIBDIR) -lBS -lzconfLS -lzcom 

LIB = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential \
      -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl



S  =	mod_blacs.f90            mod_hd.f90               add_exp.f90              \
	add_mvc.f90              bsr_hdb.f90              b_out.f90                \
	diag_hd.f90              h_out.f90                pri_mainpar.f90          \
	read_diag.f90            rsol_out.f90             read_arg.f90             \
	sub1_hd.f90              transform_mat.f90        transform_ovl.f90        \
	w_out.f90                find_channel_label.f90   find_channel_nodes.f90



O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) -mkl $(link) $(L) $(LIB) 

	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) -traceback -CB $<
 
