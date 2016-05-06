include ${HOME}/DEFS_03 
 
#  make  PROGRAM with name:
 
N  = bsr_hd3
 
# Libraries to be searched:
 
L  =  $(BS)  $(ZCONFLS) $(ZCOM) $(LAPACK) 

S  =	mod_hd.f90               bsr_hd.f90              b_out.f90                \
	diag_hd.f90              h_out.f90               rsol_out.f90             \
	r_arg.f90                sub1_hd.f90             w_out.f90                \
	find_channel_label.f90   bsr_hd_inf.f90          

 
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
 
