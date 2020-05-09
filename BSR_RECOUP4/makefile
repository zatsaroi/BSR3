# insert your compiler and flags:
include ${HOME}/DEFS_04
 
# make  PROGRAM with name:
 
N  = bsr_recoup4
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(ZCOM) $(SHELLS)
 
S  =	MOD_bsr_recoup.f90      add_block.f90      bsr_recoup.f90     \
	matrix_updates.f90      read_arg.f90       record_matrix.f90  \
	sub1.f90              
	
O = $(S:.f90=.o)  
 
$(N): $(O)
	$(FC) -o $(N)  $(O) $(link) $(L)
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
 
.f.o:
	$(FC) $(comp) $<
 
#  -fbounds-check