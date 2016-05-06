# insert your compiler and flags:
include ${HOME}/DEFS_03
 
# make  PROGRAM with name:
 
N  = bsr_prep3
 
# Libraries to be searched:
 
L  =  $(BS) $(ZCONFLS) $(ZCOM) $(LAPACK)

O  = bsr_prep.o 

bsr_prep: $(O)
	$(FC) -o $(N) $(link) $(O) $(L) 
	mv $(N) $(BINDIR)
	rm *.o 

bsr_prep.o: bsr_prep.f90
	$(FC) $(comp) bsr_prep.f90

