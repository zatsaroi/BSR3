include ${HOME}/DEFS_03 
 
N  = bsr_phot3 
 
# Libraries to be searched:
  
L  =  $(ZCOM)  $(LAPACK)

S  =	mod_phot.f90             ak_coef.f90              bsr_phot.f90             \
	phot_sec.f90             read_arg.f90             read_bound.f90   read_Ddat.f90     \
	read_Hdat.f90            read_Wdat.f90            sort_photo.f90           \
	sym_mat1.f90             zaface.f90               zopen.f90                
	
 
O = $(S:.f90=.o) crees.o
 
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
 

crees.o: crees.f
	$(FC) $(comp) crees.f
