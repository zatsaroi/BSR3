include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = bsr_prep3
 
# Libraries to be searched:

L  =  $(BS) $(ZCONFLS) $(ZCOM) $(LAPACK)

S  = mod_bsr_prep.f90        bsr_prep.f90            bsr_prep_inf.f90            \
     check_bsw_file.f90      add_sub_orbital.f90     assign_index.f90            \
     check_bsr_name.f90      read_arg.f90                                        \
     sub_check_orb.f90       sub_phys_orb.f90        sub_phys_pert.f90     

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

#-fbounds-check  -fbacktrace
#-fbounds-check  -fbacktrace -fcheck=all -Wall 