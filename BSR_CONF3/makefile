# insert your compiler and flags:
include ${HOME}/DEFS_03
 
# make  PROGRAM with name:
 
N  = bsr_conf3
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(ZCOM) 
 
S  =	mod_bsr_conf.f90         bsr_conf.f90             check_cfg.f90            \
	check_pertuber.f90       def_orth_cond.f90        def_phys_targ.f90        \
	gen_conf.f90             make_coupling.f90        record_orth.f90          \
	sub_jj.f90               sub_jk.f90               sub_ls.f90               \
        read_arg.f90             bsr_conf_inf.f90
	
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
 
.f.o:
	$(FC) $(comp) $<
 
#  -fbounds-check