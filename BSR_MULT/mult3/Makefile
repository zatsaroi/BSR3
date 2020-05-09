include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = mult3
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(SHELLS) $(ZCOM) 
 
S  =	mod_param.f90            mod_spinorb.f90          mod_termexp.f90          \
	mod_coef.f90             mod_zoef.f90             add_res.f90              \
	conf_loop.f90            def_ic.f90               det_mult.f90             \
	idet_fact.f90            mult3.f90                prepare.f90              \
	rw.f90                   r_conf.f90               term_loop.f90            \
	inf_mult.f90             
 
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

# -fbacktrace -fcheck=bounds 

 
