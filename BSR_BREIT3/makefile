include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = bsr_breit3
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(SHELLS) $(ZCOM) 
 
S  =	mod_param.f90            mod_termexp.f90          mod_spinorb.f90          \
	mod_boef.f90             mod_coef.f90             mod_zoef.f90             \
	add_it_oper.f90          add_res.f90              bsr_breit.f90            \
	conf_loop.f90            det_btreit.f90           idef_cme.f90             \
	idet_fact.f90            open_br.f90              prepare.f90              \
	read_arg.f90             r_conf.f90               term_loop.f90            \
	zno_0ee.f90              zno_1ee.f90              zno_2ee.f90              \
	zno_breit.f90

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
 
