include ${HOME}/DEFS
 
#  make  PROGRAM with name:
 
N  = coef_tab
 
# Libraries to be searched:
 
L  =  $(ZCONFLS) $(ZCOM) 
 
S  =    mod_param.f90   mod_cmdata.f90  mod_coeffs.f90  mod_det.f90  mod_detf.f90 \
	mod_dets.f90  mod_int.f90  mod_zoef.f90  coef_tab.f90  idet_simple.f90    \
	incode_int.f90 pre_iort.f90 rr_conf.f90 r_orth.f90
	
O = $(S:.f90=.o)  
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod $(N)$(exe)
 
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
 
