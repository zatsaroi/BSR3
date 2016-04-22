N      = coef_tab
 
OPTS0  = -fast -L$(MYLIBS)
OPTS1  = -c -fast 
 
LIBS   =  -lzconf90 -lzcom90
 
O  =   mod_param.o mod_int.o mod_cmdata.o mod_coeffs.o mod_zoef.o \
	mod_det.o mod_dets.o mod_detf.o \
	coef_tab.o idet_simple.o incode_int.o rr_conf.o


$(N): $(O)
	$(FC) -o $(N) $(OPTS0) $(O) $(LIBS)
	mv $(N) $(MYPROG)
	rm *.o *.mod

mod_param.o: mod_param.f90
	cp $(MYMODS)/*.mod .
	$(FC) $(OPTS1) mod_param.f90
 
mod_int.o: mod_int.f90
	$(FC) $(OPTS1) mod_int.f90
 
mod_cmdata.o: mod_cmdata.f90
	$(FC) $(OPTS1) mod_cmdata.f90
 
mod_coeffs.o: mod_coeffs.f90
	$(FC) $(OPTS1) mod_coeffs.f90
 
mod_zoef.o: mod_zoef.f90
	$(FC) $(OPTS1) mod_zoef.f90
 
mod_det.o: mod_det.f90
	$(FC) $(OPTS1) mod_det.f90
 
mod_dets.o: mod_dets.f90
	$(FC) $(OPTS1) mod_dets.f90
 
mod_detf.o: mod_detf.f90
	$(FC) $(OPTS1) mod_detf.f90
 
coef_tab.o: coef_tab.f90
	$(FC) $(OPTS1) coef_tab.f90
 
idet_simple.o: idet_simple.f90
	$(FC) $(OPTS1) idet_simple.f90
 
incode_int.o: incode_int.f90
	$(FC) $(OPTS1) incode_int.f90
 
rr_conf.o: rr_conf.f90
	$(FC) $(OPTS1) rr_conf.f90



