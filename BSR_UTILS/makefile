#     BSR  utils

include ${HOME}/DEFS_03

all:  bound_bsw  bound_tab  bsw_tab    bsw_w      bsw123  \
      photo_tab  w_bsw      zf_cc_bsr  print_bs

# extract seprate states from "bound.nnn" file:
bound_bsw: bound_bsw.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCONFLS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# list of states in the "bound.nnn" files:
bound_tab: bound_tab.f90
	$(FC) -o $@ $(link) $< $(ZCOM)
	mv $@ $(BINDIR)

# hydrogen-like bound and pseudocontinuum states for given B-soline basis:
bs_hwf: bs_hwf.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# connverts the wbs-files with orbitals in B-spline repersentation 
# into tab-files suitable for grafic display
bsw_tab: bsw_tab.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# connvert the B-spline results in the MCHF format
bsw_w:  bsw_w.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# merging  1.bsw + 2.bsw + 3.bsw + ... --> res.bsw
bsw123: bsw123.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# different photoionization paramters
photo_tab: photo_tab.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM) 
	mv $@ $(BINDIR)

# prints B-splines in separate files:
print_bs: print_bs.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# Generate B-spline representation for orbitals 
# given in MCHF package format, w-files 
w_bsw: w_bsw.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

#  f-value calculations between states in bound.nnn files
zf_cc_bsr: zf_cc_bsr.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

#  f-value calculations between set of states in the c- and bsw-files
zf_bb_bsr: zf_bb_bsr.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

