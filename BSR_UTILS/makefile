#     BSR  utils

include ${HOME}/DEFS_03

all: bound_bsw     bound_tab  bs_hwf         bsr_dd           bsw123         \
	bsw_tab       bsw_w      collect_rsol                                   \
	photo_collect photo_tab  print_bs       target_states                   \
	w_bsw         zf_bb_bsr  zf_cc_bsr    

# extract seperate states from "bound.nnn" file:
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

# calculation of dipole matrix elements between R-matrx states saved in rsol.nnn files
bsr_dd: bsr_dd.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

# merging  1.bsw + 2.bsw + 3.bsw + ... --> res.bsw
bsw123: bsw123.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# connverts the bsw-files with orbitals in B-spline repersentation 
# into tab-files suitable for grafic display
bsw_tab: bsw_tab.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# connvert the B-spline results in the MCHF format
bsw_w:  bsw_w.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# target, rsol.nnn  --> rsol.all
collect_rsol: collect_rsol.f90
	$(FC) -o $@ $(link) $< $(ZCOM)
	mv $@ $(BINDIR)

# target, threshoulds -> target.exp
convert_target: convert_target.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

# H.DAT  -->  f-values
f_values: f_values.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

# different photoionization paramters
photo_collect: photo_collect.f90
	$(FC) -o $@ $(link) $<  $(ZCOM) 
	mv $@ $(BINDIR)

# different photoionization paramters
photo_tab_new: photo_tab_new.f90
	$(FC) -o $@ $(link) $< -fbacktrace -fcheck=bounds $(ZCONFLS) $(ZCOM) 
	mv $@ $(BINDIR)

# prints B-splines in separate files:
print_bs: print_bs.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# H.DAT --> dipole matrix elements
s_values: s_values.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

# target states information
target_states: target_states.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
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

target_kpert: target_kpert.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

