#  test run for JK scattering calculations
#  for e + C problem (~ 120 min)

#------------------------------------------------------
#  Preparation  of target states:

#  Supposed we made independent Breit_Pauli CI calculations 
#  for even and odd states.  The results are given in files
#  even.c,even.j,even.w
#  odd.c,odd.j,odd.w, respectively.
#
#  The radial functions are also given
#  in ASCII frm-file, in order to be able to rum the test
#  on different computer platforms.
#  (see subfolder frm)

#  Extract the c-files and w-files for each J-state:

cfile even.j 1 0 2p2_3P0   0.00000001
cfile even.j 2 0 2p2_1S0   0.00000001
cfile even.j 1 2 2p2_3P1   0.00000001
cfile even.j 1 4 2p2_3P2   0.00000001
cfile even.j 2 4 2p2_1D2   0.00000001

cfile odd.j 1 0 2p3s_3P0   0.00000001
cfile odd.j 1 2 2p3s_3P1   0.00000001
cfile odd.j 2 2 2p3s_1P1   0.00000001
cfile odd.j 3 2 2s_2p3_3D1 0.00000001
cfile odd.j 1 4 2s_2p3_5S2 0.00000001
cfile odd.j 2 4 2p3s_3P2   0.00000001
cfile odd.j 3 4 2s_2p3_3D2 0.00000001
cfile odd.j 1 6 2s_2p3_3D3 0.00000001

#  get the B-spline representation for radial functions

w_bsw  2p2_3P0.w 
w_bsw  2p2_1S0.w 
w_bsw  2p2_3P1.w 
w_bsw  2p2_3P2.w 
w_bsw  2p2_1D2.w 

w_bsw  2p3s_3P0.w   
w_bsw  2p3s_3P1.w   
w_bsw  2p3s_1P1.w   
w_bsw  2s_2p3_3D1.w 
w_bsw  2s_2p3_5S2.w 
w_bsw  2p3s_3P2.w   
w_bsw  2s_2p3_3D2.w 
w_bsw  2s_2p3_3D3.w 

w_bsw  p_1o.w
w_bsw  p_3o.w
w_bsw  p_5o.w

#---------------------------------------------------------

#  test BSR JK-scattering calculations
#  accoding to scattering model given in
#  target 

bsr_prep3
bsr_conf3

bsr_breit3 klsp1=1 klsp2=6 oper=1111000
bsr_mat3 klsp1=1 klsp2=6
bsr_hd3 klsp1=1 klsp2=6

sum_hh 

#  Final result is the h.dat file.
#  To check if the test run is correct
#  compare the lowest R-matrix poles in 
#  the sum_hh.log and sum_hh.test files.

diff sum_hh.test sum_hh.log > diff

