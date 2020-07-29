#  test run for LS scattering calculations
#  for e + C problem (~ 5 min)

#  Preparation  of target states:

#  Each target states is represented by pair
#  of c- and w-file. It is supposed that
#  they were obtained in some separate MCHF
#  calculations. 

#  The radial functions are also given
#  in ASCII frm-file, in order to be able to rum the test
#  on different computer platforms.
#  (see subfolder frm)

#  get the B-spline represantation for
#  radial functions (bsw-files)

w_bsw 2p2_1D.w    
w_bsw 2p2_1S.w    
w_bsw 2p2_3P.w    
w_bsw 2p3s_1Po.w  
w_bsw 2p3s_3Po.w  
w_bsw 2p3_3Do.w   
w_bsw 2p3_5So.w   
w_bsw p_2Do.w     
w_bsw p_2Po.w     
w_bsw p_4So.w     

#  BSR LS-scattering calculations
#  accoding to scattering model given in
#  target 

bsr_prep3
bsr_conf3

bsr_breit3 klsp1=1 klsp2=12
bsr_mat3 klsp1=1 klsp2=12
bsr_hd3 klsp1=1 klsp2=12

sum_hh 1 12

#  Final result is the h.dat file.
#  To check if the test run is correct
#  compare the lowest R-matrix poles in 
#  the sum_hh.log and sum_hh.test files.

# fc sum_hh.log sum_hh.test > diff
diff sum_hh.log sum_hh.test > diff
