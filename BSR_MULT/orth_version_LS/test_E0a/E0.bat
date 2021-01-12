cp ../test_mult_2conf.exe  .

rm coef*

test_mult_2conf 1.c 2.c E0 tab=coef_E0

rm mult_bnk*  
mult3 1.c 2.c E0  mult_bnk  ovl=0
mult3_tab  c1=1.c c2=2.c bnk=mult_bnk  jort=1  tab=coef_E0_mult   ovl=0

