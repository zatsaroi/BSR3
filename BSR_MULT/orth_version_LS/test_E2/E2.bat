cp ../test_mult_2conf.exe  .

test_mult_2conf 1.c 2.c E2 tab=coef_E2  

rm mult_bnk*  
mult3 1.c 2.c E2
mult3_tab  c1=1.c c2=2.c bnk=mult_bnk  jort=-1  tab=coef_E2_mult

