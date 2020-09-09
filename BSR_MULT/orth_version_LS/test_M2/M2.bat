cp ../test_mult_2conf.exe  .

test_mult_2conf 4p5_5s.c 4p5_5p.c M2 tab=coef_M2

rm mult_bnk*  
mult3 4p5_5s.c 4p5_5p.c M2
mult3_tab  c1=4p5_5s.c c2=4p5_5p.c bnk=mult_bnk  jort=-1  tab=coef_M2_mult

