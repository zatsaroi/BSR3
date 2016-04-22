!======================================================================
      Subroutine ee_coef (l1,m1,ms1,l2,m2,ms2,l3,m3,ms3,l4,m4,ms4)
!======================================================================
!     computes angular coefficients for the two-electron interactio 
!     in uncoupled nlms-representation.
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1,m1,ms1,l2,m2,ms2,l3,m3,ms3,l4,m4,ms4
      Integer :: k, kl,km, kl1,kl2, km1,km2
      Real(8) :: S,C 
      Real(8), external :: Z_3jj

      if(ms1.ne.ms3.or.ms2.ne.ms4) Return

! ... define the range of multipole indeces and common multiplier:
    
      kl1=IABS(l1-l3);  kl2=IABS(l2-l4);  kl=MAX0(kl1,kl2)
      km1= l1+l3;       km2=l2+l4;        km=MIN0(km1,km2)
      if(km.lt.kl) Return

      S = (l1+l1+1)*(l2+l2+1)*(l3+l3+1)*(l4+l4+1)
      S = sqrt(S) * (-1)**(m1+m4)

! ... cycle on multipoles:
  
      Do k = kl,km,2                         
       C = S * Z_3jj(l1,0,l3,0,k,0);        if(C.eq.0.d0) Cycle
       C = C * Z_3jj(l1,-m1,l3,m3,k,m1-m3); if(C.eq.0.d0) Cycle
       C = C * Z_3jj(l2,0,l4,0,k,0);        if(C.eq.0.d0) Cycle
       C = C * Z_3jj(l2,-m2,l4,m4,k,m2-m4); if(C.eq.0.d0) Cycle
       Call Iadd_boef(C,k)
      End do

      End Subroutine ee_coef

