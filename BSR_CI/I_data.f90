!======================================================================
      Subroutine I_data(icase,kpol)
!======================================================================
!     processing of I-integrals from the module 'c_data'
!----------------------------------------------------------------------
!     we have only one structure:
!
! 1 1.0  Rk( . . . .)  ic, jc               -  bound-bound  
!
!     where .  denotes bound orbital
!
!     jcase = 2 -->  INT( p1 .;p2 . )
!----------------------------------------------------------------------
      Use c_data; Use bsr_ci

      Implicit none
      Character(1) :: sym_d, sym_r, sym_i, sym_j
      Integer :: icase, kpol, jcase
      Integer :: i,j, i1,i2, j1,j2, iii,jjj, ic,jc
      Real(8) :: C,S
      Real(8) :: xx(ns,ns), dd(ns,ns), v(ns), w(ns)
      Real(8), external :: Sum_AmB

      if(ncdata.eq.0) Return

      dd=0.d0; xx=0.d0; v=0.d0; w=0.d0

!----------------------------------------------------------------------
! ... prepare B_spline representation and define symmetries: 

      Select case(icase)
       Case( 3); Call MTK_cell(kpol); sym_i='n'; sym_j='n'    !  Tk
       Case( 4); Call MMK_cell(kpol); sym_i='s'; sym_j='s'    !  Mk
       Case( 5); Call MRK_cell(kpol); sym_i='s'; sym_j='s'    !  Rk
       Case( 8); Call MNK_cell(kpol); sym_i='s'; sym_j='s'    !  Nk
       Case( 9); Call MVK_cell(kpol); sym_i='n'; sym_j='s'    !  Vk
       Case(10); Call MNK_cell(kpol); sym_i='s'; sym_j='s'    !  Nk
       Case Default; Stop 'unknown case in I_data'
      End Select
      jcase=2; sym_d=sym_i; sym_r=sym_j

!----------------------------------------------------------------------
!                                                       Rk( . . ; . . )      
      iii=0; jjj=0; S = 0.d0
      Do j=1,ncdata;  i=IPT(j)

       if(K1(i).ne.iii) then
        i1 = K1(i)/ibi; i2 = mod(K1(i),ibi)
        Call Density(ns,ks,dd,pbs(1,i1),pbs(1,i2),sym_d)
        Call Convol(ns,ks,xx,dd,jcase,sym_i,sym_j)
       end if

       if(K1(i).ne.iii.or.K2(i).ne.jjj) then
        j1 = K2(i)/ibi; j2 = mod(K2(i),ibi)
        Call Density(ns,ks,dd,pbs(1,j1),pbs(1,j2),sym_r)
        S = SUM_AmB(ns,ks,xx,dd,sym_r)
       end if

       iii=K1(i);  jjj=K2(i)
       ic=k3(i);  jc=k4(i)

       C = S*cdata(i)
       if(icase.le.8) then
        HM(ic,jc) = HM(ic,jc) + C
       else
        write(nur) C, ic, jc, icase
       end if

      End do

      End Subroutine  I_data

