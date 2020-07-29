!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------
      Use bsr_breit,     only: pri,nui,nud, is_conf, js_conf, &
                               noper,joper,coper,JT_oper,CT_oper
      Use spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      Use term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, MLT,MST, ic_case
      Use conf_LS,       only: ne
      Use symc_list_LS,  only: JC_need, IC_need
      Use coef_list,     only: ntrm
      Use zoef_list,     only: nzoef

      Implicit none 
      Integer :: k1,k2, it,jt, ij, MLT2,MST2, is,js,ic,jc
      Real(8) :: t1,t2, C_ee,C_so,C_ss, zero=0.d0, one=1.d0  
      Integer, external :: IDEF_cme, DEF_ij
      Real(8), external :: Z_3j
      Character(80) :: conf

      Call Alloc_boef(-1)
      Call Alloc_blk (-1)

!----------------------------------------------------------------------
! ... cycle 1 over configurations:

      rewind(nud)
      Do is=1,ic_case;   is_conf = is

       Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT,MST

       if(Allocated(IP_kt1)) Deallocate(IP_kt1)
       Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

       if(Allocated(C_det1)) Deallocate(C_det1)
       Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

       if(Allocated(IM_det1)) Deallocate(IM_det1)
       Allocate(IM_det1(ne,kdt1)); Read(nud) IM_det1

       if(Allocated(IS_det1)) Deallocate(IS_det1)
       Allocate(IS_det1(ne,kdt1)); Read(nud) IS_det1

       read(nud) NNsym1(1:ne)
       read(nud) Lsym1(1:ne)

       if(IC_need(ic).eq.0) Cycle

       Call CPU_time(t1)

       Call Alloc_boef(-1)
       Call Alloc_blk (-1)

!----------------------------------------------------------------------
! ... cycle 2 over configurations:

      rewind(nud)
      Do js=1,is;     js_conf = js

       Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2

       if(Allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

       if(Allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

       if(Allocated(IM_det2)) Deallocate(IM_det2)
       Allocate(IM_det2(ne,kdt2)); Read(nud) IM_det2

       if(Allocated(IS_det2)) Deallocate(IS_det2)
       Allocate(IS_det2(ne,kdt2)); Read(nud) IS_det2

       read(nud) NNsym2(1:ne)
       read(nud) Lsym2(1:ne)

       if(MLT2.ne.MLT.or.MST2.ne.MST) Cycle
       if(MLT.ne.min(ILT1,ILT2).or.MST.ne.min(IST1,IST2)) Cycle
       ij=DEF_ij(ic,jc);  if(JC_need(ij).eq.0) Cycle      

!----------------------------------------------------------------------
! ...  define number of terms:

       ntrm = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(is.eq.js.and.it.gt.jt) Cycle;  ntrm = ntrm + 1
       End do; End do 

!----------------------------------------------------------------------
! ...  joper and JT_oper:

       if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
       Allocate(JT_oper(ntrm,noper),CT_oper(ntrm,noper))
       if(IDEF_cme(is,js).eq.0) Cycle 

!----------------------------------------------------------------------
! ... define the normalization constants for different operators:

       C_so = zero
       if(joper(4)+joper(5).ne.0) &
        C_so = Z_3j(ILT1,-MLT+2,3,1,ILT2,MLT)* &
               Z_3j(IST1,-MST+2,3,1,IST2,MST)* &
               (-1)**((ILT1-MLT+IST1-MST)/2)
       if(C_so.ne.zero) C_so=one/C_so

       C_ss = zero
       if(joper(6).ne.0) &
        C_ss = Z_3j(ILT1,-MLT+2,5,1,ILT2,MLT)* &
               Z_3j(IST1,-MST+2,5,1,IST2,MST)* &
               (-1)**((ILT1-MLT+IST1-MST)/2)
       if(C_ss.ne.zero) C_ss=one/C_ss

       C_ee = one; if(ILT1.ne.ILT2.or.IST1.ne.IST2) C_ee = zero

       if( abs(C_ee) + abs(C_so) + abs(C_ss) .eq. zero) Cycle

       coper(1) = zero; if(joper(1).gt.0) coper(1) = C_ee
       coper(2) = zero; if(joper(2).gt.0) coper(2) = C_ee
       coper(3) = zero; if(joper(3).gt.0) coper(3) = C_ee
       coper(4) = zero; if(joper(4).gt.0) coper(4) = C_so
       coper(5) = zero; if(joper(5).gt.0) coper(5) = C_so
       coper(6) = zero; if(joper(6).gt.0) coper(6) = C_ss
       coper(7) = zero; if(joper(7).gt.0) coper(7) = C_ee

!----------------------------------------------------------------------
! ...  initial allocations:

       Call Alloc_coef(-1)
       Call Alloc_ndet(-1)
       Call Alloc_ndef(-1)

!----------------------------------------------------------------------
! ...  calculations:

       Do kd1 = 1,kdt1

        Msym1(1:ne)=IM_det1(1:ne,kd1)
        Ssym1(1:ne)=IS_det1(1:ne,kd1)

        Call Det_breit1
 
       Do kd2 = 1,kdt2

        Msym2(1:ne)=IM_det2(1:ne,kd2)
        Ssym2(1:ne)=IS_det2(1:ne,kd2)

        nzoef = 0;      Call Det_breit2
        if(nzoef.gt.0)  Call Term_loop(is,js) 

       End do 
       End do

! ...  store results for given config.s:

       Call Add_res(nui,is,js); Call Add_it_oper(is,js)

      End do    ! over js

      Call CPU_time(t2)

      Call Symc_conf(ic,conf)
      write(*  ,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
        ' is=',is,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, &
          t2-t1,' sec',trim(conf)
      write(pri,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
        ' is=',is,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, &
          t2-t1,' sec',trim(conf)

      End do    ! over is

      End Subroutine Conf_loop


