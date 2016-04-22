!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!
!     run loop over configurations 
!
!-----------------------------------------------------------------------

      USE param_br; USE configs; USE inter; USE spin_orbitals
      USE term_exp; USE coef_list; USE zoef_list; Use boef_list
      USE inout_br

      Implicit none 
      Integer(4) :: i,j,k,l,m,k1,k2,is,js, it,jt, ij, MLT2,MST2
      Integer(4), External :: ML_id, MS_id, Iglq,IDEF_cme
      Real(8) :: t1,t2,tt, C_ee,C_so,C_ss, zero=0.d0, one=1.d0  
      Real(8), External :: RRTC, Z_3j

!----------------------------------------------------------------------
!                                         define possible mls orbitals:
      m = 0
      Do is = 1,nsymt
       Do i=1,noccsh(is)
        l = nocorb(i,is); if(l.gt.m) m = l
       End do
      End do

      mls_max = 4*m + 2

      if(Allocated(ml_orb)) Deallocate(ml_orb,ms_orb)
      Allocate(ml_orb(mls_max), ms_orb(mls_max))

      Do i = 1,mls_max
       ml_orb(i)=ML_id(i)
       ms_orb(i)=MS_id(i)
      End do

!----------------------------------------------------------------------
!                                                          allocations:
      if(.not.Allocated(Lsym)) then
       m = 2*ne
       Allocate(Lsym(m),Lsym1(m),Lsym2(m),  &
                Msym(m),Msym1(m),Msym2(m),  &
                Ssym(m),Ssym1(m),Ssym2(m),  &
                IPsym1(m),IPsym2(m),Ksym1(m),Ksym2(m), &
                nnsym1(m),nnsym2(m),Isym1(m),Isym2(m)  ) 
       Allocate(N1(m),N2(m),N3(m),N4(m),NP(m))
      end if

!----------------------------------------------------------------------
!                                          cycle 1 over configurations:
      rewind(nud)
      Do is=1,ic_case

       Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT,MST

       if(Allocated(IP_kt1)) Deallocate(IP_kt1)
       Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

       if(Allocated(IP_det1)) Deallocate(IP_det1)
       Allocate(IP_det1(ne,kdt1)); Read(nud) IP_det1
             
       if(Allocated(C_det1)) Deallocate(C_det1)
       Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

       if(IC_need(ic).eq.0) Cycle

       it = IP_kt1(1); no1=noccsh(it); k=1
       Do i=1,no1
        ln1(i)=nocorb(i,it); iq1(i)=nelcsh(i,it)
        Do j=k,k+iq1(i)-1
         NNsym1(j)=i; Lsym1 (j)=ln1(i)
        End do
        k=k+iq1(i)
       End do

       Call Pri_conn(no1,ln1,iq1,CONFIG)
       write(*,'(a,4i8,a)') ' ic =',ic,nsymc,kt1,kdt1,CONFIG(1:8*no1)
       write(pri,'(a,4i8,a)') ' ic =',ic,nsymc,kt1,kdt1,CONFIG(1:8*no1)

       t1=RRTC()

!----------------------------------------------------------------------
!                                          cycle 2 over configurations:
      rewind(nud)
      Do js=1,is

       Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2

       if(Allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

       if(Allocated(IP_det2)) Deallocate(IP_det2)
       Allocate(IP_det2(ne,kdt2)); Read(nud) IP_det2
              
       if(Allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

       if(MLT2.ne.MLT.or.MST2.ne.MST) Cycle
!cd
       if(MLT.ne.min(ILT1,ILT2).or.MST.ne.min(IST1,IST2)) Cycle
       i = max(ic,jc); j = min(ic,jc); ij = i*(i-1)/2 + j
       if(JC_need(ij).eq.0) Cycle      
              
       jt=IP_kt2(1); no2=noccsh(jt); k=1
       Do i=1,no2
        ln2(i)=nocorb(i,jt); iq2(i)=nelcsh(i,jt)
        Do j=k,k+iq2(i)-1
         NNsym2(j)=i; Lsym2 (j)=ln2(i)
        End do
        k=k+iq2(i)
       End do

!       ILT2 = ILT_ic(jc); IST2 = IST_ic(jc) 

!----------------------------------------------------------------------
!                                               define number of terms:
       ntrm = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle;  ntrm = ntrm + 1
       End do; End do 
  
       if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
       Allocate(JT_oper(ntrm,noper),CT_oper(ntrm,noper))

!----------------------------------------------------------------------
       if(IDEF_cme().eq.0) then
        Call DEF_IC(ic,jc); Cycle 
       end if
!----------------------------------------------------------------------

       ncoef=0; Call Alloc_coef(ntrm,icoef)
       nboef=0; Call Alloc_boef(iboef)
       nblk=0;  Call Alloc_blk(iblk)

!----------------------------------------------------------------------
! ... define the normalization constants for different operators:

       C_so = zero
       if(joper(4)+joper(5).ne.0) &
        C_so = Z_3j(ILT1,-MLt+2,3,1,ILT2,MLt)* &
               Z_3j(IST1,-MSt+2,3,1,IST2,MSt)* &
               (-1)**((ILT1-MLt+IST1-MSt)/2)
       if(C_so.ne.zero) C_so=one/C_so

       C_ss = zero
       if(joper(6).ne.0) &
        C_ss = Z_3j(ILT1,-MLt+2,5,1,ILT2,MLt)* &
               Z_3j(IST1,-MSt+2,5,1,IST2,MSt)* &
               (-1)**((ILT1-MLt+IST1-MSt)/2)
       if(C_ss.ne.zero) C_ss=one/C_ss

       C_ee = one; if(ILT1.ne.ILT2.or.IST1.ne.IST2) C_ee = zero

       if((C_ee + C_so + C_ss) .eq. zero) Cycle

       coper(1) = zero; if(joper(1).gt.0) coper(1) = C_ee
       coper(2) = zero; if(joper(2).gt.0) coper(2) = C_ee
       coper(3) = zero; if(joper(3).gt.0) coper(3) = C_ee
       coper(4) = zero; if(joper(4).gt.0) coper(4) = C_so
       coper(5) = zero; if(joper(5).gt.0) coper(5) = C_so
       coper(6) = zero; if(joper(6).gt.0) coper(6) = C_ss
       coper(7) = zero; if(joper(7).gt.0) coper(7) = C_ee

!----------------------------------------------------------------------
! ...  calculations:

       Do kd1 = 1,kdt1

        Do i=1,ne
          m = IP_det1(i,kd1)
          Msym1(i)=(ML_orb(m)-1)/2; Ssym1(i)=MS_orb(m)
        End do

       Do kd2 = 1,kdt2

        Do i=1,ne
          m = IP_det2(i,kd2)
          Msym2(i)=(ML_orb(m)-1)/2; Ssym2(i)=MS_orb(m)
        End do

        Call Det_breit

        if(nzoef.gt.0) Call Term_loop 

       End do;  End do

!----------------------------------------------------------------------
! ...  store results for given config.s:

       if(ncoef.gt.0) Call Add_res(nui); Call DEF_IC(ic,jc)

      End do    ! over jc

      t2=RRTC();  tt=t2-t1
      if(tt.gt.100) write(pri,'(a,F12.2,a)') ' time: ',tt,' sec.'
      if(tt.gt.100) write(*,'(a,F12.2,a)') ' time: ',tt,' sec.'

      End do    ! over ic

      Do ic=1,nsymc;  Do jc=1,ic; Call DEF_IC(ic,jc); End do; End do 

      End Subroutine Conf_loop


!======================================================================
      Subroutine Pri_conn (n,ln,iq,CONFIG)
!======================================================================
!
!     gives the configuration without terms
!
!----------------------------------------------------------------------

      Character(*) :: CONFIG
      Character(4), External :: ELF4
      Integer(4), Intent(in) :: n
      Integer(4), Intent(in), Dimension(n) :: ln,iq
      Integer(4) :: k,i,ii 

      ii = LEN_TRIM(CONFIG)
      Do i=1,ii;  Config(i:i)=' ';  End do

      k=1
      Do i=1,n
       write(CONFIG(k:k+7),'(a4,''('',i2,'')'')') &
             ELF4(0,ln(i),0),iq(i)
       k=k+8
      End do

      End Subroutine Pri_conn

