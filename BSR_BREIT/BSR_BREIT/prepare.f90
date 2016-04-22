!======================================================================
      Subroutine Pre_detexp 
!======================================================================
!
!     define the det. expansions and write the information
!     in scratch files 'nud' and 'nua'
!
!-----------------------------------------------------------------------

      USE configs;  USE inter;  USE term_exp; Use spin_orbitals 
      USE inout_br

      Implicit none 

      Integer(4) :: i,j, ij, k,m,kt,kdt,kti,it,it1,it2, jt
      Integer(4) :: IL_dif,IS_dif, mso,msoo,mss
      Integer(4), Allocatable, Dimension(:) :: IPT_jc, IP_kt
      Real(8), Allocatable :: C_det(:)
      Integer(4), Allocatable, Dimension(:,:) ::IP_det
      Real(8), Allocatable :: CC_det(:,:)
      Integer(4), External :: Iglq 

!----------------------------------------------------------------------
! ... IPT_jc - pointer on unconsidered configurations:

      if(Allocated(IPT_jc)) Deallocate(IPT_jc); Allocate(IPT_jc(nsymc))

! ... allocations:

      if(.not.allocated(Idet)) Allocate(Idet(ne),Jdet(ne))

!----------------------------------------------------------------------
! ... loop over conf. symmeteries:

      ic_case = 0;  rewind(nud); rewind(nua)

      Do ic=1,nsymc;  if(IC_need(ic).eq.0) Cycle

!----------------------------------------------------------------------
! ... define configuration ic:

       k=IC_term(ic); it1=mod(k,ibc); it2=k/ibc; kti=it2-it1+1
       it = IP_term(it1); no1=noccsh(it); k = 1
       Do i=1,no1
        ln1(i)=nocorb(i,it); iq1(i)=nelcsh(i,it)
        in(i)=k; k=k+iq1(i)
        md(i)=Iglq(ln1(i),iq1(i))
       End do

!----------------------------------------------------------------------
! ... define max. values for shell and intermidiate terms (LS1)
! ... and the number of different angular symmetries (kt):

       LS = 0; kt = 0
       if(Allocated(IP_kt)) Deallocate(IP_kt);  Allocate(IP_kt(kti))
       Do k =it1,it2
        it=IP_term(k); if(IT_stat(it).eq.0) Cycle
        kt = kt + 1; IP_kt(kt) = it
        Do i=1,no1
         if(LSL(i,it).gt.LS(i,2)) LS(i,2)=LSL(i,it)
         if(LSS(i,it).gt.LS(i,3)) LS(i,3)=LSS(i,it)
         if(LPL(i,it).gt.LS(i,4)) LS(i,4)=LPL(i,it)
         if(LPS(i,it).gt.LS(i,5)) LS(i,5)=LPS(i,it)
        End do
       End do

      if(kt.eq.0) then; IC_need(ic)=0; Cycle; end if

      ILT1 = ILT_ic(ic); IST1 = IST_ic(ic)

      if(Allocated(C_det)) Deallocate(C_det);  Allocate(C_det(kt))
!----------------------------------------------------------------------       
! ... cycle over all config. symmetries:
   
      IPT_jc = 1   
      Do jc = 1,nsymc

      i=max(ic,jc); j=min(ic,jc); ij=i*(i-1)/2+j

      if(JC_need(ij).eq.0) Cycle      
 
!----------------------------------------------------------------------
! ... define configuration jc (to define orthogonality):

      i=IC_term(jc); m=mod(i,ibc); jt=IP_term(m)
      no2=noccsh(jt)
      Do i=1,no2
       ln2(i)=nocorb(i,jt); iq2(i)=nelcsh(i,jt)
      End do
      ILT2 = ILT_ic(jc); IST2 = IST_ic(jc)

!----------------------------------------------------------------------
! ... orthogonality on l:

      if(ic.ne.jc) then
       nn1=iq1; nn2=iq2
       Do i=1,no1; Do j=1,no2
         if(nn2(j).eq.0) Cycle
         if(ln1(i).eq.ln2(j)) then
          k=min(nn1(i),nn2(j)); nn1(i)=nn1(i)-k; nn2(j)=nn2(j)-k
          if(nn1(i).eq.0) Exit
         end if
       End do; End do
       k = SUM(nn2(1:no2))
       if(k.gt.2) then
        Call DEF_IC(ic,jc); JC_need(ij)=0; IPT_jc(jc)=0
        Cycle  ! over jc
       end if
      end if

!----------------------------------------------------------------------
! ... angular orthogonality:

       IL_dif = iabs(ILT1-ILT2); IS_dif = iabs(IST1-IST2) 

       mso = ioper(4); msoo = ioper(5); mss = ioper(6)

       m = 1
       if(ic.ne.jc) then
        if(IL_dif.gt.4.or.IS_dif.gt.4) m = 0
        if(mss.eq.0.and.(IL_dif.gt.2.or.IS_dif.gt.2)) m = 0
        if(mss+mso+msoo.eq.0.and.(IL_dif.gt.0.or.IS_dif.gt.0)) m = 0
       end if

       if(m.eq.0) then
        Call DEF_IC(ic,jc); JC_need(ij)=0; IPT_jc(jc)=0
        Cycle  ! over jc
       end if

!----------------------------------------------------------------------
      if(IPT_jc(jc).eq.0) Cycle    ! here because we should check
                                   ! the above othogonalities
!----------------------------------------------------------------------
! ... define the det. expansion:

       MLT = min(ILT1,ILT2); MST = min(IST1,IST2)

       rewind(nua);  Call DET_EXPN (nua,kti,kt,kdt,IP_kt,C_det)

       if(kdt.eq.0) Stop 'Pre_detexp: kdt = 0'

!----------------------------------------------------------------------
! ... record results:

       ic_case=ic_case+1
       write(nud) ic,kt,kdt,ILT1,IST1,MLT,MST
       write(nud) IP_kt(1:kt)
       rewind(nua)      
       Allocate(CC_det(kt,kdt),IP_det(ne,kdt))
       Do i = 1,kdt;  read(nua) IP_det(:,i),CC_det(:,i); End do
       write(nud) IP_det
       write(nud) CC_det
       Deallocate(CC_det,IP_det)

! ... mark the configurations with the same term:

       Do i = jc,nsymc
        if(MLT.eq.min(ILT1,ILT_ic(i)).and. &
           MST.eq.min(IST1,IST_ic(i))) IPT_jc(i)=0
       End do

      End do    ! over jc
      End do    ! over ic

      if(Allocated(IP_kt)) Deallocate(IP_kt)
      if(Allocated(IPT_jc)) Deallocate(IPT_jc)
      if(Allocated(C_det)) Deallocate(C_det)

      End Subroutine Pre_detexp 




!======================================================================
      Subroutine Det_expn (nua,kti,kt,kdt,IP_kt,Cdet) 
!======================================================================
!
!     procedure of exaustion of all possible determinants
!
!----------------------------------------------------------------------

      USE configs;  USE inter; USE term_exp; Use spin_orbitals 

      Implicit none

      Integer(4) :: nua,kt,kti,kd,kdt, i,j,k,m, it, ii
      Integer(4), Dimension(kti) :: IP_kt
      Real(8), Dimension(kt) :: Cdet
      Real(8) :: C
      Real(8), External :: Clebsh,DETC_sh

      kd=0; kdt=0; i=1; nd(i)=1              
    1 kd = kd + 1
      ii = in(i)
      Call DET_sh(ln1(i),iq1(i),nd(i),ML(i),MS(i),Idet(ii)) 

      m = iabs(ML(i)); if(ML(i).lt.0) m = m + 2
      if(m.gt.LS(i,2)) go to 2
      m = iabs(MS(i)); if(MS(i).lt.0) m = m + 2
      if(m.gt.LS(i,3)) go to 2

      if(i.eq.1) then
       MLp(1)=ML(1); MSp(1)=MS(1)
      else
       MLp(i) = MLp(i-1)+ML(i)-1
       MSp(i) = MSp(i-1)+MS(i)-1

       m = iabs(MLp(i)); if(MLp(i).lt.0) m = m + 2
       if(m.gt.LS(i,4)) go to 2
       m = iabs(MSp(i)); if(MSp(i).lt.0) m = m + 2
       if(m.gt.LS(i,5)) go to 2
      end if

      if(i.lt.no1) then; i = i + 1; nd(i) = 1;  go to 1; end if

      if(MLp(no1).ne.MLt) go to 2
      if(MSp(no1).ne.MSt) go to 2


!--------------------------------------------------------------------
!                                            coefficient calculation:
      Cdet = 0.d0
      Do k=1,kt; it = IP_kt(k); C=1.d0
       Do j=1,no1
        ii = LSA(j,it); C=C*DETC_sh(ln1(j),iq1(j),ii,nd(j))
        if(C.eq.0.d0) Exit
       End do
       if(C.eq.0.d0) Cycle
       Do j = 2,no1
        C = C *  Clebsh(LPL(j-1,it),MLp(j-1), &
                          LSL(j,it),ML(j),LPL(j,it),MLp(j))
        if(C.eq.0.d0) Exit
        C = C *  Clebsh(LPS(j-1,it),MSp(j-1), &
                          LSS(j,it),MS(j),LPS(j,it),MSp(j))
        if(C.eq.0.d0) Exit
       End do
       Cdet(k) = C
      End do

      kdt=kdt+1; write(nua) Idet,Cdet

!--------------------------------------------------------------------
    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue


      End Subroutine Det_expn

