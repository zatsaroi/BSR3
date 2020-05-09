!=======================================================================
      Integer Function Iort_conf()
!=======================================================================
! ... orthogonality on l between config.1 and config.2
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,j,k 

      Iort_conf = 0
      np1=iq1; np2=iq2
      Do i=1,no1
       Do j=1,no2
        if(np2(j).eq.0) Cycle
        if(ln1(i).ne.ln2(j)) Cycle
        k=min(np1(i),np2(j)); np1(i)=np1(i)-k; np2(j)=np2(j)-k
        if(np1(i).eq.0) Exit
       End do
      End do
      k = SUM(np2(1:no2));  if(k.gt.2) Iort_conf = 1

      End Function Iort_conf


!=======================================================================
      Subroutine Jort_conf(joper)
!=======================================================================
! ... orthogonality on l between config.1 and config.2 
! ... for different operators
!-----------------------------------------------------------------------

      Use conf_LS, L1=>Ltotal1, L2=>Ltotal2, S1=>Stotal1, S2=>Stotal2 

      Implicit none
      Integer :: joper(7)
	     Integer :: i,j,k     
      Integer, external :: Iort_conf, ITRI

      i=Iort_conf(); joper=i; if(i.eq.1) Return

      if(L1.ne.L2.or.S1.ne.S2) then
       joper(1:3)=1; joper(7)=1
      end if

      if(ITRI(L1,L2,3).eq.0) joper(4:5)=1
      if(ITRI(S1,S2,3).eq.0) joper(4:5)=1

      if(ITRI(L1,L2,5).eq.0) joper(6)=1
      if(ITRI(S1,S2,5).eq.0) joper(6)=1

      End Subroutine Jort_conf


!=======================================================================
      Integer Function Kort_conf()
!=======================================================================
! ... orthogonality on l between config.1 and config.2
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,j,k 

      np1=iq1; np2=iq2
      Do i=1,no1
       Do j=1,no2
        if(np2(j).eq.0) Cycle
        if(ln1(i).ne.ln2(j)) Cycle
        k=min(np1(i),np2(j)); np1(i)=np1(i)-k; np2(j)=np2(j)-k
        if(np1(i).eq.0) Exit
       End do
      End do
      Kort_conf = SUM(np2(1:no2))

      End Function Kort_conf


!=======================================================================
      Integer Function Iort_dipol(kpol,ktype,L1,S1,P1,L2,S2,P2)
!=======================================================================
! ... orthogonality for dipole transitions
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: kpol,L1,S1,P1,L2,S2,P2
      Character(1) :: ktype 
      Integer :: k,m  
      Integer, external :: ITRA, ITRI

      Iort_dipol = 1

       if(S1.ne.S2.and.ktype.eq.'E') Return

       if(P1.eq.P2) then
        if(ktype.eq.'E'.and.mod(kpol,2).ne.0) Return
        if(ktype.eq.'M'.and.mod(kpol,2).eq.0) Return
       else
        if(ktype.eq.'E'.and.mod(kpol,2).ne.1) Return
        if(ktype.eq.'M'.and.mod(kpol,2).eq.1) Return
       end if

       if(S1.eq.0) then  !  L = 2J
        k=kpol+kpol; 
        if(ITRA (L1,L2,k).eq.0) Return
       else
        k=kpol+kpol+1
        if(ktype.eq.'E'.and.ITRI (L1,L2,k).eq.0) Return
        m=kpol+kpol-1
        if(ktype.eq.'M'.and. &
         ITRI (L1,L2,k)+ITRI (L1,L2,m).eq.0) Return
      end if

      Iort_dipol = 0

      End Function Iort_dipol
