!======================================================================
      Subroutine DET_breit2
!======================================================================
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     coefficients between possible combinations of symmetries
!
!     Isym - new disribution (perutatiom) of orbitals
!     IPsym - pointer on last same orbutal in new list 
!----------------------------------------------------------------------

      USE bsr_breit,     only: joper
      USE conf_LS,       only: ne
      USE spin_orbitals

      Implicit none
      Integer :: i,i1,i2, j,j1,j2, k,k2, N1(2*ne), N2(2*ne)
      Integer, External :: Isort

      NSYM = nsym1
!----------------------------------------------------------------------
! ... check for the same orbitals in the 2-nd configuration:      

      ksym2 = 1;  k2 = 0; IPsym2=0

      Do i = 1,NSYM
       Do j = 1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(i).ne.Lsym2(j)) Cycle
        if(Msym(i).ne.Msym2(j)) Cycle
        if(Ssym(i).ne.Ssym2(j)) Cycle
        k2=k2+1; IPsym2(i)=k2; Isym2(k2)=j; ksym2(j)=0
       End do        
      End do

! ... exzaust the 2-st configuration:

      Do i = 1,ne 
       if(ksym2(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym2(i); Msym(Nsym)=Msym2(i); Ssym(Nsym)=Ssym2(i)
       k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=i; ksym2(i)=0

! ... check for the same orbitals rest of 2-st configuration:      
       
       Do j = i+1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym2(j)) Cycle
        if(Msym(Nsym).ne.Msym2(j)) Cycle
        if(Ssym(Nsym).ne.Ssym2(j)) Cycle
        k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=j; ksym2(j)=0
       End do        
      End do

      if(k2.ne.ne) Stop 'Det_breit: k2 <> ne '
      IN=Isym2(1:ne);  kz2 = Isort(ne,IN)


      if(Nsym.gt.Nsym1) IPsym1(nsym1+1:nsym)=IPsym1(nsym1)
      Do i=2,nsym
       if(IPsym2(i).eq.0) IPsym2(i)=IPsym2(i-1)
      End do 
!----------------------------------------------------------------------
!                              define the number of different orbitals:
      Ksym1(1)=ipsym1(1)
      Ksym2(1)=ipsym2(1)
      Do i = 2,NSYM
       Ksym1(i)=ipsym1(i)-ipsym1(i-1)
       Ksym2(i)=ipsym2(i)-ipsym2(i-1)
      End do

! ... how much symmetries differ:

      k = 0
      Do i = 1,NSYM
       N1(i) = KSYM1(i)-KSYM2(i)
       N2(i) = KSYM2(i)-KSYM1(i)
       if(N1(i).gt.0) k = k + N1(i)
      End do

      if(k.gt.2) Return

!---------------------------------------------------------------------
!                                                         k = 2  case:
      Select case(k)

      Case(2)

       if(joper(3)+joper(5)+joper(6)+joper(7).eq.0) Return

       Do i=1,NSYM
        if(N1(i).le.0) Cycle; i1=i; N1(i)=N1(i)-1; Exit
       End do
       Do i=i1,NSYM
        if(N1(i).le.0) Cycle; i2=i; Exit
       End do
       Do i=1,NSYM
        if(N2(i).le.0) Cycle; j1=i; N2(i)=N2(i)-1; Exit
       End do
       Do i=j1,NSYM
        if(N2(i).le.0) Cycle; j2=i; Exit
       End do

       Call Zno_2ee(i1,i2,j1,j2)
  
!---------------------------------------------------------------------
!                                                         k = 1  case:
      Case(1)

       if(joper(3)+joper(5)+joper(6)+joper(7).eq.0) Return

       Do i=1,NSYM; if(N1(i).le.0) Cycle; i1 = i; Exit; End do
       Do i=1,NSYM; if(N2(i).le.0) Cycle; j1 = i; Exit; End do

       Do i = 1,Nsym
        if(Ksym1(i).eq.0) Cycle
        if(i.eq.i1.and.Ksym1(i).le.1) Cycle
        if(i.eq.j1.and.Ksym2(i).le.1) Cycle

        if(i.le.i1.and.i.le.j1)  then
          Call Zno_2ee(i,i1,i,j1)
        elseif(i.gt.i1.and.i.le.j1) then
          Call Zno_2ee(i1,i,i,j1)
        elseif(i.gt.i1.and.i.gt.j1) then
          Call Zno_2ee(i1,i,j1,i)
        elseif(i.le.i1.and.i.gt.j1) then
          Call Zno_2ee(i,i1,j1,i)
        end if

       End do

!---------------------------------------------------------------------
!                                                         k = 0  case:
      Case(0)

       if(joper(1).gt.0)           Call ZNO_0ee
       if(joper(2)+joper(4).gt.0)  Call ZNO_1ee

       if(joper(3)+joper(5)+joper(6)+joper(7).eq.0) Return

       Do i = 1,Nsym
        Do j = i,Nsym
         if(i.eq.j.and.Ksym1(i).le.1) Cycle
         Call Zno_2ee(i,j,i,j)
        End do
       End do

      End Select

      End Subroutine DET_breit2


!======================================================================
      Subroutine DET_breit1
!======================================================================
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     coefficients between possible combinations of symmetries
!
!     Isym - new disribution (perutatiom) of orbitals
!     IPsym - pointer on last same orbutal in new list 
!----------------------------------------------------------------------

	     USE spin_orbitals
	     USE conf_LS,         only: ne

      Implicit none
      Integer :: i, j, k1
      Integer, External :: Isort

      NSYM = 0

!----------------------------------------------------------------------
! ... exzaust the 1-st configuration:

      ksym1=1; k1=0

      Do i = 1,ne 
       if(ksym1(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym1(i); Msym(Nsym)=Msym1(i); Ssym(Nsym)=Ssym1(i)
       k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=i; ksym1(i)=0

! ... check for the same orbitals rest the 1-st configuration:      
 
       Do j = i+1,ne
        if(ksym1(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym1(j)) Cycle
        if(Msym(Nsym).ne.Msym1(j)) Cycle
        if(Ssym(Nsym).ne.Ssym1(j)) Cycle
        k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=j; ksym1(j)=0
       End do        

      End do

      if(k1.ne.ne) Stop 'Det_breit: k1 <> ne '
      IN=Isym1(1:ne);  kz1 = Isort(ne,IN)
      NSYM1 = NSYM

      End Subroutine DET_breit1