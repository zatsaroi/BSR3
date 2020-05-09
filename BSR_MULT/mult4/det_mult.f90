!======================================================================
      Subroutine DET_mult1
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

      End Subroutine DET_mult1


!======================================================================
      Subroutine DET_mult2
!======================================================================
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     coefficients between possible combinations of symmetries
!
!     Isym - new disribution (perutatiom) of orbitals
!     IPsym - pointer on last same orbutal in new list 
!----------------------------------------------------------------------

      USE mult_par
      USE conf_LS,       only: ne
      USE spin_orbitals

      Implicit none
      Integer :: i,i1, j,j1, k,k2, N1(2*ne), N2(2*ne)
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

      if(k.gt.1) Return
!----------------------------------------------------------------------

      Select case (k)

      Case(1)

       if(kpol.eq.0) Return

       Do i=1,NSYM; if(N1(i).le.0) Cycle; i1=i; Exit; End do
       Do i=1,NSYM; if(N2(i).le.0) Cycle; j1=i; Exit; End do

       Call ZNO_001(i1,j1)

      Case(0)

      if(kpol.eq.0) then
        Call ZNO_overlap
      else
        Call ZNO_000
      end if

      End Select

      End Subroutine DET_mult2



!====================================================================
      Subroutine ZNO_001(is,js)
!====================================================================
!    angular part of electric multipole operator between two det.w.f
!--------------------------------------------------------------------

      USE mult_par  
      Use spin_orbitals

      Implicit none
      Integer, Intent(in) :: is,js
      Integer :: i,j,i1,i2,j1,j2,k1,k2,idf,int,kz

      Integer, External :: Idet_fact, Incode_mult

      Call Radi_matr(Lsym(is),Msym(is),Ssym(is), &
                     Lsym(js),Msym(js),Ssym(js))

      if(ktype.eq.'E'.and.CA.eq.0.d0) Return
      if(ktype.eq.'M'.and.abs(CA)+abs(CB).eq.0.d0) Return

      i1 = 1; if(is.gt.1) i1=IPsym1(is-1)+1; i2=IPsym1(is)
      j1 = 1; if(js.gt.1) j1=IPsym2(js-1)+1; j2=IPsym2(js)
      Do i=i1,i2; k1=nnsym1(Isym1(i))
      Do j=j1,j2; k2=nnsym2(Isym2(j))

       idf = Idet_fact(i,0,j,0);  kz = (-1)**(kz1+kz2+i+j)

       Select case(ktype)
        Case('E')
         int = Incode_mult(1,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
        Case('M')
         int = Incode_mult(2,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
         int = Incode_mult(3,k1,k2);  Call Iadd_zoef(CB*kz,int,idf)
       End Select     

      End do; End do

      END Subroutine ZNO_001


!====================================================================
      Subroutine ZNO_000
!====================================================================
!    angular part of electric multipole operator between two det.w.f
!--------------------------------------------------------------------

      USE mult_par
      USE spin_orbitals

      Implicit none
      Integer :: i,j,i1,i2,k,k1,k2,is,idf,int,kz
      Integer, External :: Idet_fact, Incode_mult

      Do is = 1,NSYM

       Call Radi_matr(Lsym(is),Msym(is),Ssym(is), &
                      Lsym(is),Msym(is),Ssym(is))

       if(ktype.eq.'E'.and.CA.eq.0.d0) Cycle
       if(ktype.eq.'M'.and.abs(CA)+abs(CB).eq.0.d0) Cycle

       i1 = 1; if(is.gt.1) i1=IPsym1(is-1)+1; i2=IPsym1(is)
       Do i=i1,i2; k=Isym1(i); k1=nnsym1(k)
       Do j=i1,i2; k=Isym2(j); k2=nnsym2(k)

       idf = Idet_fact(i,0,j,0); kz = (-1)**(kz1+kz2+i+j)

       Select case(ktype)
        Case('E')
         int = Incode_mult(1,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
        Case('M')
         int = Incode_mult(2,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
         int = Incode_mult(3,k1,k2);  Call Iadd_zoef(CB*kz,int,idf)
       End Select     

       End do; End do

      End do

      END Subroutine ZNO_000


!====================================================================
      Subroutine ZNO_overlap
!====================================================================
!    computes overlap integral between two determinants
!    Calls: Idet_fact, Iadd_zoef, Incode_mult
!--------------------------------------------------------------------

      USE spin_orbitals

      Implicit none
      Integer(4) :: idf,int
      Real(8) :: C
      Integer, External :: Idet_fact, Incode_mult

      C = (-1)**(kz1+kz2)
      idf = Idet_fact (0,0,0,0)
      int = Incode_mult (0,1,1)

      Call Iadd_zoef (C,int,idf)

      End Subroutine ZNO_overlap


!====================================================================
      Subroutine Radi_matr(l1,m1,s1,l2,m2,s2)
!====================================================================
!     angular part of electric or magnetic transition operator 
!     between 'nlms' orbitals:
!
!            <n1,l1,m1,s1| T(kq) | n2,l2,m2,s2>
!--------------------------------------------------------------------

      Use mult_par 

      Implicit none

      Integer, Intent(in) :: l1,m1,s1,l2,m2,s2
      Real(8), External :: Z_3j, Z_3jj, Z_6jj, ZCLKL, CLEBCH

      CA = 0.d0; CB = 0.d0

      Select case(ktype)

      Case('E')
       CA = (-1)**(l1-m1) * Z_3jj(l1,-m1,kpol,qpol,l2,m2) 
       if(mltpol.eq.0) CA = CA * ZCLKL(l1,kpol,l2)
      Case('M')

       CA = (-1)**(l1-m1) * Z_3jj(l1,-m1,kpol,qpol,l2,m2) 
       if(mltpol.eq.0) CA = CA &
          * ZCLKL(l1,kpol-1,l2) * sqrt(1.d0*l2*(l2+1)*(l2+l2+1)) &
          * sqrt(kpol+kpol+1.d0) * (-1)**(l1+kpol+l2) &
          * Z_6jj(kpol-1,1,kpol,l2,l1,l2)

       CB = (-1)**(l1-m1) * Z_3jj(l1,-m1,kpol-1,qpol,l2,m2) * &
            (-1)**((2-s1)/2) * Z_3j(2,-s1+2,3,spol+spol+1,2,s2) * &
            CLEBCH(kpol-1,qpol,1,spol,kpol,mpol)

       if(mltpol.eq.0) CB = CB * ZCLKL(l1,kpol-1,l2) * sqrt(1.5d0)

      End Select

      if(s1.ne.s2) CA=0.d0

      End Subroutine Radi_matr

