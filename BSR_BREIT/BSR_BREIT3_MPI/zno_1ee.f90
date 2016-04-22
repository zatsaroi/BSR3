!====================================================================
      Subroutine ZNO_1ee
!====================================================================
!    angular part of one-electron operator between two det.w.f
!    Calls: Idet_fact, Incode_int, Iadd_zoef.
!--------------------------------------------------------------------

      Use bsr_breit,     only: joper
	     Use spin_orbitals, only: nsym, kz1,kz2,Msym,Ssym,IPsym1, &
                               Isym1,Isym2,nnsym1,nnsym2

      Implicit none
      Integer :: i,j,i1,i2,k,k1,k2,is,idf,int
      Real(8) :: C,CFF
      Integer, External :: Idet_fact, Incode_int

      Do is = 1,NSYM
       i1 = 1; if(is.gt.1) i1=IPsym1(is-1)+1; i2=IPsym1(is)

       Do i=i1,i2; k=Isym1(i); k1=nnsym1(k)
       Do j=i1,i2; k=Isym2(j); k2=nnsym2(k)

       idf = Idet_fact(i,0,j,0)

       C=(-1)**(kz1+kz2+i+j)

       if(joper(2).gt.0) then                          ! L-integrals
        int = Incode_int (6,0,k1,k1,k2,k2)
        CFF=-0.5;  Call Iadd_zoef(C*CFF,int,idf)
       end if

       if(joper(4).gt.0.and.Msym(is).ne.0) then        ! Z-integrals
        C = C * Msym(is) * (Ssym(is)-1)
        int = Incode_int (7,0,k1,k1,k2,k2)
        Call Iadd_zoef(C,int,idf)
       end if

       End do
       End do
      End do

      END Subroutine ZNO_1ee
