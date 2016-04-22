!======================================================================
      SUBROUTINE ZNO_2ee (i1,i2,j1,j2)
!======================================================================
!     angular part of matrix elements between two det.w.f.
!     for two-electron operator
!     Calls: Check_boef, Idet_fact, Iadd_zoef.
!----------------------------------------------------------------------

	     Use spin_orbitals, only: nsym, kz1,kz2,Lsym,Msym,Ssym, &
                               IPsym1,IPsym2,Isym1,Isym2,nnsym1,nnsym2

	     USE BOEF_list,     only: kblk,ncblk, ib_int, boef

      Implicit none
      Integer, Intent(in) :: i1,i2,j1,j2
      Integer :: i11,i12,i21,i22,j11,j12,j21,j22
      Integer :: io1,io2,io3,io4, ib1,ib2,ib3,ib4
      Integer :: i,j,k, k1,k2,k3,k4, met,int,idf,kz,ii1,ii2
      Integer :: ibint(4)
      Integer, External :: Idet_fact, Incode_int

!----------------------------------------------------------------------
      if(mod(Lsym(i1)+Lsym(i2)+Lsym(j1)+Lsym(j2),2).ne.0) Return
      if(Msym(i1)+Msym(i2).ne.Msym(j1)+Msym(j2)) Return
      if(Ssym(i1)+Ssym(i2).ne.Ssym(j1)+Ssym(j2)) Return
!----------------------------------------------------------------------

      Call Check_boef(Lsym(i1),Msym(i1),Ssym(i1), & 
                      Lsym(i2),Msym(i2),Ssym(i2), &
                      Lsym(j1),Msym(j1),Ssym(j1), &
                      Lsym(j2),Msym(j2),Ssym(j2))

!----------------------------------------------------------------------
!
      i11 = 1;  if(i1.gt.1) i11 = IPsym1(i1-1)+1;  i12 = IPsym1(i1)
      i21 = 1;  if(i2.gt.1) i21 = IPsym1(i2-1)+1;  i22 = IPsym1(i2)
      j11 = 1;  if(j1.gt.1) j11 = IPsym2(j1-1)+1;  j12 = IPsym2(j1)
      j21 = 1;  if(j2.gt.1) j21 = IPsym2(j2-1)+1;  j22 = IPsym2(j2)

      Do k1 = i11,i12;  i=Isym1(k1); ibint(1)=nnsym1(i)
      Do k2 = i21,i22;  j=Isym1(k2); ibint(2)=nnsym1(j)  
                         if(k2.le.k1) Cycle 

      Do k3 = j11,j12;  i=Isym2(k3);  ibint(3)=nnsym2(i)
      Do k4 = j21,j22;  j=Isym2(k4);  ibint(4)=nnsym2(j)
                         if(k4.le.k3) Cycle 

       idf = Idet_fact(k1,k2,k3,k4)

       kz = (-1)**(kz1+kz2+k1+k2+k3+k4)

       ii1 = 1; if(kblk.gt.1) ii1=ncblk(kblk-1)+1; ii2=ncblk(kblk)
       Do i = ii1,ii2
        Call Decode_int (met,k,ib1,ib2,ib3,ib4,IB_int(i))
        io1 = ibint(ib1); io2 = ibint(ib2) 
        io3 = ibint(ib3); io4 = ibint(ib4)
        Call Jsym_int(met,io1,io2,io3,io4)
        int = Incode_Int (met,k,io1,io2,io3,io4)
        Call Iadd_zoef(Boef(i)*kz,int,idf)

       End do

      End do;  End do;  End do;  End do

      END SUBROUTINE ZNO_2ee


!======================================================================
      Subroutine Jsym_int(icase,j1,j2,j3,j4)
!======================================================================

      Implicit none
      Integer :: icase,j,j1,j2,j3,j4

      if(icase.eq.3.or.icase.eq.4.or.icase.eq.5) then
       if(j1.gt.j2) then
        j=j1;j1=j2;j2=j; j=j3;j3=j4;j4=j
       elseif(j1.eq.j2.and.j3.gt.j4) then
        j=j3;j3=j4;j4=j
       end if
      end if

      End Subroutine Jsym_int