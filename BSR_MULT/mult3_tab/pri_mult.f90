!======================================================================
      SUBROUTINE PRI_MULT(nu,C,itype,k,i1,i2,idf)
!======================================================================
!     Prints one angular integral
!----------------------------------------------------------------------

      Use orb_LS,     only: ELF
      Use mult_tab
      USE ndet_list
      Use ndef_list

      Implicit real(8) (A-H,O-Z)

      Character(1), Dimension(3) :: AI
      Data AI /'d','m','M'/
      Character(3) :: EL,EL1,EL2,EL3,EL4
      Character(220) :: A
      Character(1) :: a1,b1,a2

      if(abs(C).lt.Eps_C) Return

!----------------------------------------------------------------------

      Select case(itype)

      Case(0)                             ! Overlap
       ia=0
      Case(1,2,3)                         ! One-electron integrals
       write(A,'(a1,i1,a2,2(a4,a1))') &
             AI(itype),k,' (',ELF(i1),',',ELF(i2),')'
       ia=LEN_TRIM(A)
      Case default
       write(*,*) ' itype=',itype
       Stop ' Pri_coef: unknown case'
      End select

! ... determinant factor:

      if(idf.gt.0) then
       kd=KPF(idf);  ip=IPF(idf)
       Do i=1,kd
        ns=mod(NPF(ip+i),ibf); ni=NPF(ip+i)/ibf; 
        nd=KPD(ni); np=IPD(ni)
        ia=ia+1; A(ia:ia)='<'
       Do j=1,nd
        ni=NPD(np+j)/ibd
        write(A(ia+1:ia+5),'(a4,a1)') ELF(ni),','
        ia=ia+5
       End do
       A(ia:ia)='|'
       Do j=1,nd
        ni=mod(NPD(np+j),ibd)
        write(A(ia+1:ia+5),'(a4,a1)') ELF(ni),','
        ia=ia+5
       End do
       if(ns.gt.1) then
        A(ia:ia+1)='>^'; ia=ia+2; write(A(ia:ia),'(i1)') ns; ia=ia+1
       else
        A(ia:ia+1)='> '; ia=ia+2
       end if
      End do       ! Over kd
      end if       ! idf > 0

      Call Num(C,k1,k2,999999,1.d-9)

      a1='['; b1=':';a2=']'
      k1 = iabs(k1); m1=sqrt(1.*k1)+0.1;m2=sqrt(1.*k2)+0.1
      if(k1.eq.m1*m1.and.k2.eq.m2*m2) then
      k1=m1;k2=m2;a1='(';a2=')';end if
      write(nu,'(f12.6,1x,a1,i6,a1,i6,a1,2x,220a1)') &
                C,A1,iabs(k1),B1,k2,A2,(a(i:i),i=1,ia)

      END SUBROUTINE PRI_MULT



!====================================================================
      Subroutine NUM(z,k1,k2,i_max,acr)
!====================================================================
!
!    finds  k1 and k2 such that abs(z) = sqrt ( k1 / k2 )
!    and sign(z) = sign(k1)
!
!    k1, k2 < i_max
!    acr - tollerance for fitting
!--------------------------------------------------------------------

      Implicit none

      Real(8), intent(in) :: z, acr
      Integer(4), intent(in) :: i_max
      Integer(4), intent(out) :: k1,k2

      Integer(4) :: i_min, i1,i2, m1,m2 
      Real(8) :: zz, s,ss, s1,s2

      k1=0
      k2=1
      zz=z*z
      if(zz.lt.1.d0/i_max) Return

      s=zz
      i_min=zz
      if(i_min.lt.1) i_min=1
      Do i1=i_min,i_max
       s1=DBLE(i1)
       m1=s1/zz
       if(m1.lt.1) m1=1
       m2=m1+1
      Do i2=m1,m2
       s2=DBLE(i2)
       ss=abs(s1/s2-zz)
       if(ss.lt.s) then
        s = ss
        k1 = i1
        if(z.lt.0.d0) k1 = -k1
        k2 = i2
        end if
       if(ss.lt.acr) return
      End do
      End do

      End Subroutine NUM


