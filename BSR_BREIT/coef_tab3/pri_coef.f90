!====================================================================
      SUBROUTINE PRI_COEF (nu,ic1,ic2,mo,mi,mso,mee,msoo,mss,moo)
!====================================================================
!     Prints the resulting angular integrals
!--------------------------------------------------------------------
      USE rk4_data; Use int_list

! ... print HEADER:

      write(nu,'(/70(''=''))')
      write(nu,'(12x,''< state'',i2,'' || (H-E) || state'',i2,''>'')') &
            ic1,ic2
      write(nu,'(70(''=''))') 
      Call Pri_conf (nu,ic1,0.d0);  write(nu,'(70(''-''))')
      Call Pri_conf (nu,ic2,0.d0);  write(nu,'(70(''-''))')

      if(nrk.eq.0) Return

!----------------------------------------------------------------------
!                                                             Overlaps:
      if(mo.eq.0) then 
       write(nu,'(/a/)') 'Overlaps:  no calculations!'
      elseif(mo.gt.0) then
       write(nu,'(/12x,'' Overlaps''/)')
       Do i=1,nrk; Call Pri_coef1(nu,i,11,ic1,ic2); End do
      end if
!---------------------------------------------------------------------
!                                                 Coulomb interaction:
      if(mi+mee.eq.0) then 
       write(nu,'(/a/)') 'Coulomb interactions:  no calculations!'
      elseif(mi+mee.gt.0) then
       write(nu,'(/12x,'' Coulomb interaction''/)')
       Do i=1,nrk; Call Pri_coef1(nu,i,5,ic1,ic2); End do
       Do i=1,nrk; Call Pri_coef1(nu,i,6,ic1,ic2); End do
      end if
!----------------------------------------------------------------------
!                                              Orbit-orbit interaction:
      if(moo.eq.0) then 
       write(nu,'(/a/)') ' Orbit-orbit:          no calculations!'
      elseif(moo.gt.0) then
       write(nu,'(/12x,''  Orbit-orbit interaction ''/)')
       Do i=1,nrk; Call Pri_coef1(nu,i,3,ic1,ic2); End do
       Do i=1,nrk; Call Pri_coef1(nu,i,4,ic1,ic2); End do
      end if
!----------------------------------------------------------------------
!                                               Spin-orbit interaction:
      if(mso.eq.0) then 
       write(nu,'(/a/)') ' Spin-orbit:           no calculations!'
      elseif(mso.gt.0) then
       write(nu,'(/12x,''  Spin-orbit interaction ''/)')
       Do i=1,nrk; Call Pri_coef1(nu,i,7,ic1,ic2); End do
      end if
!----------------------------------------------------------------------
!                                         Spin-other-orbit interaction:
      if(msoo.eq.0) then 
       write(nu,'(/a/)') ' Spin-other-orbit:     no calculations!'
      elseif(msoo.gt.0) then
       write(nu,'(/12x,''  Spin-other-orbit interaction ''/)')
       Do i=1,nrk; Call Pri_coef1(nu,i,8,ic1,ic2); End do
       Do i=1,nrk; Call Pri_coef1(nu,i,9,ic1,ic2); End do
      end if
!---------------------------------------------------------------------
!                                               Spin-spin interaction:
      if(mso.eq.0) then 
       write(nu,'(/a/)') ' Spin-spin:            no calculations!'
      elseif(mso.gt.0) then
       write(nu,'(/12x,''  Spin-spin interaction ''/)')
       Do i=1,nrk; Call Pri_coef1(nu,i,10,ic1,ic2); End do
      end if
!----------------------------------------------------------------------
      write(nu,'(70(''-'')/)')

      End SUBROUTINE PRI_COEF
 

!======================================================================
      SUBROUTINE PRI_COEF1 (nu,icoef,met,ic1,ic2)
!======================================================================
!     Prints one angular integral
!----------------------------------------------------------------------
      USE param_br
      USE conf_LS
      USE orb_LS
      USE rk4_data, iz_int => kr1, iz_df => kr2
      USE int_list
      USE ndet_list
      USE ndef_list

      Implicit none

      Character(1), Dimension(10) :: AI
      Data AI /' ',' ','T','M','R','L','Z','N','V','N'/
      Character(3) :: EL,EL1,EL2,EL3,EL4
      Character(220) :: A
      Character(1) :: a1,b1,a2
      Integer, Intent(in) :: nu,icoef,met,ic1,ic2
      Integer :: i,j, i1,i2,i3,i4, k,k1,k2, m1,m2
      Integer :: ia,ip,jp,idf,int,kd,nd,ns, ii

      if(abs(crk(icoef)).lt.Eps_C) Return

      if(ic1.ne.kr3(icoef).or.ic2.ne.kr4(icoef)) Return
!----------------------------------------------------------------------

      i=IZ_int(icoef); int=Jcase(i); if(int.ne.met) Return
      k = Jpol(i)
      if(int.eq.4.or.int.eq.8.or.int.eq.9.or.int.eq.10) k=k-1
      i1 = J1(i); i2 = J2(i); i3 = J3(i); i4 = J4(i)

      A = ' '

      Select case(int)

      Case(11)                            ! Overlap
      ia=0

      Case(6,7)                           ! One-electron integrals
      EL1 = ELF(I1)(2:4); EL2 = ELF(I3)(2:4)
      Write(A(1:15),'(a1,a4,2(a3,a1))') AI(int),'   (',EL1,',',EL2,')'
      ia=15

      Case(3,4,5,8,9,10)                   ! Two-electron integrals
      EL1 = ELF(I1)(2:4); EL2 = ELF(I2)(2:4);
      EL3 = ELF(I3)(2:4); EL4 = ELF(I4)(2:4)
      Write(A(1:23),'(a1,i2,a2,4(a3,a1))') &
        AI(int),k,' (',EL1,',',EL2,';',EL3,',',EL4,')'
      ia=23
 
      Case default
        write(*,*) ' int=',int
        Stop ' Pri_coef: unknown case'
      End select

! ... determinant factor:

      idf=IZ_df(icoef)
      if(idf.gt.0) then
       kd=KPF(idf);  ip=IPF(idf)
       Do i=1,kd
        ns=mod(NPF(ip+i),ibf); ii=NPF(ip+i)/ibf; 
        nd=KPD(ii); jp=IPD(ii)
        ia=ia+1; A(ia:ia)='<'
       Do j=1,nd
        ii=NPD(jp+j)/ibd;  EL = ELF(ii)(2:4)
        Write(A(ia+1:ia+5),'(10(a3,'',''))') EL
        ia=ia+4
       End do
       A(ia:ia)='|'
       Do j=1,nd
        ii=mod(NPD(jp+j),ibd); EL = ELF(ii)(2:4)
        Write(A(ia+1:ia+5),'(10(a3,'',''))') EL
        ia=ia+4
       End do
       if(ns.gt.1) then
        A(ia:ia+1)='>^'; ia=ia+2; write(A(ia:ia),'(i1)') ns; ia=ia+1
       else
        A(ia:ia+1)='> '; ia=ia+2
       end if
      End do       ! Over kd
      end if

      Call Num(crk(icoef),k1,k2,999999,1.d-9)

      a1='['; b1=':';a2=']'
      k1 = iabs(k1); m1=sqrt(1.*k1)+0.1;m2=sqrt(1.*k2)+0.1
      if(k1.eq.m1*m1.and.k2.eq.m2*m2) then
      k1=m1;k2=m2;a1='(';a2=')';end if
      Write(nu,'(f12.6,1x,a1,i6,a1,i6,a1,2x,220a1)') &
                crk(icoef),A1,iabs(k1),B1,k2,A2,(a(i:i),i=1,ia)


      End SUBROUTINE PRI_COEF1


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
      Integer, intent(in) :: i_max
      Integer, intent(out) :: k1,k2

      Integer :: i_min, i1,i2, m1,m2 
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





!======================================================================
      Subroutine Check_int(C,icase,k,j1,j2,j3,j4)
!======================================================================

      USE param_br

      Implicit none
      Integer :: icase,k,j1,j2,j3,j4,j
      Real(8) :: C

      Select case(icase)

      Case(4,5)   ! Mk,Rk

       if(j1.gt.j2) then
        j=j1;j1=j2;j2=j; j=j3;j3=j4;j4=j
       elseif(j1.eq.j2.and.j3.gt.j4) then
        j=j3;j3=j4;j4=j
       end if

      Case(3)     ! Tk

       if(j1.gt.j2) then
        j=j1;j1=j2;j2=j; j=j3;j3=j4;j4=j
       elseif(j1.eq.j2.and.j3.gt.j4) then
        j=j3;j3=j4;j4=j
       end if

       if(j1.eq.j2.and.j2.eq.j3.and.j3.eq.j4) then
        C=-C*(k-1)*(k+2)/(k+k+1)/2 ! /2 ???
        icase=4
       end if

!----------------------------------------------------------------------
!                                                        T[k] --> M[k]:
!     Here are used the relations:
!
!     U[k](a,b;c,d) + U[k](c,d;a,b) =    A(a,b;c,d) -
!
!     - (k-1)(k+2)/(2k+1) { N[k-1](a,b;c,d) + N[k-1](b,a;d,c) }
!
!     and   U[k](a,b;c,d) = T[k](a,b;c,d) + T[k](c,b;a,d).
!
!     then a=b=c=d the integrals A(a,a,a,a) disapear because
!     the T and U integrals enter only as I[k+1] - I[k-1].
!
!     In particular:  Uk(a,a,a,a)=2Tk(a,a,a,a)
!                                =-(k-1)(k+2)/(2k+1) (1/2) M[k-1]
!---------------------------------------------------------------------

      Case(9)     ! Vk(a,b;a,b)

! ???
!----------------------------------------------------------------------
!                                                        V[k] --> N[k]:
!
!     elseif(met.eq.9.and.j1.eq.j3.and.j2.eq.j4) then
!
!      S=-C*(k+1)/2
!      int=0
!      Call INT_pack(8,k-1,j1,j2,j3,j4,int)
!      i=Iadd_boef(S,int)
!      S= C*(k+0)/2
!      int=0
!      Call INT_pack(8,k+1,j1,j2,j3,j4,int)
!      i=Iadd_boef(S,int)
!
!      here is used the relation:   (R.Glass, Z.Physik,A292,131(1979))
!
!      2*V[k](r,s;r,s) = -(k+2)*N[k-1](s,r;s,r) + (k+1)N[k+1](r,s;r,s)
!
!      that is consequence of more general relation
!
!      V[k](r,s;r's') + V[k](r's';r,s) =
!                 -(k+2)*N[k-1](s,r;s',r') + (k+1)N[k+1](r,s;r',s')
!
!----------------------------------------------------------------------

       End Select

       End Subroutine Check_int



 
