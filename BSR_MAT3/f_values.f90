!======================================================================
      Subroutine f_values
!======================================================================
!     define the f-values between target states based on the
!     given asimptotic coefficients ACF for k=1 
!----------------------------------------------------------------------
      Use bsr_mat
      Use target
      Use channel
      Use bsr_matrix
      Use zconst, only: c_au, time_au  

      Implicit none
      Real(8) :: AK(ntarg,ntarg), RDME(ntarg,ntarg)      
      Integer :: IP(ntarg,ntarg)      
      Real(8), external :: CLEBCH, Z_6jj, Reduce_factor
      Real(8) :: S,SS, g1,g2, de, a,f
      Integer :: i,j, i1,i2,nt 

      AK=0.d0; IP=0; RDME=0.d0
      Do i=1,nch-1; Do j=i+1,nch; if(i.eq.j) Cycle
       S=ACF(i,j,1)/2.d0; 
       i1=iptar(i); i2=iptar(j)
       if(abs(S).lt.eps_acf) Cycle
       SS = Reduce_factor(i,j,1)
       if(abs(SS).lt.eps_acf) Cycle
       S = S / SS
       de=Etarg(i2)-Etarg(i1)
       RDME(i1,i2) = RDME(i1,i2) + S

       if(istarg(i1).ne.0) then
        S = S*S*istarg(i1)
        g1 = (2*ltarg(i1)+1) * istarg(i1)
        g2 = (2*ltarg(i2)+1) * istarg(i2)
       else
        S = S*S        
        g1 = jtarg(i1)
        g2 = jtarg(i2)
       end if

       f = 2.d0/3.d0*de*S /g1
       a = 4.d0/3.d0*de**3*S/c_au**3/time_au /g2

       AK(i1,i2) = AK(i1,i2) + f
       AK(i2,i1) = AK(i2,i1) + a 
       IP(i1,i2) = IP(i1,i2) + 1
      End do; End do

      nt = 0
      Do i1=1,ntarg; Do i2=i1,ntarg
       if(IP(i1,i2).eq.0) Cycle
       RDME(i1,i2) = RDME(i1,i2) / IP(i1,i2)
       AK(i1,i2) = AK(i1,i2) / IP(i1,i2)
       AK(i2,i1) = AK(i2,i1) / IP(i1,i2)
       nt = nt + 1
      End do; End do

! ... total decay probabilities:
  
      Do i=2,ntarg;  AK(i,i) = SUM(AK(i,1:i-1)); End do  

! ... print results:

      if(nt.gt.0) &
      write(pri,'(/a,a/)') 'Target radiative data:  '
      write(pri,'(/a,a/)') 'f-value, A-value, branching ratio,',&
                           ' dipole reduced matrix elements'

      Do i=1,ntarg-1; Do j=i+1,ntarg
       if(AK(i,j).eq.0.d0) Cycle

       write(pri,'(2i5,1PE12.3,E12.3,0Pf10.5,E15.5,5x,a,a)') &
        i,j,AK(i,j),AK(j,i),AK(j,i)/AK(j,j), RDME(i,j),BFT(i),BFT(j)
      End do; End do

! ... polarizability of the ground state:

      a = 0.d0
      Do i = 2,ntarg
       de = Etarg(i)-Etarg(1)
       a = a + AK(1,i)/(de*de) 
      End do

      if(a.ne.0.d0) &
      write(pri,'(/a,f10.3/)') 'Polarizability of the ground state =',a

      End  Subroutine f_values

!======================================================================
      Real(8) Function Reduce_factor(ich,jch,k)
!======================================================================
!     define factor connecing reduced dipole matrix element with
!     asymptotic coefficient for multipole index k
!----------------------------------------------------------------------
      Use target
      Use channel

      Implicit none
      Integer, intent(in) :: ich,jch,k
      Integer :: ll1,ll2,jj1,jj2,it,jt,L1,L2,S1,S2,J1,J2,LT,JJ, kz
      Real(8) :: S,SS, zero = 0.d0
      Real(8), external :: ZCLKL, Z_6jj, Z_6j
       
      Reduce_factor = zero

      ll1 = lch(ich);  ll2 = lch(jch) 
      S = ZCLKL(ll1,k,ll2)
      if(S.eq.zero) Return
      jj1 = jkch(ich); jj2 = jkch(jch) 

      it = iptar(ich); jt = iptar(jch)
      L1 = ltarg(it);  L2 = ltarg(jt)
      S1 = istarg(it); S2 = istarg(jt)
      J1 = jtarg(it);  J2 = jtarg(jt) 

      LT = lpar; JJ =jpar

      if(coupling.eq.'LS') then

       if(S1.ne.S2) Return
       S = S * Z_6jj(L1,ll1,LT,ll2,L2,k) 
       kz = L2+ll1+LT
       S = S * (-1) ** kz

      elseif(coupling.eq.'JK') then

       if(jj1.ne.jj2) Return       
       S = S * Z_6j(ll2+ll2+1,ll1+ll1+1,k+k+1,J1,J2,jj1) 
       kz = J2+JJ+JJ-ll1-ll1-jj1+1; kz=kz/2
       S = S * (-1) ** kz

      elseif(coupling.eq.'JJ') then

       SS = jj1*jj2;   S = S * sqrt(SS)
       S = S * Z_6j(ll1+ll1+1,2,jj1,jj2,k+k+1,ll2+ll2+1)
       S = S * Z_6j(J1,jj1,JJ,jj2,J2,k+k+1)
       kz = J2+jj1+JJ+ll1+ll1+1+jj2+k+k; kz = kz/2
       S = S * (-1) ** kz

      end if

      Reduce_factor = S

      End Function Reduce_factor