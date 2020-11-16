!=======================================================================
      Module hf_energy
!=======================================================================
!     This module defines the energy expression as sum  of
!     one-electron       Sum_i  qsum(i) valueI (i)   
!     and two-electron   Sum_i,j,k  coef(i,j,k)  value(i,j,k)  
!     integrals.   value(i,j,k) ->  Fk(i,j) or Gk(i,j) integrals,
!     which is defined by subroutines a(i,j,k) and b(i,j,k). 
!     Note that coef(i,j,k) with i > = j refer to direct Fk-integrals
!     whereas coef(i,j,k) with i < j refer to exchage Gk-integrals. 
!----------------------------------------------------------------------
      Implicit none

      Real(8), allocatable :: coef (:,:,:), value(:,:,:), valueI(:)

      End Module hf_energy


!=======================================================================
      Real(8) Function a(i,j,k)
!=======================================================================
!     coefficient to the potential for electron i of y^k (i,j)
!-----------------------------------------------------------------------
      Use hf_energy
      Use hf_orbitals, only: lbs

      Implicit none
      Integer, intent(in)  :: i,j,k
      a = 0.d0
      if(mod(k,2).ne.0) Return
      if(k > (lbs(i)+lbs(j)) ) Return
      if(k < 0) Return
      a = coef(max(i,j),min(i,j),k)

      End Function a


!=======================================================================
      Real(8) Function b(i,j,k) 
!=======================================================================
!     determine the coefficient of the y^k (i,j) p(j) term in the exchange
!     expression for electron i
!-----------------------------------------------------------------------
      Use hf_energy
      Use hf_orbitals, only: lbs

      Implicit none
      Integer, Intent(IN) :: i,j,k

      b = 0.d0
      if(i.eq.j) Return
      if(mod(k+lbs(i)+lbs(j),2) /= 0) Return
      if( k < abs(lbs(i)-lbs(j)) ) Return
      if( k >    (lbs(i)+lbs(j)) ) Return 
      b = coef(min(i,j),max(i,j),k)
    
      End Function b


!======================================================================
      Subroutine Def_energy_coef
!======================================================================
!     define angular coefficients for energy expression
!----------------------------------------------------------------------
      Use hf_energy
      Use bsr_hf
      Use hf_orbitals, only: nbf,lbs

      kmax = 2*maxval(lbs(1:nbf))
      if(allocated(coef)) Deallocate(coef);  Allocate(coef(nbf,nbf,0:kmax))
      coef = 0.d0
      if(allocated(value)) Deallocate(value);  Allocate(value(nbf,nbf,0:kmax))
      value = 0.d0
      if(allocated(valueI)) Deallocate(valueI);  Allocate(valueI(nbf))
      valueI = 0.d0


      if(term.ne.'AV') then                                 ! ???
       if(ncore.gt.0) Call av_energy_coef(ncore)
       Call LS_energy_coef(0)
      else
       Call av_energy_coef(nbf)
      end if

      End Subroutine Def_energy_coef


!======================================================================
      Subroutine av_energy_coef(jj)
!======================================================================
!     Angular coefficients for energy in average-term approximation
!
!     Into shell, only direct interaction:  
!
!     q(q-1)/2 [ F0(a,a)  -  (2j+1)/2j  {l_a k l_a} ^2  Fk(a,a)  ]
!                                       { 0  0  0 }
!
!     Between shells: direct only k=0
!
!     q_a*q_b/2 [ F0(a,b) - 1/2  {l_a  k l_b} ^2  Gk(a,b)   ]
!                                { 0   0  0 }            
!
!     We devide on (2la+1)(2lb+1) due to using Function Cjkj instead 3J-symbol 
!
!     Input parameter jj equal to ncore or nbf, depending what we need:
!     core energy or total energy 
!----------------------------------------------------------------------
      Use hf_energy
      Use hf_orbitals

      Implicit none
      Integer, intent(in) :: jj
      Integer :: i,j, k
      Real(8) :: c
      Real(8), external :: zclkl, WW 
	  
      coef = 0.d0
      if(jj.le.0) Return

      Do i = 1,nbf
      Do j = 1,i
       if(j.gt.jj) Exit
       if(i.eq.j) then
        c = WW(i,j)/2.d0       ! qsum(i)*(qsum(i)-1.d0)/2
        coef(i,i,0) = c
        Do k=2,lbs(i)+lbs(i),2
         coef(i,i,k)= -c * zclkl(lbs(i),k,lbs(i))**2 /((2*lbs(i)+1)*(4*lbs(i)+1))
        End do
       else
        c = WW(i,j)            ! qsum(i)*qsum(j)
        coef(i,j,0)=c

        Do k=iabs(lbs(i)-lbs(j)),(lbs(i)+lbs(j)),2

         coef(j,i,k)= -c * zclkl(lbs(i),k,lbs(j))**2 / ((2*lbs(i)+1)*(2*lbs(j)+1)) /2            
        End do
       end if

      End do; End do

	  
      End Subroutine av_energy_coef 


!======================================================================
      Real(8) Function WW(i,j)
!======================================================================
! ... generalized value of WC(ic)*WC(jc) - product of expansion
! ... coeficients for configurations ic and jc
!======================================================================
      Use bsr_hf
      Use hf_orbitals
      
      Implicit none
      Integer, Intent(in) :: i,j 
      Integer :: ic

      if(nconf.eq.1) then
       if(i.eq.j) WW = qsum(i)*(qsum(i)-1.d0)       
       if(i.ne.j) WW = qsum(i)*qsum(j)
       Return
      end if

      WW = 0.d0
      if(i.eq.j) then
       Do ic=1,nconf
        WW = WW + iqconf(ic,i)*(iqconf(ic,i)-1)*weight(ic)
       End do
      else
       Do ic=1,nconf
        WW = WW + iqconf(ic,i)*iqconf(ic,j)*weight(ic)
       End do
      end if

      End Function WW


!======================================================================
      Subroutine update_int(ii)
!======================================================================
!     Update the integral list for orbital "ii"
!     if ii=0  -  all orbitals
!----------------------------------------------------------------------
      Use hf_energy
      Use hf_orbitals, only: nbf, l => lbs

      Implicit none
      Integer, intent(in) :: ii
      Integer :: i,j,k
      Real(8), external :: I_int, a, b, rk

      Do i = 1,nbf

       if(ii.ne.0.and.i.ne.ii) Cycle

       valueI(i) = I_int (i,i)
       Do j = 1,i   ! nbf 
        Do  k = 0,2*min(l(i),l(j)),2
         if(a(i,j,k).eq.0.d0) Cycle
         value(i,j,k) = rk(i,j,i,j,k)
        End do
       End do 
 
       Do  j = 1,i-1  ! nbf 
        Do k = abs(l(i)-l(j)),l(i)+l(j),2
         if(b(i,j,k).eq.0.d0) Cycle
         value(j,i,k) = rk(i,j,j,i,k)
        End do
       End do

      End do

      End Subroutine update_int


!======================================================================
      Subroutine energy
!======================================================================
!     Determine the total energy of the state, together with                                                                
!     one-electron and two-electron (direct and exchange) parts
!----------------------------------------------------------------------
      Use hf_energy
      Use bsr_hf,     only: Etotal, E1body, E2body
      Use hf_orbitals, only: nbf, l => lbs, qsum

      Implicit none
      Integer :: i,j,k
      Real(8) :: c, S
      Real(8), external :: a,b

! ... direct interaction and one-electron energies

      Etotal = 0.d0; E1body =0.d0; E2body =0.d0
      Do i = 1,nbf
       S = valueI(i)
       etotal = etotal + qsum(i)*S
       E1body = E1body + qsum(i)*S
       Do j = 1,i
        Do  k = 0,2*min0(l(i),l(j)),2
         c = a(i,j,k); if(c.eq.0.d0) Cycle
         S = value(i,j,k)                        
         etotal = etotal + c*S
         E2body = E2body + c*S
        End do
       End do 

! ... exchange interaction
 
       Do  j = 1,i-1
        Do k = abs(l(i)-l(j)),l(i)+l(j),2
         c = b(i,j,k); ; if(c.eq.0.d0) Cycle
         S = value(j,i,k)                         
         etotal = etotal + c*S
         E2body = E2body + c*S
        End do
       End do

      End do  !  over orbitals, i

      End Subroutine energy


!======================================================================
      Subroutine LS_energy_coef(jc)
!======================================================================
!     This routine gets agular coefficients for set of atomic states
!     in LS coupling (atomic states are taken from c-file, unit nuc)
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_energy
      Use hf_orbitals

      Implicit none
      Integer :: i,j,k, ic,jc, ii,jj
      Real(8) :: c
      Real(8), allocatable :: coefs(:,:,:)
      Integer, external :: Ifind_orb 

      rewind(nus)
      ic = 0
      Do 
       read(nus,'(a)') CONFIG
       if(CONFIG(5:5).ne.'(') Cycle
       read(nus,'(a)') COUPLE
       ic = ic + 1
       if(jc.ne.0.and.ic.ne.jc) Cycle

       Call Decode_cn(CONFIG,COUPLE,no,nn,ln,iq,in,LS) 

       Allocate(coefs(no,no,0:kmax))

       Call Coef_ee_1conf(no,ln,iq,LS,kmax,coefs) 

       Do i=1,no;  ii = Ifind_orb(nn(i),ln(i),0)
       Do j=1,no;  jj = Ifind_orb(nn(j),ln(j),0)
       Do k=0,kmax
        C = coefs(i,j,k); if(abs(C).lt.eps_c) C=0.d0
        if(C.eq.0.d0) Cycle
        coef(ii,jj,k) = coef(ii,jj,k) + C * weight(ic)  
       End do; End do; End do
       Deallocate(coefs)

       if(jc.ne.0.and.ic.eq.jc) Exit
       if(ic.eq.nconf) Exit
      End do  ! over configurations
  
      End Subroutine LS_energy_coef 



!======================================================================
      Subroutine core_energy
!======================================================================
!     Determine the total energy of the state, together with                                                                
!     one-electron and two-electron (direct and exchange) parts
!----------------------------------------------------------------------
      Use hf_energy
      Use bsr_hf,      only: ncore,Ecore
      Use hf_orbitals, only: nbf, l => lbs, qsum

      Implicit none
      Integer :: i,j,k
      Real(8) :: C,S
      Real(8), external :: a,b

! ... direct interaction and one-electron energies

      Ecore = 0.d0
      Do i = 1,ncore
       S = valueI(i)
       Ecore = Ecore + qsum(i)*S
       Do j = 1,i
        Do  k = 0,2*min0(l(i),l(j)),2
         c = a(i,j,k); if(c.eq.0.d0) Cycle
         S = value(i,j,k)                        
         Ecore = Ecore + c*S
        End do
       End do 

! ... exchange interaction
 
       Do  j = 1,i-1
        Do k = abs(l(i)-l(j)),l(i)+l(j),2
         c = b(i,j,k); ; if(c.eq.0.d0) Cycle
         S = value(j,i,k)                         
         Ecore = Ecore + c*S
        End do
       End do

      End do  !  over orbitals, i

      End Subroutine core_energy
