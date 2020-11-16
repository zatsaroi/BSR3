!=======================================================================
      Real(8) Function a(i,j,k)
!=======================================================================
!     coefficient to the potential for electron i of y^k (i,j)
!-----------------------------------------------------------------------
      Use bsr_mchf

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
      Use bsr_mchf

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
      Use bsr_mchf

      Implicit none
      Integer, intent(in) :: jj
      Integer :: i,j, k
      Real(8) :: c
      Real(8), external :: zclkl   !, WW 
	  
      coef = 0.d0
      if(jj.le.0) Return

      Do i = 1,nbf
      Do j = 1,i
       if(j.gt.jj) Exit
       if(i.eq.j) then
        c =  1.d0    !qsum(i)*(qsum(i)-1.d0)/2   ! WW(i,j)/2.d0     
        coef(i,i,0) = c
        Do k=2,lbs(i)+lbs(i),2
         coef(i,i,k)= -c * zclkl(lbs(i),k,lbs(i))**2 /((2*lbs(i)+1)*(4*lbs(i)+1))
        End do
       else
        c = 1.d0     !qsum(i)*qsum(j)   !  WW(i,j)      
        coef(i,j,0)=c
        Do k=iabs(lbs(i)-lbs(j)),(lbs(i)+lbs(j)),2
         coef(j,i,k)= -c * zclkl(lbs(i),k,lbs(j))**2 / ((2*lbs(i)+1)*(2*lbs(j)+1)) /2           
        End do
       end if

      End do; End do
	  
      End Subroutine av_energy_coef 

