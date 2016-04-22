!======================================================================
      Subroutine ZNO_breit( no1,l1,m1,ms1,no2,l2,m2,ms2,  &
                            no3,l3,m3,ms3,no4,l4,m4,ms4,  ide)
!======================================================================
!
!     computes angular coefficients for the two-electron Breit-Pauli
!     operators in uncouple nlms-representation.
!
!     ide = +1  ->  direct interaction
!     ide = -1  ->  exchange interaction
!
!----------------------------------------------------------------------

      USE bsr_breit, only: mk, joper

      Implicit none

      Integer, Intent(in) :: no1,l1,m1,ms1,no2,l2,m2,ms2,  &
                             no3,l3,m3,ms3,no4,l4,m4,ms4,  ide

      Integer :: q,qq, k,kk, k1,k2, kq,kz
      Integer :: KL,KM, KL1,KL2, KM1,KM2, met

      Real(8) :: A1(0:mk+3),A2(0:mk+3),B1(0:mk+3),B2(0:mk+3)

      Real(8) :: C, C1,C2, CQ, S,SS,S1,S2,A,AB, FS1,FS2 
      Real(8) :: D0 = 0.d0
      Real(8), External :: Z_3jj

!---------------------------------------------------------------------
! ... define the range of multipole indeces and common multipliers:
    
      KL1=IABS(l1-l3);  KL2=IABS(l2-l4);  KL=MAX0(KL1,KL2)
      KM1= l1+l3;       KM2=l2+l4;        KM=MIN0(KM1,KM2)
      if(KM.lt.KL) Return

      S = (l1+l1+1)*(l2+l2+1)*(l3+l3+1)*(l4+l4+1)
      S = sqrt(S) * ide * (-1)**(m3+m4)
      q = m1-m3; kq = (-1)**q

!---------------------------------------------------------------------
! ... define all relevant 3j-symbols: 

      DO kk = 0,KM-KL+3
       k = KL + kk - 1; if(k.gt.mk) Cycle
       A1(kk)=D0; B1(kk)=D0
       if(k.ge.KL1.and.k.le.KM1) then
        if(mod(l1+l3+k,2).eq.0) A1(kk) = Z_3jj(l1,0,l3,0,k,0)
        B1(kk) = Z_3jj(l1,-m1,l3,m3,k,m1-m3)
       end if
       A2(kk)=D0; B2(kk)=D0
       if(k.ge.KL2.and.k.le.KM2) then
        if(mod(l2+l4+k,2).eq.0) A2(kk) = Z_3jj(l2,0,l4,0,k,0)
        B2(kk) = Z_3jj(l2,-m2,l4,m4,k,m2-m4)
       end if
      End do
!----------------------------------------------------------------------
! ... cycle on multypoles:
  
      Do kk = 1,km-kl+1,2                         

       k = kl + kk - 1; if(k.gt.mk) Cycle

       AB = A1(kk)*B1(kk)*A2(kk)*B2(kk)*S *kq      ! dk * dk * (-1)^q

!----------------------------------------------------------------------
! ... e-e interaction:

      C = AB
      if(joper(3).gt.0.and.ms1.eq.ms3.and.ms2.eq.ms4) then

       met=5; Call Add_ci(C,met,k,no1,no2,no3,no4)  ! R[k](r,s;r',s')

      end if

!---------------------------------------------------------------------
! ... o-o interaction (with Uk --> Tk + Tk'):

      if(joper(7).ne.0.and.ms1.eq.ms3.and.ms2.eq.ms4) then
 
       C = -AB
       if(k.gt.0.and.C.ne.D0) then               

       S1=l1*(l1+1)-l3*(l3+1)-k*(k+1)
       S2=l2*(l2+1)-l4*(l4+1)-k*(k+1)
       SS=S1*S2
 
       met = 3; A = 2*k*(k+1)*C

       CALL ADD_ci( A,met,k+1,no1,no2,no3,no4)      ! T[k+1](r,s;r',s')
       CALL ADD_ci(-A,met,k-1,no1,no2,no3,no4)      ! T[k-1](r,s;r',s')
 
!      for Uk-inegrals, we use the replacement:
!      U(1,2;3,4) = T(1,2;3,4) + T(3,2;1;4)
!      [W.Dankwort, J.Phys.B10,L369(1977)]

       A = S1*C
       CALL ADD_ci( A,met,k+1,no1,no2,no3,no4)      ! U[k+1](r,s;r',s')
       CALL ADD_ci( A,met,k+1,no3,no2,no1,no4)
       CALL ADD_ci(-A,met,k-1,no1,no2,no3,no4)      ! U[k-1](r,s;r',s')
       CALL ADD_ci(-A,met,k-1,no3,no2,no1,no4)
 
       A = S2*C
       CALL ADD_ci( A,met,k+1,no2,no1,no4,no3)      ! U[k+1](s,r;s',r')
       CALL ADD_ci( A,met,k+1,no4,no1,no2,no3)   
       CALL ADD_ci(-A,met,k-1,no2,no1,no4,no3)      ! U[k-1](s,r;s',r')
       CALL ADD_ci(-A,met,k-1,no4,no1,no2,no3)   
 
       met = 4;  A = C*SS*(k-2)/k/(k+k-1)   ! /2 ?

       CALL ADD_ci( A,met,k-1,no1,no2,no3,no4)      ! M[k-2](r,s;r',s')
 
       A = -C*SS*(k+3)/(k+1)/(k+k+3)  ! /2 ?

       CALL ADD_ci( A,met,k+1,no1,no2,no3,no4)      ! M[k](r,s;r',s')
 
       end if   ! (k>0)
 
       met=4
       C = -A1(kk)*B1(kk+1)*A2(kk)*B2(kk+1)*S*kq
       if(C.ne.D0) then
        A =(l1+l3+k+2)*(l3-l1+k+1)*(l1-l3+k+1)*(l1+l3-k)*   &
           (l2+l4+k+2)*(l4-l2+k+1)*(l2-l4+k+1)*(l2+l4-k)
        A=C*sqrt(A)/((k+1)*(k+2))*2
        CALL ADD_ci(A,met,k+1,no1,no2,no3,no4)      ! M[k](r,s;r',s')
       end if

      end if

!----------------------------------------------------------------------
! ... s-s interaction:

      if(joper(6).ne.0) then

       met=10

       kz=(ms1+ms2)/2;  kz=kq*(-1)**(kz+1)
       k1=k+1;  k1=k1*k1;  k2=k+2;  k2=k2*k2;  qq=q*q
       A=(k1-qq)*(k2-qq); A=sqrt(A); A=A*S*kz

       C  = A1(kk+2)*B1(kk+2)*A2(kk)*B2(kk)*A       ! N[k](r,s;r',s')
       CALL ADD_ci(C,met,k+1,no1,no2,no3,no4)

       C  = A1(kk)*B1(kk)*A2(kk+2)*B2(kk+2)*A       ! N[k](s,r;s',r')
       CALL ADD_ci(C,met,k+1,no2,no1,no4,no3)

      end if

!----------------------------------------------------------------------
! ... s-o-o interaction :                                   

      if(joper(5).ne.0.and.ms1.eq.ms3.and.ms2.eq.ms4) then

       FS1=(ms1+ms2+ms2-3); FS2=(ms2+ms1+ms1-3)

       met=9                                        ! V[k-1](r,s;r',s')
       CQ = D0;  if(iabs(q).le.k)  CQ = q

       if(k.gt.0) then
        C  = AB*CQ
        C1 = -C*FS1; CALL ADD_ci(C1,met,k,no1,no2,no3,no4)
        C2 =  C*FS2; CALL ADD_ci(C2,met,k,no2,no1,no4,no3)
       end if

       met=8                                        ! N[k](r,s;r',s')

       if(k*l3.gt.0) then                            
        C = C/(k+1)/2
        C1 = l1*(l1+1)-l3*(l3+1)-k*(k+1)      
        C1 = -C1*C*FS1
        CALL ADD_ci(C1,met,k+1,no1,no2,no3,no4)     ! N[k]
        C1 = -C1*(k+1)/k
        CALL ADD_ci(C1,met,k-1,no2,no1,no4,no3)     ! N[k-2]
       end if

       if(k*l4.gt.0) then                            
        C2 = l2*(l2+1)-l4*(l4+1)-k*(k+1)   
        C2 = C2*C*FS2
        CALL ADD_ci(C2,met,k+1,no2,no1,no4,no3)     ! N[k]
        C2= -C2*(k+1)/k
        CALL ADD_ci(C2,met,k-1,no1,no2,no3,no4)     ! N[k-2]
       end if
 
!----------------------------------------------------------------------

       CQ = D0; if(iabs(q).le.k) CQ = (k+q+1)*(k-q+1)
       CQ = sqrt(CQ)* kq /(2*(k+1))

       C  = A1(kk)*B1(kk+1)*A2(kk)*B2(kk)*S*CQ
       C1 = (l1+l3+k+2)*(l1+l3-k)*(l3-l1+k+1)*(l1-l3+k+1)
       C1 = sqrt(C1)*C*FS1
       if(l3.gt.0) CALL ADD_ci(C1,met,k+1,no1,no2,no3,no4) ! N[k]

       C  = A1(kk)*B1(kk)*A2(kk)*B2(kk+1)*S*CQ
       C2 = (l2+l4+k+2)*(l2+l4-k)*(l4-l2+k+1)*(l2-l4+k+1)
       C2 = sqrt(C2)*C*FS2
       if(l4.gt.0) CALL ADD_ci(C2,met,k+1,no2,no1,no4,no3) ! N[k]

       if(k.gt.0) then

       CQ = D0
       if(iabs(q).le.k) CQ = (k+q)*(k-q); CQ = -sqrt(CQ)*kq /(2*k)

       C  = A1(kk)*B1(kk-1)*A2(kk)*B2(kk)*S*CQ
       C1 = (l1+l3+k+1)*(l1+l3-k+1)*(l3-l1+k)*(l1-l3+k)
       C1 = sqrt(C1)*C*FS1
       CALL ADD_ci(C1,met,k-1,no2,no1,no4,no3)             ! N[k-2]

       C  = A1(kk)*B1(kk)*A2(kk)*B2(kk-1)*S*CQ
       C2 = (l2+l4+k+1)*(l2+l4-k+1)*(l4-l2+k)*(l2-l4+k)
       C2 = sqrt(C2)*C*FS2
       CALL ADD_ci(C2,met,k-1,no1,no2,no3,no4)             ! N[k-2]

       end if   ! on k>0

      end if

!---------------------------------------------------------------------
      End do   ! over k

      End Subroutine ZNO_breit


!======================================================================
      Subroutine Add_ci(C,met,k,j1,j2,j3,j4)
!======================================================================
!     control the recording of one resulting integral
!
!     Calls: Incode_int, Iadd_boef
!----------------------------------------------------------------------

      USE bsr_breit, only: Eps_c, mk

      Implicit none
      Integer, Intent(in) :: met,k,j1,j2,j3,j4
      Real(8), Intent(in) :: C
      Integer :: int
      Integer, External :: Incode_int

      if(abs(C).lt.Eps_c) Return;  if(k.gt.mk) Return

      int = Incode_int(met,k,j1,j2,j3,j4); Call Iadd_boef(C,int)

      End Subroutine Add_ci
