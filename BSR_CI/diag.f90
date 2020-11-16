!======================================================================
      Subroutine DIAG(jot,nu)
!======================================================================
!
!     this subroutine calls the LAPACK (htmm://www.netlib.org)
!
!     routines (DSYGV or DSYGVX) to solve 
!     the generalized eigenvalue problem
!
!                Hx = E Sx,    NORT > -1,
!
!     or simple eigenvalue problems (DSYEV or DSYEVX)
!
!                Hx = E x,     NORT = -1,
!
!     and write the result expansions on unit 'nu'
!
!     jot = 2J+1 (or = 0  for LS calculations)
!
!     MEIV - number of needed eigenvectors (in module ci_data)
!
!----------------------------------------------------------------------
      Use bsr_ci

      Implicit none
      Integer :: j,k,kk,jot,nu,info
      Real(8) :: CM, gvj, gvls
      Real(8) :: eval(ncfg)

!----------------------------------------------------------------------
!   be careful:  DSYEV (DSYGV) does not like the case of degenerate
!   states with zero interaction. In this case, we get accidentaly
!   rotated solutions! 
!
!    C = 1.d-8; EPS = 1.d-12
!    Do i=1,nzero
!     m=0; Do j=1,i-1; if(H(i,j).eq.0.d0) Cycle; m=1; Exit; End do
!     if(m.ne.0) Cycle      
!     Do ii=i+1,nzero
!      if(abs(H(ii,ii)-H(i,i)).gt.EPS) Cycle
!      mm=0; Do j=1,ii-1; if(H(ii,j).eq.0.d0) Cycle; mm=1; Exit; End do
!      if(mm.eq.0) then
!       H(ii,ii) = H(ii,ii)+C*H(ii,ii); C=C*2.d0
!      end if
!     End do
!    End do 
!
!----------------------------------------------------------------------
!     solve the  HX = E SX  problem with LAPACK routines:
!----------------------------------------------------------------------

      if(NORT.gt.-1) then   ! generalized eigenvalue problem

       if(meiv.lt.ncfg) then
        Call LAP_DSYGVX('V','L',ncfg,ncfg,HM,HS,eval,meiv,info)
       else

        Call LAP_DSYGV ('V','L',ncfg,ncfg,HM,HS,eval,info)
       end if
 
      else                  ! simple eigenvalue problem

       if(meiv.lt.ncfg) then
        Call LAP_DSYEVX('V','L',ncfg,ncfg,HM,eval,meiv,info)
       else
        Call LAP_DSYEV ('V','L',ncfg,ncfg,HM,eval,info)
       end if

      end if

!----------------------------------------------------------------------
!     PRINT OUT THE EIGENVALUES AND EIGENVECTORS:

! ... number of solusions:

      WRITE (nu,'(//A8,I4,2X,A8,I4)') '  2*J = ',jot-1,'NUMBER =',meiv      

      Do k = 1,meiv

       CM=0.d0; kk=0                    ! find the biggest component
       DO J=1,ncfg
        if(abs(HM(J,K)).gt.CM) then; CM=abs(HM(J,K)); kk=j; end if
       End do

       Call g_factor(jot,k,kk,gvJ,gvLS) 

       write(nu,'(3x,a,f15.10,4x,a,f15.10,2x,a,f15.10)') &
        'Ssms=',0.d0,'g_J=',gvj,'g_JLS=',gvls
       write(nu,'(I6,F16.8,3x,a)')  k, EVAL(k) + EC, trim(LABEL(kk))
       write(nu,'(7F11.8)') (HM(J,K),J=1,ncfg)

      End do

      End Subroutine DIAG



