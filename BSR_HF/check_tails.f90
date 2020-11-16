!======================================================================
      Subroutine Check_tails(jo)
!======================================================================
!     nulify too small coefficients in the end for orbital jo
!     if jo=0 - check all orbitals
!     only large component is used to define small B-splines
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals
      
      Implicit none
      Integer :: io,jo,i
      Real(8) :: cm

      Do io=1,nbf; if(jo.ne.0.and.io.ne.jo) Cycle
       cm = maxval( abs( p(:,io) ) )
       Do i=ns-1,1,-1
        mbs(io)=i     
        if(abs(p(i,io))/cm.lt.end_tol) Cycle
        Exit
       End do
       p(i+1:ns,io)=0.d0
      End do
      Call Boundary_conditions 

      End Subroutine Check_tails

       