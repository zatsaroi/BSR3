!======================================================================
      Subroutine  Get_dvl(i1,i2,cl,cv,v,w)
!======================================================================
!     provide d-vectors
!----------------------------------------------------------------------
      Use spline_param,    only: ns
      Use spline_orbitals, only: lbs
      Use bsr_dmat,        only: ktype,kpol,rbs,dbsl,tbs

      Implicit none
      Integer, intent(in)  :: i1,i2
      Real(8), intent(in)  :: cl,cv
      Real(8), intent(out) :: v(ns),w(ns)
      Real(8) :: fl

      v(1:ns) = cl*rbs(1:ns,i1)

      if(ktype.eq.'M') then; w(1:ns)= cv*rbs(1:ns,i1); Return; end if
       
      Call FL_kpol(kpol,lbs(i1),lbs(i2),fl); fl=fl*cv

      w(1:ns)=cv*dbsl(1:ns,i1)+fl*tbs(1:ns,i1)
     
      End Subroutine Get_dvl

