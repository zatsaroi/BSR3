!======================================================================
      Subroutine fl_kpol(kpol,li,lj,fl)
!======================================================================
!     provide l-dependent factor for V-form element:
!
!     kpol=1:  <i|dv|j> = INT{ P_j(r) [d/dr + fl/r] P_i(r); dr}
!                         fl = (l_j*(l_j+1) - l_i*(l_i+1))/2
!
!     kpol=2:  <i|dv|j> = INT{ P_j(r) [r d/dr + fl] P_i(r); dr}
!                         fl = l_i + 2 for lj-li = 2
!                         fl = 1/2     for lj-li = 0
!                         fl =-l_i + 1 for lj-li =-2
!======================================================================
      Implicit none
      Integer, intent(in) :: kpol,li,lj
      Real(8), intent(out) :: fl 
      Integer :: ll       

      Select case(kpol)

      Case(1)                   ! dipole operator

       if(iabs(li-lj).ne.1) Stop 'fl_kpol: inconsistent li,lj for kpol=1'
       fl = DBLE( (lj*(lj+1)-li*(li+1))/2 )

      Case(2)                   ! quadrupole operator

       ll = lj-li;            fl = DBLE(li)
       if(ll.eq.2) then;      fl = fl + 2.d0
       elseif(ll.eq.0) then;  fl = 1.5d0
       elseif(ll.eq.-2) then; fl =-fl + 1.d0
       else; Stop 'fl_kpol: inconsistent li,lj for kpol=2'              
       end if

      Case default

       fl = 0.d0

      End Select      

      End Subroutine fl_kpol

