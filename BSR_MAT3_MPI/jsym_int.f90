!======================================================================
      Subroutine Jsym_int(met,i1,i2,i3,i4)
!======================================================================
!     use integral symmetry to obtaine the 'canonical' form
!----------------------------------------------------------------------

      m=1; ii=min(i1,i2,i3,i4)

      Select case(met)
       Case(4,5)                              !   M-,R-integrals
        if(ii.eq.i4) m = 4
        if(ii.eq.i3) m = 3
        if(ii.eq.i2) m = 2
        if(ii.eq.i1) m = 1
       Case(3)                                !   T - integrals
        if(ii.eq.i2) m = 2
        if(ii.eq.i1) m = 1
       Case(6,7,8,10)                         !   L, Z, N integrals
        if(ii.eq.i3) m = 3
        if(ii.eq.i1) m = 1
      End select

      if(m.eq.1) return

      j1 = i1; j2 = i2; j3 = i3; j4 = i4

      if(m.eq.2) then
        i1 = j2; i2 = j1; i3 = j4; i4 = j3
      elseif(m.eq.3) then
        i1 = j3; i2 = j4; i3 = j1; i4 = j2
      elseif(m.eq.4) then
        i1 = j4; i2 = j3; i3 = j2; i4 = j1
      end if

      End Subroutine Jsym_int

