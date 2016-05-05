!======================================================================
      Subroutine Check_det(kdn,N1,N2,iext)
!======================================================================
!     Expend the total overlap determinants so that to extract 
!     matrix elements with continuum orbitals.
!     Converted overlaps are stored in module new_dets with Iadd_ndets.
!----------------------------------------------------------------------
      Use spline_orbitals, only: iech, IBORT 
      Use conf_LS,         only: ne

      Implicit none
      Integer, intent(in) :: kdn, iext
      Integer, intent(in) :: N1(kdn),N2(kdn)
      Integer :: N3(ne),N4(ne),N5(ne),N6(ne)
      Integer :: i,j, j1,j2,j3,j4, io,jo, kdn1,kdn2
      Real(8) :: S      
      Real(8), external :: VDET

!----------------------------------------------------------------------
      if(kdn.gt.ne) Stop 'Check_det: dimension > ne'

! ... define the cont.orb.

      j1=0
      j2=0
      Do i=1,kdn
       if(iech(N1(i)).ne.0) j1=i
       if(iech(N2(i)).gt.0) j2=i
      End do

      if(j1.eq.0.and.j2.eq.0) then               ! no cont.orb

       S = VDET(kdn,N1,N2)**iext
       Call Iadd_ndets(0,0,S)  

      elseif(kdn.eq.1) then                      ! <kl|nl>

       if(iext.gt.1) Stop ' Check_det: iext > 1'
          
       io=IBORT(N1(1),N2(1))

       if(io.ne.0) Call Iadd_ndets(io,0,1.d0)

      elseif(j1.gt.0.and.j2.eq.0) then           ! <kl ...|...>

       if(iext.gt.1) Stop ' Check_det: iext > 1'

        Do j=1,kdn
         io=IBORT(N1(j1),N2(j))        
         if(io.eq.0) Cycle
         Call Shift(kdn,j1,N1,N3)
         Call Shift(kdn,j ,N2,N4)
         kdn1=kdn-1
         S = VDET(kdn1,N3,N4)*(-1)**(j1+j)
         Call Iadd_ndets(io,0,S)
        End do

      elseif(j2.gt.0.and.j1.eq.0) then          ! <...|... kl>

        if(iext.gt.1) Stop ' Check_det: iext > 1'

         Do j=1,kdn
          io=IBORT(N1(j),N2(j2))        
          if(io.eq.0) Cycle
          Call Shift(kdn,j ,N1,N3)
          Call Shift(kdn,j2,N2,N4)
          kdn1=kdn-1
          S = VDET(kdn1,N3,N4)*(-1)**(j2+j)
          Call Iadd_ndets(io,0,S)
         End do

      elseif(j1.gt.0.and.j2.gt.0) then          !  < kl ... | ... kl>

        if(iext.gt.1) Stop ' Check_det: iext > 1'

         Do j=1,kdn

          Call Shift(kdn,j1,N1,N3)
          Call Shift(kdn, j,N2,N4)
          io=IBORT(N1(j1),N2(j))

          if(j.eq.j2) then

            kdn1=kdn-1
            if(io.eq.0) Stop ' Orthogonal kl ?'
            S=VDET(kdn1,N3,N4)*(-1)**(j1+j2)
            Call Iadd_ndets(io,0,S)

          else

            j3=j2;if(j.lt.j2) j3=j2-1
            kdn1=kdn-1
            Do j4=1,kdn1
             jo=IBORT(N3(j4),N4(j3))
             if(io.eq.0.or.jo.eq.0) Cycle
             Call Shift(kdn1,j3,N4,N6)
             Call Shift(kdn1,j4,N3,N5)
             kdn2=kdn-2
             S=VDET(kdn2,N5,N6)*(-1)**(j1+j+j3+j4)
             Call Iadd_ndets(io,jo,S)
            End do

           end if

          End do

      end if

      End Subroutine Check_det


!======================================================================
      Subroutine Shift (n,m,N1,N2)
!======================================================================
!     N2 obtained from N1 by deleting element 'm'
!----------------------------------------------------------------------
      Implicit none 
      Integer, intent(in)  :: n,m, N1(n)
      Integer, intent(out) :: N2(n)
      Integer :: i,k

      k=0; Do i=1,n; if(i.eq.m) Cycle; k=k+1; N2(k)=N1(i); End do

      End Subroutine Shift
