!======================================================================
      Real(8) Function VDET(kd,N1,N2)
!======================================================================
!     calculate the value of overlap determinant 
!----------------------------------------------------------------------
      Use spline_orbitals, only: obs
      Use conf_LS,         only: ne
      
      Implicit none
      Integer, intent(in) :: kd, N1(kd), N2(kd)
      Integer :: i,j
      Real(8) :: ADET(ne*ne)      
      Real(8), external :: DET

      if(kd.gt.ne) Stop 'problem with kd in VDET'

! ... trivial overlap

      if(kd.eq.0) then; VDET=1.d0; Return;  End if

! ... define the overlap matrix

      Do i=1,kd; Do j=1,kd
       ADET((i-1)*kd+j) = obs(N1(i),N2(j))
      End do; End do

      VDET = DET(kd,ADET)

      End Function VDET
