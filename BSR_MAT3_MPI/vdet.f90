!======================================================================
      Real(8) FUNCTION VDET (kd,N1,N2)
!======================================================================
!     calculate the value of overlap determinant for given orbitals
!     Calls:  DET
!----------------------------------------------------------------------

      USE spline_orbitals, ONLY: OBS
      USE target, only: nelc
      
      Implicit none

      Integer, Intent(in) :: kd, N1(kd), N2(kd)
      Integer :: i,j 
      Real(8) :: ADET((nelc+1)*(nelc+1))
      Real(8), External :: DET

      if(kd.eq.0) then                ! ... trivial case

       VDET = 1.d0

      elseif(kd.eq.1) then

       VDET = OBS(N1(1),N2(1))

      elseif(kd.eq.2) then

       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))
      
      elseif(kd.eq.3) then

       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2))*OBS(N1(3),N2(3)) +  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(3))*OBS(N1(3),N2(1)) +  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(1))*OBS(N1(3),N2(2)) -  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(2))*OBS(N1(3),N2(1)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))*OBS(N1(3),N2(3)) -  &
              OBS(N1(1),N2(1))*OBS(N1(2),N2(3))*OBS(N1(3),N2(2)) 

      else                
      
       Do i=1,kd;  Do j=1,kd
         adet((i-1)*kd+j)=OBS(N1(i),N2(j))
       End do; End do

       VDET = DET(kd,adet)      

      end if

      END FUNCTION VDET
