!======================================================================
      Integer Function Incode_INT (met,k,I1,I2,I3,I4)
!======================================================================
!     incode the integral
!----------------------------------------------------------------------
 
      Implicit none
      Integer :: ib4=2**4, ib5=2**5, ib6=2**6
      Integer, intent(in) :: met,k,I1,I2,I3,I4       
      if(max(i1,i2,i3,i4).ge.ib5) Stop 'Incode_INT: i >= 2**5'
      if(met.gt.ib4) Stop 'Incode_INT: met >= 2**4'
      if(k.gt.ib6) Stop 'Incode_INT: met >= 2**6'

      Incode_INT = ((((k*ib4+met)*ib5+i4)*ib5+i3)*ib5+i2)*ib5+i1

      End Function Incode_INT

!======================================================================
      Subroutine Decode_INT (met,k,I1,I2,I3,I4,int)
!======================================================================
!     decode the integral
!----------------------------------------------------------------------
 
      Implicit none
      Integer :: ib4=2**4, ib5=2**5
      Integer, Intent(in)  :: int
      Integer, Intent(out) :: met,k,I1,I2,I3,I4

      k   = int
      I1  = mod(k,ib5);  k = k/ib5
      I2  = mod(k,ib5);  k = k/ib5
      I3  = mod(k,ib5);  k = k/ib5
      I4  = mod(k,ib5);  k = k/ib5
      met = mod(k,ib4);  k = k/ib4

      End Subroutine Decode_INT


!======================================================================
      Subroutine Decode_met(met,int)
!======================================================================
!     decode the integral
!----------------------------------------------------------------------
 
      Implicit none
      Integer :: ib4=2**4, ib20=2**20
      Integer, Intent(in)  :: int
      Integer, Intent(out) :: met

      met = int/ib20; met = mod(met,ib4) 

      End Subroutine Decode_met
