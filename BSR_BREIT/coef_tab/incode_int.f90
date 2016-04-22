!======================================================================
      Integer(4) Function Incode_INT (met,k,I1,I2,I3,I4)
!======================================================================
!
!     incode the integral
!
!----------------------------------------------------------------------
 
      Use param_br

      Implicit none
      Integer(4), intent(in) :: met,k,I1,I2,I3,I4       

      if(max(i1,i2,i3,i4).ge.jb) Stop ' INT_pack: i > ibase '
      if(k.ge.mk) Stop 'INT_pack: k > kmax'
      if(met.ge.21) Stop 'INT_pack: met too large'

      Incode_INT = ((i1*jb+i2)*jb+i3)*jb+i4 + k*jb4 + met*jb8

      End Function Incode_INT


!======================================================================
      Subroutine Decode_INT (met,k,I1,I2,I3,I4,int)
!======================================================================
!
!     decode the integral
!
!----------------------------------------------------------------------
 
      Use param_br

      Implicit none
      Integer(4), Intent(in) :: int
      Integer(4), Intent(out) :: met,k,I1,I2,I3,I4

      K  = int;   met = k/jb8;  k = mod(k,jb8)
      I4 = mod(K,jb);  K  =  K/jb
      I3 = mod(K,jb);  K  =  K/jb
      I2 = mod(K,jb);  K  =  K/jb
      I1 = mod(K,jb);  K  =  K/jb

      End Subroutine Decode_INT
