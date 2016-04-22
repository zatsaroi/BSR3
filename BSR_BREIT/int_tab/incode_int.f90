!======================================================================
      Integer(4) Function Incode_INT (itype,k,I1,I2,I3,I4)
!======================================================================
!
!     incode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb,jb4,jb8

      Implicit none
      Integer(4), intent(in) :: itype,k,I1,I2,I3,I4       

      if(max(i1,i2,i3,i4).ge.jb) Stop ' Incode_int: i > ibase '
      if(itype.ge.21) Stop 'Incode_int: itype too large'

      Incode_INT = ((i1*jb+i2)*jb+i3)*jb+i4 + k*jb4 + itype*jb8

      End Function Incode_INT



!======================================================================
      Subroutine Decode_INT (itype,k,I1,I2,I3,I4,int)
!======================================================================
!
!     decode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb,jb8

      Implicit none
      Integer(4), Intent(in) :: int
      Integer(4), Intent(out) :: itype,k,I1,I2,I3,I4

      K  = int;   itype = k/jb8;  k = mod(k,jb8)
      I4 = mod(K,jb);  K  =  K/jb
      I3 = mod(K,jb);  K  =  K/jb
      I2 = mod(K,jb);  K  =  K/jb
      I1 = mod(K,jb);  K  =  K/jb

      End Subroutine Decode_INT



!======================================================================
      Integer(4) Function Incode_mult (itype,i1,i2)
!======================================================================
!
!     incode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb

      Implicit none
      Integer(4), intent(in) :: itype,i1,i2       

      if(max(i1,i2).ge.jb) then
       write(*,'(a,2i5,a,i3)')  &
               ' Incode_mult: i1,i2 =',i1,i2,'  > base =',jb
       Stop
      end if
      if(itype.gt.3) then
       write(*,'(a,i5,a)')  &
               ' Incode_mult: itype =',itype,'  is out of range (3)'
       Stop 
      end if

      Incode_mult = (itype*jb+i1)*jb+i2

      End Function Incode_mult


!======================================================================
      Subroutine Decode_mult (itype,i1,i2,int)
!======================================================================
!
!     decode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb

      Implicit none
      Integer(4), Intent(in) :: int
      Integer(4), Intent(out) :: itype,i1,i2 
      Integer(4) :: k

      k  = int
      i2 = mod(k,jb);  k = k/jb
      i1 = mod(k,jb);  k = k/jb
      itype = mod(k,jb)

      End Subroutine Decode_mult
