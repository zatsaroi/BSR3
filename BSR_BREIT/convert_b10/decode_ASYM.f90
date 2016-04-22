!======================================================================
      Subroutine Decode_ASYM(ASYM)
!======================================================================
!     decode the config. from 'sym' format to ZOI format
!----------------------------------------------------------------------

      USE conf_LS

      Implicit none

      Character(26), Intent(in) ::  ASYM  
      Integer :: i,k
      Integer, External :: LA
      
      no=LEN_TRIM(ASYM(1:24))/3

      k=1
      Do i=1,no
       ln(i) = LA(ASYM(k:k))
       read(ASYM(k+1:k+2),*) iq(i) 
       k=k+3
      End do

      End Subroutine Decode_ASYM
