!======================================================================
      Subroutine Decode_conf(BSYM,mo)
!======================================================================
!     decode the ang. symmetry from 'sym' format to ZOI format
!----------------------------------------------------------------------

      USE conf_LS

      Implicit none
      Character(40), Intent(in) ::  BSYM  ! angular symmetry
      Integer :: mo,i,j
      Integer, External :: LA
      
      j=1
      Do i=1,mo
       LS(i,3) = LA(BSYM(j:j))
       LS(i,2)=2*LA(BSYM(j+1:j+1))+1
       read(BSYM(j+2:j+2),'(i1)') LS(i,1)
       LS(i,5) = LA(BSYM(j+3:j+3))
       LS(i,4)=2*LA(BSYM(j+4:j+4))+1
       j=j+5
      End do

      End Subroutine Decode_conf
