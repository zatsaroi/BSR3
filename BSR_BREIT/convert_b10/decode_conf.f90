!======================================================================
      Subroutine Decode_conf(ASYM,BSYM)
!======================================================================
!     decode the config. from 'sym' format to ZOI format
!----------------------------------------------------------------------

      USE conf_LS

      Character(26), Intent(in) ::  ASYM  ! config. symmetry
      Character(40), Intent(in) ::  BSYM  ! angular symmetry
      Integer :: i,j,k
      Integer, External :: LA
      
       no=LEN_TRIM(ASYM(1:24))/3

       k=1; j=1
       Do i=1,no
        ln(i) = LA(ASYM(k:k))
        read(ASYM(k+1:k+2),*) iq(i) 
        k=k+3
        LS(i,3) = LA(BSYM(j:j))
        LS(i,2)=2*LA(BSYM(j+1:j+1))+1
        read(BSYM(j+2:j+2),'(i1)') LS(i,1)
        LS(i,5) = LA(BSYM(j+3:j+3))
        LS(i,4)=2*LA(BSYM(j+4:j+4))+1
        j=j+5
       End do

       End Subroutine Decode_conf
