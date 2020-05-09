!======================================================================
      Subroutine Do_coef(jtype)
!======================================================================
!     process the coefficients for given type
!----------------------------------------------------------------------
      Use bsr_breit, only : nub
      Use c_blocks

      Implicit none
      Integer, intent(in) :: jtype
      Integer :: i,icase,kpol,itype

      if(ncdata.le.0) Return
      Call Decode_type(jtype,icase,kpol,itype)

      write(nub)  icase,kpol,itype,ncdata

      write(nub)  (cdata(ipt(i)),i=1,ncdata)
      write(nub)  (k1(ipt(i)),i=1,ncdata)
      write(nub)  (k2(ipt(i)),i=1,ncdata)
      write(nub)  (k3(ipt(i)),i=1,ncdata)
      write(nub)  (k4(ipt(i)),i=1,ncdata)

      End Subroutine Do_coef
