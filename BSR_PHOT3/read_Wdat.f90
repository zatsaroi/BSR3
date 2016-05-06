!======================================================================
      Subroutine Read_Wdat (nu)
!======================================================================
!     read W.DAT file (unit nu) 
!----------------------------------------------------------------------      
      Use bsr_phot

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i,j,kch, kcp 

      if(nwt.lt.0) Return

      read(nu) kch, kcp, khm
 
      if(kch.ne.nch) Stop ' other nch in W.DAT'
      if(khm.ne.nhm) Stop ' other nhm in W.DAT'
      ncp=kcp
 
      nwt = nch + ncp;  Allocate (WT(nwt,nhm),AK(nwt,nch),WTch(nwt))
      Do i = 1,nhm; read(nu) (WT(j,i),j=1,nwt); End do
 
      Close(nu)

      End Subroutine Read_Wdat
