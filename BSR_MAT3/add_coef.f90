!======================================================================
      Subroutine Add_coef(C,j1,j2,j3,j4)
!======================================================================
!     add new coefficient to the list in cmdata module
!----------------------------------------------------------------------
      Use cmdata

      Implicit none
      Real(8) :: C
      Integer :: j1,j2,j3,j4
      Integer :: i,k,m,n, ip,jp, jtype,jpol

! ... add coefficient to list:

      i = iblk(kpol,itype); ip = ipblk(i); jp = jpblk(i) 
      Call Add_cdata(ip,jp,C,j1,j2,j3,j4)                    
      jpblk(i) = jp

! ... check if the block full:

      if(jp-ip+1.lt.mblocks) Return

! ... find new block:

      if(nblk(kpol,itype).lt.kblocks) then
      m = 0
      Do i=1,nblocks
       if(jpblk(i).ge.0) Cycle
       iblk(kpol,itype) = i; jpblk(i) = 0
       nblk(kpol,itype) = nblk(kpol,itype) + 1
       kblk(kpol,itype,nblk(kpol,itype))=i
       m = 1; Exit
      End do
      if(m.ne.0) Return
      end if
     
! ... everything is full - it is time to generate matrix:

      jtype=itype; jpol=kpol; n = nblk(kpol,itype)
      Do i = 1,ntype; Do k = -1,npol
       if(nblk(k,i).le.n) Cycle
       jtype = i; jpol = k
      End do; End do
      Call Gen_matrix(jtype,jpol)
      if(jtype.eq.itype.and.jpol.eq.kpol) Return

! ... find again new block:

      m = 0
      Do i=1,nblocks
       if(jpblk(i).ge.0) Cycle
       iblk(kpol,itype) = i; jpblk(i) = 0
       nblk(kpol,itype) = nblk(kpol,itype) + 1
       kblk(kpol,itype,nblk(kpol,itype))=i
       m = 1; Exit
      End do
      if(m.eq.0) Stop ' Add_coef: problems with new block '
      
      End Subroutine Add_coef 
