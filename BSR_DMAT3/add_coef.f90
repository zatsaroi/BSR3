!======================================================================
      Subroutine Add_coef(CL,CV,j1,j2,j3)
!======================================================================
!     add new coefficient to the list for given itype
!     Calls: Add_cdata, Gen_matrix
!----------------------------------------------------------------------
      Use cmdata

      Real(8), intent(in) :: CL,CV
      Integer, intent(in) :: j1,j2,j3

! ... add coefficient to list:
      
      i = iblk(itype); ip = ipblk(i); jp = jpblk(i) 
      Call Add_cdata(ip,jp,CL,CV,j1,j2,j3)                    
      jpblk(i) = jp

! ... check if the block full:

      if(jp-ip+1.lt.mblocks) Return

! ... find new block:

      m = 0
      Do i=1,nblocks
       if(jpblk(i).ge.0) Cycle
       iblk(itype) = i; jpblk(i) = 0
       nblk(itype) = nblk(itype) + 1
       kblk(itype,nblk(itype))=i
       m = 1; Exit
      End do
      if(m.ne.0) Return
     
! ... everything is full - it is time to generate matrix for
! ... the current and most extended cases:

      Call Gen_matrix

      m = 0
      Do i=1,ntype
       if(nblk(i).lt.m) Cycle
       m=nblk(i); itype=i	    
      End do
      Call Gen_matrix

      End Subroutine Add_coef 
