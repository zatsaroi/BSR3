!======================================================================
      Subroutine make_coupling
!======================================================================
!     prepare coupling scheme for state in module conf_LS: J1_coupling
!----------------------------------------------------------------------
      Use bsr_conf,  JJ_coupling => J1_coupling
      Use conf_LS,   n=>no

      if(n.gt.msh) Stop 'make_coupling: no > mshells'

! ... JJ-coupling:

      JJ_coupling(1,1) = 1
      JJ_coupling(2,1) =   2
      JJ_coupling(3,1) = n+2
	  
      Do i=2,n-1
       JJ_coupling(1,i) = n+i
       JJ_coupling(2,i) =   i+1
       JJ_coupling(3,i) = n+i+1
      End do

! ... moments

      Do i=1,n
       momentS(i)  =  LS(i,3)        
       momentS(i+n) = LS(i,5)
       momentL(i)  =  LS(i,2)        
       momentL(i+n) = LS(i,4)
      End do

      ncup = n-1
      nmom = 3*n

      End Subroutine make_coupling


!======================================================================
      Subroutine make_coupling_insert(ii)
!======================================================================
!     prepare coupling scheme for state in module conf_jj
!     with outer orbital is going to position ii: J2_coupling
!----------------------------------------------------------------------

      Use bsr_conf,  JJ_coupling => J2_coupling   

      Use conf_LS,   n=>no

      Integer :: i,ii, ipos(msh)
 
      if(n.gt.msh) Stop 'make_coupling: no > mshells'

      Do i=1,n
       if(i.lt.ii) then;       ipos(i)=i
       elseif(i.eq.ii) then;   ipos(i)=n
       else;                   ipos(i)=i-1
       end if
      End do

! ... JJ-coupling:

      m = n + n
      JJ_coupling(1,1) = ipos(1)
      JJ_coupling(2,1) = ipos(2)
      JJ_coupling(3,1) = m+2
	  
      Do i=2,n-1
       JJ_coupling(1,i) = i+m
       JJ_coupling(2,i) = ipos(i+1)
       JJ_coupling(3,i) = i+m+1
      End do

! ... moments

      Do i=1,n
       momentS(i+m) = LS(i,5)
       momentL(i+m) = LS(i,4)
      End do

      JJ_coupling(3,ncup)=J1_coupling(3,ncup)
 
      End Subroutine make_coupling_insert


!======================================================================
      Subroutine make_coupling_trap
!======================================================================
!     prepare coupling scheme for state in module conf_jj
!     with outer orbital is going to position ii
!----------------------------------------------------------------------

      Use bsr_conf,  JJ_coupling => J2_coupling   

      Use conf_LS

      Integer :: i,ii
 
      n = ncup + 1
      m = n + n

! ... JJ-coupling:

      ii = iabs(insert)

      JJ_coupling(1,1) = ii
      JJ_coupling(2,1) = n
      JJ_coupling(3,1) = m+1

      JJ_coupling(1,2) = 1
      JJ_coupling(2,2) = 2
      JJ_coupling(3,2) = m+2
	  
      if(ii.eq.1) then
       JJ_coupling(1,2) = m+1
       JJ_coupling(2,2) = 2
       JJ_coupling(3,2) = m+2
      end if

      Do i=2,no-1
       JJ_coupling(1,i+1) = m+i
       JJ_coupling(2,i+1) =   i+1
       JJ_coupling(3,i+1) = m+i+1
      End do

      if(ii.gt.1) JJ_coupling(2,ii) = m+1

      JJ_coupling(3,ncup)=J1_coupling(3,ncup)

! ... moments

      momentS(m+1) = LS(ii,3)
      momentL(m+1) = LS(ii,2)
      Do i=2,no
       momentS(m+i) = LS(i,5)
       momentL(m+i) = LS(i,4)
      End do

      End Subroutine make_coupling_trap


