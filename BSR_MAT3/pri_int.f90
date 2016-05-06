!=======================================================================
      Subroutine pri_int (nu,icase,k,i1,i2,i3,i4,S)
!=======================================================================
!     print the integral for debug purposes
!-----------------------------------------------------------------------
      Use spline_orbitals,  only: EBS

      Implicit none
      Integer, intent(in) :: nu,icase,k,i1,i2,i3,i4
      Real(8), intent(in) :: S
      Character AINT(11)/'F','G','T','S','R','L','Z','N','V','N','O'/

      if(icase.eq.6.or.icase.eq.7) then
       write(nu,'(a,a,a,a,a,a,D14.7,f12.1)') &        
        AINT(icase),'(',EBS(i1),',',EBS(i2),')=',S,S*219474   
      else
       write(nu,'(a,i2,a,a,a,a,a,a,a,a,a,D14.7,f12.1)') &
        AINT(icase),k,'(',EBS(i1),',',EBS(i2),';', &    
        EBS(i3),',',EBS(i4),')=',S,S*219474  
      end if

      End Subroutine pri_int