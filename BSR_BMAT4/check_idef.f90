!======================================================================
      Integer Function Check_idef(icase,jtype,ich,jch)
!======================================================================
      Use bsr_mat

      Implicit none
      Integer, intent(in) :: icase,jtype,ich,jch

      Check_idef = 1

      Select case (icase)

       Case(1,2,3,4,5,8,9,10)  ! Call I_data(jtype,jpol) 

        if(jtype.ne.4) Return 
        if(coupling.eq.'LS'.and.icase.gt.7) Return
        if(iitar.gt.0.and.ich.ne.jch) Return
        check_idef = 0

       Case(6)                 ! Call L_data(jtype,jpol)  

        if(jtype.eq.4) then
          if(iitar.gt.0.and.ich.ne.jch) Return
          check_idef = 0
          Return 
        end if
        if(jtype.ge.13.and.jtype.le.16) then
          if(ich.eq.jch.or.iitar.gt.1) Return
          check_idef = 0
          Return 
        end if

       Case(7)                 ! Call Z_data(jtype,jpol)  

        if(jtype.ne.4) Return 
        if(coupling.eq.'LS') Return
        check_idef = 0

       Case(11)                ! Call O_data(jtype,jpol)

        if(jtype.ne.4) Return 
        if(iitar.gt.1.and.ich.ne.jch) Return
        check_idef = 0

      End Select

      End Function Check_idef



