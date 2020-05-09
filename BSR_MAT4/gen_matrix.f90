!======================================================================
      Subroutine Gen_matrix(jtype,jpol)
!======================================================================
!     generate the interaction matrix for current "itype" and "kpol"
!----------------------------------------------------------------------
      Use bsr_mat
      Use c_data, nc => ncdata

      Implicit none
      Integer, intent(in) :: jtype,jpol

! ... generate matrix:

      Select case (icase)
       Case(1,2,3,4,5,8,9,10); Call I_data(jtype,jpol) 
       Case(6);                Call L_data(jtype,jpol)  
       Case(7);                Call Z_data(jtype,jpol)  
       Case(11);               Call O_data(jtype,jpol)
      End Select

      End Subroutine Gen_matrix
