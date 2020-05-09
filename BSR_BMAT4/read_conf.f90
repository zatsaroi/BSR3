!======================================================================
      Subroutine Read_conf(nuc)
!======================================================================
!     Read the configuration list from c-file and alloc auxiliary arrays 
!----------------------------------------------------------------------
      Use bsr_mat, only:  eps_soo

      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS
      Use orb_LS  

      Implicit none
      Integer, intent(in) :: nuc
 
! ... remove data from previous calculations:

      Call Alloc_cfg_LS (0)
      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)
      Call Alloc_orb_LS(0)
      Call R_closed(nuc)

! ... define configurations in c-file:

      Call Add_conf_LS(nuc,0)

! ... sorting states according terms (define IT_state1 and IT_state2):

      Call Sort_states  

! ... sorting the terms according to configurations,
! ... (define IC_term1 and IC_term2):

      Call Sort_terms  

! ... define IT_oper, IT_need, IC_need:

      Call Prepare_it_oper(eps_soo)

      End Subroutine Read_conf



