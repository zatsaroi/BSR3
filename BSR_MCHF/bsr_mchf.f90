!======================================================================
!     PROGRAM BSR_MCHF
!======================================================================
!     MULTICONFIGURATION HARTREE-FOCK PROGRAM
!----------------------------------------------------------------------
!                   C O P Y R I G H T -- 2020
!     Written by:   Oleg Zatsarinny with help of Charlotte Froese Fischer
!     email:        oleg_zoi@yahoo.com
!----------------------------------------------------------------------
!     This program is a part of the BSR complex and computes 
!     radial one-electron functions in B-spline basis for the 
!     multi-configuration Hartree-Fock problem. 
!!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Real(8) :: t1,t2
      Integer :: ILT,IST,IPT
      Character, external :: AL
      Real(8), external :: Ecore_hf

      Call CPU_time(t1)

! ... read main parameters:

      Call Get_case

! ... read configurations from c-file:

      Call Read_conf

! ... B-spline parameters:

      write(log,'(/a,T20,a)') 'atom:',atom
      Call Def_term_BSR(nuc,ILT,IST,IPT)
      write(term,'(i1,a1)') IST,AL(ILT,2)
      write(log,'(/a,T20,a)') 'term:',term
      Call Get_spline_param

! ... prepare radial-function arrays:

      Call Def_orbitals

      Call Get_estimates

! ... read angular coefficients: 

      Call Read_ang_coef

      Ecore = Ecore_hf()

! ... Computing:

      if(all.eq.0)   Call scf_mchf
!      if(all.eq.1)   Call scf_mchf_all

! ... output the results:

      Call Write_bsw

      Call Summry

      Call CPU_time(t2)
      write(scr,'(/a,f10.2,a)') 'time:',t2-t1,' sec'
      write(log,'(/a,f10.2,a)') 'time:',t2-t1,' sec' 

      write(log,'(/80(''-'')/)')

!      Call Debug_time

      End ! Program BSR_MCHF

