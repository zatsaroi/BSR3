!======================================================================
!     Program BSR_HF
!======================================================================
!                   C O P Y R I G H T -- 2020
!     Written by:   Oleg Zatsarinny 
!     email:        oleg_zoi@yahoo.com
!----------------------------------------------------------------------
!     This program computes the radial functions for simple HF 
!     problem in B-spline basis. In addition to standard one-electron
!     case, the program allows also simultaineous calculation of
!     several configurations in different optimization schemes.
!     Another possibility - calculation whole Rydberg series in
!     frozen-core approximation.   
!
!     For short instructions, run dbsr_hf ? and look in name.inp.
!----------------------------------------------------------------------
      Use bsr_hf

      Implicit none
      Real(8) :: t1,t2,t3,t4

      Call CPU_time(t1)

! ... atomic paramters:

      Call Get_case

! ... set up the energy expression:

      Call Def_energy_coef

! ... B-spline parameters:

      Call Get_spline_param

! ... computing:

      Call CPU_time(t3)

      Call Get_estimates

      Call CPU_time(t4)
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'Get_estimates:',t4-t3,'  sec' 

      Call CPU_time(t3)

      Call Solve_HF

      Call CPU_time(t4)
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'Solve_HF:',t4-t3,'  sec' 

! ... output of results:

      Call CPU_time(t3)
      Call Write_bsw
      Call Write_nl
      Call Plot_bsw
      Call Summry
      Call CPU_time(t4)
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'Summry:',t4-t3,'  sec' 

      Call CPU_time(t2)
      write(scr,'(/a,f10.2,a)') 'time:',t2-t1,'  sec' 
      write(log,'(/a,f10.2,a)') 'time:',t2-t1,'  sec' 

      End ! Program BSR_HF

