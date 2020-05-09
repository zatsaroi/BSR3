!======================================================================
!     PROGRAM       B S R _ R E C O U P
!
!                   C O P Y R I G H T -- 2017
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!                        
!======================================================================
!     Make recoupling the LS Hamiltionian/Overlap matrixes 
!     to the JK coupling based on the target expansions over LS-states
!======================================================================
!
!     INPUT FILES:
!
!     target_JK      -  description of target states and partial waves                        
!                       in JK-coupling
!     target_LS      -  description of target states and partial waves
!                          in LS-coupling                                
!     bsr_mat_LS.nnn -  LS Hamiltonian matrix for nnn-partion wawe 
!                    
!     target_exp     -  LSJ target expansions over LS-states
!
!     OUTPUT FILES:  
!                    
!     recoup.nnn     -  running information
!
!     bsr_mat_JK.nnn -  LS Hamiltonian matrix for nnn-partion wawe
!                    
!     INPUT ARGUMENT:  see subroutine read_arg.f90
!
!----------------------------------------------------------------------
      Use bsr_recoup

      Implicit real(8) (A-H,O-Z)
      Character(80) :: AS

!      Call bsr_recoup_inf
!---------------------------------------------------------------------- 
! ... read target information:

      Call Check_file(AF_tar)
      Open(nut,file=AF_tar,status='OLD',action='READ')
      Call R_target(nut)

      Call Check_file(BF_tar)
      Open(mut,file=BF_tar,status='OLD',action='READ')
      Call R_target_ion(mut)
      Call R_channels_ion(mut)

! ... read recoupling coefficients:

      Call Check_file(AF_expn)
      Open(nue,file=AF_expn,status='OLD',action='READ')
      Call Read_ipar(nue,'mpert',mpert)
      Allocate(ip_expn(0:ntarg),it_expn(mpert),c_expn(mpert))
      ip_expn = 0
      rewind(nue)
      k = 0
      Do it = 1,ntarg
       read(nue,*)  AS,jt,AS,AS,npert
       if(it.ne.jt) Stop 'it <> jt  in target_LS_expn'
       Do i = 1,npert
        k = k + 1
        read(nue,*) it_expn(k),c_expn(k) 
       End do
       ip_expn(it) = k      
      End do
      if(k.ne.mpert) Stop 'mpert ??? on target_LS_expn'

! ... read arguments:

      Call Read_arg

!----------------------------------------------------------------------
! ... loop over partial waves:

      Do klsp = klsp1,klsp2

       i=len_trim(AF_p)
       write(AF_p(i-2:i),'(i3.3)') klsp
       open(pri,file=AF_p,action='WRITE')

       write(*,  '(/a,i3/)') 'BSR_RECOUP:  klsp =', klsp
       write(pri,'(/a,i3/)') 'BSR_RECOUP:  klsp =', klsp

       Call CPU_TIME(t1);  Call SUB1;  Call CPU_TIME(t2)

       write(pri,'(/a,5x,f8.2,a)') 'Total time:', (t2-t1)/60, ' min'
       write(*  ,'( a,5x,f8.2,a)') 'Total time:', (t2-t1)/60, ' min'

      End do  ! over klsp

      End  ! program bsr_mat
