!----------------------------------------------------------------------
      Module bsr_recoup
!----------------------------------`------------------------------------
!     main parameters for the BSR_MAT program
!----------------------------------------------------------------------
      Use target
      Use channel

      Implicit none

! ... main files: 

      Integer, parameter :: ma=80;  Character(ma) :: AF

      Integer :: pri = 66; Character(ma) :: AF_p     = 'recoup.nnn'
      Integer :: nut = 10; Character(ma) :: AF_tar   = 'target_JK'
      Integer :: mut = 11; Character(ma) :: BF_tar   = 'target_LS'
      Integer :: nue = 12; Character(ma) :: AF_expn  = 'target_LS_expn'
      Integer :: nui = 13; Character(ma) :: AF_mat   = 'bsr_mat_JK.nnn'

      Integer :: kui = 20; Character(ma) :: BF_mat   = 'bsr_mat_LS.nnn'

! ... range of partial waves:

      Integer :: klsp=0, klsp1=0, klsp2=0

      Integer :: pri_ac=0
      Real(8) :: eps_acf  = 1.0D-5   !  tolerance for asympt. coeff.s
      Character(200) :: line

! ... LSJ target expansions: 

      Integer, allocatable :: ip_expn(:), it_expn(:)
      Real(8), allocatable :: c_expn(:)

! ... Hamiltonian matrix information:

      Integer :: mk =  7   ! max. multipole index
      Integer :: km =  7   ! max. multipole index

      Integer :: ns =  0   !  block size
      Integer :: maxblk =  0   !  number of blocks

      Real(8), allocatable :: hcc(:,:,:),hcb(:,:),hbb(:)
      Real(8), allocatable :: ACF(:,:,:), BCF(:,:,:)
      Real(8), allocatable :: x(:,:)
      Integer, allocatable :: icc(:,:), icb(:,:), ibb(:,:),imycase(:,:)

      Integer :: iicc, iicb, iibb  !  dimensions

! ... MPI insert:

      Integer, allocatable :: my_channel(:)
      Integer :: nprocs=1, myid=0, ierr=0
      Integer, allocatable :: ip_proc(:)
      Integer :: debug=0

! ... Calculations in parts:     ???

      Integer :: ij_block = 0
      Integer :: I1_channel, I2_channel, J1_channel, J2_channel 
      Integer :: I1_chan=0, I2_chan=0, J1_chan=0, J2_chan=0 

      End Module bsr_recoup

     