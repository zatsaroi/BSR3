!======================================================================
      MODULE bsr_hf
!======================================================================
!     main parameters for the DBSR_HF program
!----------------------------------------------------------------------
      Use zconst
      Use spline_atomic
      Use spline_param
      Use spline_galerkin
      Use spline_grid
      Use spline_hl
      Use spline_integrals
      Use spline_moments
      Use spline_slater

      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma = 80   ! limit for file-name length
      Character(ma) :: AF     

      Integer :: inp = 5;   Character(ma) :: AF_dat  = 'name.inp'
                            Character(ma) :: BF_dat  = '.inp'
      Integer :: log = 3;   Character(ma) :: AF_log  = 'name.log'                            
                            Character(ma) :: BF_log  = '.log'
      Integer :: scr = 0
      Integer :: nuw = 11;  Character(ma) :: AF_inp  = 'name.bsw'
                            Character(ma) :: BF_inp  = '.bsw'
                            Character(ma) :: AF_out  = 'name.bsw'
                            Character(ma) :: BF_out  = '.bsw'
                            Character(ma) :: AF_w    = 'name.w'
                            Character(ma) :: BF_w    = '.w'
                            Character(ma) :: AF_nl   = 'name.nl'
                            Character(ma) :: BF_nl   = '.nl'
      Integer :: nup = 12;  Character(ma) :: AF_plt  = 'name.plot'
                            Character(ma) :: BF_plt  = '.plot'
      Integer :: nuc = 13;  Character(ma) :: AF_conf = 'name.conf'
                            Character(ma) :: BF_conf = '.conf'
      Integer :: nus = 14;  Character(ma) :: AF_cfg  = 'name.cfg'
                            Character(ma) :: BF_cfg  = '.cfg'
      Integer :: nuk = 15;  Character(ma) :: CF_grid = 'name.knot'
                            Character(ma) :: BF_grid = '.knot'
      Integer :: nua = 21;  ! scratch file

! ... name of case:

      Character(ma) :: name = ' '
      Character(ma) :: knot = 'knot.dat'

! ... atomic parameters: 

!      Real(8) :: z = 0.d0     ! in spline_atomic
      Integer :: an = 0, ai = 0
      Real(8) :: atw = 0.d0
      Real(8) :: Etotal, E1body, E2body, Ecore

      Character(2)   :: atom = ' '
      Character(2)   :: ion  = ' '
      Character(160) :: configuration = ' ', conf_AV = ' ', conf_LS = ' '
      Character(2)   :: term = 'AV'
      Character(80)  :: anit = 'all'

      Integer :: nelc  = 0
      Integer :: nconf = 0
      Integer :: ncfg  = 0
      Integer :: eal = 5
      Real(8), allocatable :: weight(:)
      Integer, allocatable :: iqconf(:,:)

! ... convergence:

      Real(8) :: scf_tol=1.d-11,  scf_diff
      Real(8) :: orb_tol=1.d-08,  orb_diff
      Real(8) :: end_tol=1.d-08
      Real(8) :: eps_ovl=1.d-08
      Integer :: max_it = 75

      Integer :: ac = 0
      Real(8) :: aweight = 0.7
      Real(8) :: bweight = 0.7

! ... debug options:

      Integer :: debug  = 0
      Integer :: rotate = 0
      Character :: meth = 'c'
      
! ... core 

      Integer,parameter :: mcore = 50
      Character(250) :: core=' '
      Integer :: ncore = 0
      Integer :: n_core(mcore), l_core(mcore) 
      Character(4) :: e_core(mcore)

! ... description of 1 conf.w.function:  

      Integer, parameter :: msh = 25 ! max. number of shells behind core
      Integer :: no
      Integer, dimension(msh) :: nn,ln,iq,in
      Integer :: LS(msh,5)
      Integer :: Ltotal=0, Stotal=0

! ... Storing configuration as character strings:

      Character(200) :: CONFIG, COUPLE

      Integer :: nwf = 0
      Integer :: nit = 0
      Integer :: ilzero = 0
      Integer :: ibzero = 2
      Integer :: kmin = 0
      Integer :: kmax = 0

      Real(8) :: eps_c = 1.d-6
      Integer :: ibi = 2**16
      Integer :: mdiag = 0

! ... solutions for Rydberg series:

      Integer :: out_nl = 0
      Integer :: nsol_nl = 0
      Real(8), allocatable :: p_nl(:,:), e_nl(:)

      Integer :: out_w = 0     ! output in the GRASP format
      Integer :: out_plot = 0  ! output in table form

! ... frequently called functions: (instead interface)  

      Integer, external :: Icheck_file

! ... debuging time:

      Real(8) :: time_hf_eiv=0.d0, time_hf_matrix=0.d0, &
                 time_hf_matrix_breit=0.d0, time_update_int=0.d0
      Real(8) :: au_cm, au_eV

      End Module bsr_hf

