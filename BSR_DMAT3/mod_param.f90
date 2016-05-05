!======================================================================
      Module  bsr_dmat
!======================================================================
!     containes some common parameters, variables and arrays
!     specific for this program
!----------------------------------------------------------------------
      Implicit none

! ... file units:

      Integer :: in1 = 1   ! 1-st c-file
      Integer :: in2 = 2   ! 2-nd c-file
      Integer :: inb1 = 11 ! 1-st set of expansions
      Integer :: inb2 = 12 ! 2-nd set of expansions

      Integer, parameter :: ma = 40

      Integer :: nub = 3;  Character(ma) :: AF_bnk = 'mult_bnk'
      Integer :: nud = 4;  Character(ma) :: AF_d   = 'd.nnn'
      Integer :: pri =66;  Character(ma) :: AF_log = 'bsr_dmat.log'
      Integer :: nuw = 7;  Character(ma) :: AF_bsw = 'target.bsw'
      Integer :: nut = 8;  Character(ma) :: AF_tar = 'target'
      Integer :: nuo = 9;  Character(ma) :: AF_rsol= 'rsol.nnn'
      Integer :: nuv = 10; Character(ma) :: AF_v   = 'dv.nnn'
      Integer :: nur = 13; Character(ma) :: AF_res = 'zf_res'

      Character(ma) :: name1,name2,AF

! ... transition type:

      Integer :: kpol = 1
      Character(1) :: ktype = 'E'

! ... LSJ-mode:

      Integer :: jmode=0, jot1=0, jot2=0, jout=0
      Character(1) :: ctype1='c', ctype2='c'

! ... number of configurations, channels, perturters:

      Integer :: ilsp1, nch1, ncp1, npert1, ipert1
      Integer :: ilsp2, nch2, ncp2, npert2, ipert2
      Character(3) :: ALS1, ALS2

! ... number of solutions:

      Integer :: nstate1, mstate1=0, istate1=0
      Integer :: nstate2, mstate2=0, istate2=0

! ... expansion coefficients:

      Real(8), allocatable :: C1(:), C2(:)

! ... terms:

      Integer :: ILT1,IST1,IPT1
      Integer :: ILT2,IST2,IPT2

! ... tolerences for expansion coefficients, oneelctron overlaps,
! ... and tota; determinant factors:

      Real(8), parameter :: Eps_c = 1.d-10      
      Real(8), parameter :: Eps_ovl = 1.d-10      
      Real(8), parameter :: Eps_ndet = 1.0D-10

! ... atomic constants:

      Real(8) :: time_au = 2.4189d-17
      Real(8) :: c_au = 137.03599976d0
      Real(8) :: au_eV = 27.2113834d0
      Real(8) :: au_cm = 219471.62d0
      Real(8) :: ge = 2.0023193043622d0

! ... packing number for overlaps:

      Integer :: ibo = 2**15    ! 32768

! ... additional output for polarizabilities:

      Integer :: ialfa = 0

! ... output gf or f-values:

      Character(1) :: GF = 'g'

! ... output gf or f-values:

      Integer :: debug = 0

! ... transition matrixes (1:ns,1:ks+ks-1) in B-spline basis:

      Real(8), allocatable :: bbbs(:,:)    !  <.|1   |.>
      Real(8), allocatable :: rrbs(:,:)    !  <.|r^k |.>
      Real(8), allocatable :: ddbs(:,:)    !  <.|d/dr|.>
      Real(8), allocatable :: ttbs(:,:)    !  <.|1/r |.>

! ... one-electron transition vectors (1:ns,1:nbf):

      Real(8), allocatable :: rbs (:,:)    !  <.|r^k |p>
      Real(8), allocatable :: dbsr(:,:)    !  <.|d/dr|p>
      Real(8), allocatable :: dbsl(:,:)    !  <p|d/dr|.>
      Real(8), allocatable :: tbs (:,:)    !  <.|1/r |p>

! ... one-electron integrals (1:nbf,1:nbf):

      Real(8), allocatable :: dipL(:,:)    !  <p|d_L|q>
      Real(8), allocatable :: dipV(:,:)    !  <p|d_V|q>
      
      End Module  bsr_dmat

