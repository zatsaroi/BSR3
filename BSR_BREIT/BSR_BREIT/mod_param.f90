!======================================================================
      MODULE  param_br
!======================================================================
!
!     main parameters used in the program
!
!     REMARKS:
!
! 1.  Parking basis for orbitals 'jb' should be > nsh,
!     where nsh - max. number of shells, as define in module CONFIGS.
!     We use spectroscopic description for configurations in so-called
!     c-files (MCHF complex) with nsh=8, so it restricts the possible
!     configurations by 8 open shells above some closed-shells core. 
!     It seems not to imply strong restrictions in real calculations
!     except inner-shell excitation in many-electron atoms.   
!
! 2.  The parking basis for det.overlaps and det.factors (ibd,ibf)
!     seems not to restrict any possible calculations.
!
! 3.  The initial dimensions for some array are found to be enough for
!     simple calculations, do not take much memory and the corresponding
!     arrays will be reallocated dynamically if it is necessary.
!
! 4.  There is hidden parameter 'ibc=1^15~32000' from module CONFIGS
!     which determines here the max.number of differerent symmetries.
!     (in int_bnk we use parking it*ibc+jt). The number of symmetries
!     usually is much less then the number of configurations, so it 
!     it does not seem as a restriction.
!
! 5.  Due to mode of parking the integrals, the max. multipole index
!     should be less then jb^4-1 = 9999, that is not a restriction.
!     The default value of mk can be change from the input.
! 
! 6.  Parameters jbl,jbm is used for encoding (parking) of two-electron
!     integrals in uncouple lms-reprizantation. It restricts orbitals 
!     by l < 64. To reduce this restriction, we can use another parking
!     with 4 indentifiers, instead 3 as now.
!
!----------------------------------------------------------------------

      Implicit none
      Save

! ... tolerence for coefficients:

      Real(8) :: Eps_c = 1.d-7      

! ... parking basis for orbitals in the data bank:

      Integer, parameter :: jb=10, jb4=jb**4, jb8=jb**8

! ... maximum multipole index:

      Integer :: mk = 7 

! ... other parameters:

      Integer, parameter :: isd = 20000  ! initial dimension for det.overlaps
      Integer, parameter :: jsd = 3      ! avarage size of overlap det.
      Integer, parameter :: ibd = 2**15  ! parking basis for det.overlaps

      Integer, parameter :: isf = 200000 ! initial dimension for det.factors
      Integer, parameter :: jsf = 5      ! avarage number of dets. in def.
      Integer, parameter :: ibf = 16     ! parking basis for det.factors

! ... initial (supposed) number of coef.s:

      Integer, parameter :: iszoef = 2000    ! in module ZOEF
      Integer, parameter :: iscoef = 25000   ! in module COEF
      Integer, parameter :: isboef = 50000   ! in module BOEF
      Integer, parameter :: isblk  = 5000    ! in module BOEF   

      Integer, Parameter :: jbl  = 2**7      ! in Check_BOEF
      Integer, Parameter :: jbm  = 2**6      ! in Check_BOEF

! ... switch for Vk -> V'k in soo interaction:

      Integer :: is_soo = 0

      End MODULE  param_br

