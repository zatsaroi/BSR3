!----------------------------------------------------------------------
      Module bsr_mat
!----------------------------------`------------------------------------
!     main parameters for the BSR_MAT program
!----------------------------------------------------------------------
      Use target
      Use channel

      Use spline_param
      Use spline_atomic
      Use spline_grid
      Use spline_galerkin
      Use spline_orbitals

      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma=40    ! limit for file-name length
      Character(ma) :: AF,BF     

      Integer :: nut = 10; Character(ma) :: AF_tar = 'target'
      Integer :: nup = 11; Character(ma) :: AF_par = 'bsr_par'
      Integer :: nuw = 14; Character(ma) :: AF_bsw = 'target.bsw'
      Integer :: nuo = 17; Character(ma) :: AF_orb = 'target_orb'

! ... relativistic corrections:

      Integer :: mlso  =   5         !  max. l for so-interaction
      Real(8) :: eps_soo  = 1.0D-02  !  tolerance for two-electron
                                     !  relativistic integrals
! ... dimension limits:

      Integer :: nb   =   2000   ! number of blocks in c_data       
      Integer :: mb   =   2000   ! size of blocks
      Integer :: kb   =    100   ! max.nb for the given type of integrals

      Integer,parameter :: mtypes = 1000
      Integer :: ntypes = 0
      Integer :: ktypes (mtypes)   
	                               
! ... tolerence parameters:

      Real(8) :: eps_det  = 1.0D-10   !  tolerance for determinants
      Real(8) :: eps_ovl  = 1.0D-10   !  tolerance for overlaps

! ... packing basis for orbitals in the bsr_mat (as i*ibo+j):

      Integer, parameter :: ibo = 2**15

      Integer, allocatable :: ip_channel(:)

! ... structure of data:

      Integer, parameter :: ibi = 2**15 ! packing basis for orbitals

      Integer :: ipert_ch = 1
      Integer :: check_target = 1
      Integer :: iitar = 0

      End Module bsr_mat


!======================================================================
      Integer Function Ifind_type(icase,kpol,itype)
!======================================================================
      Use bsr_mat
      
      Integer :: icase,kpol,itype, k

      k = 100 * (100*(kpol+1) + icase) + itype
      
      m = 0
      Do i = 1,ntypes
       if(k.ne.ktypes(i)) Cycle
       m = i; Exit
      End do

      Ifind_type = m
      if(m.gt.0) Return

      if(ntypes.gt.mtypes) Stop 'ntypes > mtypes'
      ntypes=ntypes+1
      ktypes(ntypes) = k
      Ifind_type = ntypes

      End Function Ifind_type


!======================================================================
      Subroutine Decode_type (jtype,icase,kpol,itype)
!======================================================================
      Use bsr_mat
      Integer :: jtype,icase,kpol,itype,k

      k = ktypes(jtype)
      kpol = k / 10000
      k = k - kpol*10000
      kpol = kpol - 1
      icase = k /100
      itype = k - icase*100

      End  Subroutine Decode_type      


   
 

      


     