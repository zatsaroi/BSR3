!======================================================================
      Module bsr_conf
!======================================================================
! ... main parameters and variables for BSR_CONF program
!----------------------------------------------------------------------
      Implicit none

! ... files:

      Integer, parameter :: ma=80
      Integer :: nut = 1;  Character(ma) :: AF_tar ='target'
      Integer :: nup = 2;  Character(ma) :: AF_par ='bsr_par'
      Integer :: pri = 3;  Character(ma) :: AF_log ='bsr_conf.log'
      Integer :: nuo = 7;  Character(ma) :: AF_orb ='target_orb'
      Integer :: nuc = 8;  Character(ma) :: AF      ! c-files
      Integer :: nuw = 9;  Character(ma) :: AF_wfn ='target.bsw'
      
       
! ... default values for parameters:    

      Real(8) :: c_comp = 1.010  ! tolerance for compensation configurations          
      Real(8) :: c_conf = 0.333  ! tolerance for physical configurations        
      Real(8) :: c_pert = 0.500  ! tolerance for physical configurations in perturber        
      Real(8) :: c_norm = 0.100  ! tolerance for target states normalization        
      Real(8) :: c_orth = 0.250  ! compensations < c_orth are ignored        

      Integer :: max_it = -1     ! upper limit for number of target states
      Integer :: max_ll = -1     ! upper limit for small "l"
      Integer :: min_ll = -1     ! lower limit for small "l"
      Integer :: max_LT = -1     ! upper limit for big "L" (2L+1 - value) 
      Integer :: max_ST = -1     ! upper limit for total spin (2S+1 - value)
      Integer :: kort   = -1     ! if kort > 0, program will read additional
                                 ! orth. conditions from cfg.nnn file

! ... values used in different subroutions:

      Integer :: ilsp
      Integer :: ncfg_phys, ncfg_targ, ncfg_sct, ncfg_pert, ncfg_comp
      Integer :: lcfg_phys, lcfg_targ, lcfg_sct, lcfg_pert, lcfg_comp
      Integer :: nwf_targ, nwf_sct, nwf_pert

      Real(8) :: wc_comp = 0.d0
      Integer :: ic_comp = 0
      Integer :: ie_comp = 0
      Integer :: ii_comp = 0
      Integer :: insert  = 0
      Integer :: debug   = 0
      Integer :: igen_conf = 0
      Integer :: iread_targ = 0

! ... coupling scheemes:

      Integer, parameter :: mshells = 25
      Integer, parameter :: mcup=2*mshells   
      Integer, parameter :: mmom=6*mshells   
      Integer :: J1_coupling(3,mcup)
      Integer :: J2_coupling(3,mcup)
      Integer :: momentS(mmom)
      Integer :: momentL(mmom)
      Integer :: ncup, nmom
      Integer :: ILT,IST            ! total term
      Integer :: IAP,ILP,ISP        ! parent term
      Integer :: IL_trap,IS_trap    ! trap_term
      Real(8) :: S_cfp 
      Real(8) :: S_recup 

! ... additional arrays:

      Integer, allocatable :: ip_phys_conf(:)  
      Integer, allocatable :: jj_phys_conf(:)  
      Integer, allocatable :: ic_targ(:)  
      Integer, allocatable :: jc_targ(:)  
      Integer, allocatable :: ic_pert(:)  
      Real(8), allocatable :: WC_pert(:)

! ... channel-delete information if any:

      Character(ma) :: AF_del = 'target_del'
      Integer :: ndel = 0  
      Integer, allocatable :: dlsp(:),dlch(:),dtar(:),djk(:)

      End module bsr_conf


!======================================================================
      Subroutine Read_del
!======================================================================
!     read "chanel-delete" conditions if any
!---------------------------------------------------------------------- 
      Use bsr_conf
      Implicit none
      Integer :: i,j,m,nlsp,nch,lch,iptar,ipconf,jkch,idel,nud
      Integer, external :: Icheck_file, Ifind_position

      if(Icheck_file(AF_del).eq.0) Return
      Call Find_free_unit(nud); Open(nud,file=AF_del)

      Call Read_ipar(nud,'nlsp',nlsp); if(nlsp.le.0) Return
      i = Ifind_position(nud,'channels')
      read(nud,*) AF

      ndel=0
      Do i = 1,nlsp
       read(nud,*) AF
       read(nud,*) AF,AF,AF,AF,nch
       Do j = 1,nch
        read(nud,*) AF,lch,iptar,m,ipconf,jkch
        if(ipconf.eq.0) ndel=ndel+1
       End do
      End do
      if(ndel.eq.0) Return

      Allocate(dlsp(ndel),dlch(ndel),dtar(ndel),djk(ndel))

      i = Ifind_position(nud,'channels')
      read(nud,*) AF
      idel=0
      Do i = 1,nlsp
       read(nud,*) AF
       read(nud,*) AF,AF,AF,AF,nch
       Do j = 1,nch
        read(nud,*) AF,lch,iptar,m,ipconf,jkch
        if(ipconf.gt.0) Cycle
        idel=idel+1
        dlsp(idel)=i
        dlch(idel)=lch
        dtar(idel)=iptar
        djk (idel)=jkch
       End do
      End do

      Close(nud)

      End Subroutine Read_del


!======================================================================
      Integer Function Icheck_del(k1,k2,k3,k4)
!======================================================================
!     check if channel with identifiers k1,k2,k3,k4 is in "delete" list
!----------------------------------------------------------------------
      Use bsr_conf
      Implicit none
      Integer :: i,k1,k2,k3,k4

      Icheck_del = 0
      if(ndel.eq.0) Return
      Do i = 1,ndel
       if(dlsp(i).ne.k1) Cycle
       if(dlch(i).ne.k2) Cycle
       if(dtar(i).ne.k3) Cycle
       if(djk(i) .ne.k4) Cycle
       Icheck_del = 1
       Exit
      End do

      End Function Icheck_del




