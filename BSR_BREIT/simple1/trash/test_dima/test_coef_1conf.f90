!=====================================================================
!     PROGRAM   test_coef_ee_1conf                      
!
!               C O P Y R I G H T -- 2013
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!    it is a debug program to check subroutine     "coef_ee_1conf" 
!    to generate angular coefficient for Breit_Pauli calculations 
!    for one-configuration in case of orthogonal one-electron 
!    radial functions
!
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:  name  or   cfg=...  tab=...  
!    
!----------------------------------------------------------------------
!
!    INPUT FILE:     AF_cfg    (default  - cfg.inp) 
!    OUTPUT FILES:   AF_tab    (default  - coef.tab)
!    
!---------------------------------------------------------------------     

      Implicit none 

      Integer :: nuc=1; Character(40) :: AF_cfg = 'cfg.inp'
      Integer :: out=2; Character(40) :: AF_tab = 'coef.tab'
      Character(40) :: name = ' '

      Integer, parameter :: msh = 8
      Integer :: no, nn(msh),ln(msh),iq(msh),kn(msh), LS(msh,5)

      Character(8*msh) :: CONFIG, COUPLE
      Real(8) :: C, eps_C = 1.d-7
      Integer :: i,j, k, ic, ncfg
      Integer, external :: Idef_ncfg 

! ... result coefficients:

      Real(8), allocatable :: coefs(:,:,:)
      Integer :: kmax
!----------------------------------------------------------------------
      Call Read_name(name)
      if(len_trim(name).ne.0) then
       AF_cfg=trim(name)//'.c'
       AF_tab=trim(name)//'.tab'
      else
       Call Read_aarg('cfg',AF_cfg)
       Call Read_aarg('tab',AF_tab)
      end if

      Call Check_file(AF_cfg)
      Open(nuc,file=AF_cfg)
      ncfg = Idef_ncfg(nuc)
      if(ncfg.eq.0) Stop 'ncfg=0: nothing to do'
      Open(out,file=AF_tab)

      rewind(nuc)
      Do ic = 1,ncfg

      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(5:5).ne.'(') Cycle
       read(nuc,'(a)') COUPLE
       Exit
      End do

      Call Decode_c (CONFIG,COUPLE,msh,no,nn,ln,iq,kn,LS)

      kmax = 2*maxval(ln(1:no))
      if(allocated(coefs)) Deallocate(coefs)
      Allocate(coefs(no,no,0:kmax))

      Call coef_ee_1conf(msh,no,ln,iq,LS,kmax,coefs)

      write(out,'(a)') trim(CONFIG) 
      write(out,'(a)') trim(COUPLE)
      write(out,*)

      Do i=1,no
       Do j=i,no
        Do k=0,kmax
		       C = coefs(i,j,k); if(abs(C).lt.eps_c) Cycle
         write(out,'(a,i2,a,i2,a,i2,a,f10.5)') &
		            'F',k,'(',i,',',j,')=',C        
        End do
        if(i.eq.j) Cycle
        Do k=0,kmax
		       C = coefs(j,i,k); if(abs(C).lt.eps_c) Cycle
         write(out,'(a,i2,a,i2,a,i2,a,f10.5)') &
		            'G',k,'(',i,',',j,')=',C        
        End do
       End do; End do
       write(out,*)

      End do  ! over ic

      END ! Program ...


!======================================================================
       Subroutine Decode_c(CONFIG,COUPLE,msh,no,nn,ln,iq,kn,LS)
!======================================================================
!      decode the config. from c-file format to ZOI format
!----------------------------------------------------------------------
       Implicit none
       Character(*), intent(in) :: CONFIG,COUPLE
       Integer, intent(in) :: msh
       Integer :: no,nn(msh),ln(msh),iq(msh),kn(msh),LS(msh,5)
       Integer :: ii, k,i,j 
       Integer, external :: LA

       no=0; ii=LEN_TRIM(CONFIG); ii=ii/8;  k=1; j=2
       Do i=1,ii
        if(CONFIG(k+4:k+4).ne.'(') Exit
        no=no+1
        Call EL4_nlk(CONFIG(k:k+3),nn(i),ln(i),kn(i))
        read(CONFIG(k+5:k+6),'(i2)') iq(i)
        k=k+8
        read(COUPLE(j:j),'(i1)') LS(i,3)
        LS(i,2)=2*LA(COUPLE(j+1:j+1))+1
        read(COUPLE(j+2:j+2),'(i1)') LS(i,1)
        j=j+4
       End do

       LS(1,4)=LS(1,2)
       LS(1,5)=LS(1,3)
       Do i=2,no
        read(COUPLE(j:j),'(i1)') LS(i,5)
        LS(i,4)=2*LA(COUPLE(j+1:j+1))+1
        j=j+4
       End do

       End Subroutine Decode_c



