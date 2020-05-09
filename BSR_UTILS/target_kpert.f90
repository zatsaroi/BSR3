!=====================================================================
!     UTILITY  target_kpert
!=====================================================================
!
!     target_LS -> target in LSJ scheme with the kpert list 
!
!---------------------------------------------------------------------
      Use target; Use channels

      Implicit real(8) (A-H,O-Z)

      Character(80) ::  AF

      Integer :: nu1=1;  Character(80) :: AF_LS   = 'target_LS'
      Integer :: nu2=2;  Character(80) :: AF_LSJ  = 'target'
      Integer :: nu3=3;  Character(80) :: AF_pert = 'thresholds_pert'
      Integer :: nu4=4;  Character(80) :: AF_thresholds = 'thresholds'

      Integer, allocatable :: index_kpert(:,:)

!----------------------------------------------------------------------
! ... original LS target file:

      Call Check_file(AF_LS)
      Open(nu1,file=AF_LS,status='OLD')
      Call R_target  (nu1)
      Call R_channels(nu1)
      rewind(nu1); read(nu1,'(a)') title
      Close(nu1)

      Z = nz
      Call Conv_au (Z,0.d0,au_cm,au_eV,0)

!----------------------------------------------------------------------
! ... output target file:

      open(nu2,file=AF_LSJ) 

      write(nu2,'(a)') 'target file for pertuber calculations' 
      write(nu2,'(80(''-''))')
      write(nu2,'(a,a2,4x,a)') &
                'coupling = ',coupling, ' !   coupling scheme'
      write(nu2,'(a,i4,5x,a)') &
                'nz     = ',nz,   '!   nuclear charge' 
      write(nu2,'(a,i4,5x,a)') &
                'nelc   = ',nelc-1, '!   number of electrons'

      write(nu2,'(80(''-''))')
      write(nu2,'(a,i4,5x,a)') &
                'ntarg = ',1,' !   number of target states'
      write(nu2,'(80(''-''))')
      write(nu2,'(a)') 'any ion state'
      write(nu2,'(80(''-''))')
      
!----------------------------------------------------------------------
! ... find all possible J-values:

      JJ_min = 1000; JJ_max = 0
      Do it = 1,ntarg 
       JJ = abs(2*ltarg(it)+1-istarg(it))
       if(JJ.lt.JJ_min) JJ_min=JJ
       JJ = 2*ltarg(it)+istarg(it)-1
       if(JJ.gt.JJ_max) JJ_max=JJ
      End do

      nlsp1 = JJ_max-JJ_min + 2

      write(nu2,'(a,i4,5x,a)') &
           'nlsp  = ',nlsp1,' !   number of partial waves' 
      write(nu2,'(80(''-''))')
      i = 0
      Do JJ = JJ_min,JJ_max,2
       i = i+1; write(nu2,'(i3,3i4)')  i,JJ,0,-1
       i = i+1; write(nu2,'(i3,3i4)')  i,JJ,0,+1
      End do
      write(nu2,'(80(''-''))')

!----------------------------------------------------------------------
! ... find all kpert:

      open(nu3,file=AF_pert) 

      kpert = 0; ilsp = 0
      Do JJ = JJ_min, JJ_max, 2
       Do ip = -1,1,2
        ilsp = ilsp + 1
        Do it=1,ntarg; if(ip.ne.iptarg(it)) Cycle
         if(ITRI(JJ+1,2*ltarg(it)+1,istarg(it)).eq.0) Cycle
         kpert = kpert + 1
         write(nu3,'(i5,3x,a)') ilsp, trim(BFT(it)),it
        End do
       End do
      End do

      write(nu2,'(a,i4,5x,a)') 'kpert = ',kpert,' !    number of perturbers'
      write(nu2,'(80(''-''))')
                                                                                     
      Allocate(index_kpert(nlsp1,ntarg)); index_kpert=0

      rewind(nu3)
      Do i = 1,kpert

       read(nu3,*) ilsp, AF,it

       Do j=1,ntarg
        if(index_kpert(ilsp,j).ne.0) Cycle
        index_kpert(ilsp,j)=it
        Exit
       End do

       write(nu2,'(i5,3x,a)') ilsp, trim(AF)

      End do

      write(nu2,'(80(''-''))')

!----------------------------------------------------------------------
! ... create thrsholds_pert

      rewind(nu3)

      Do it = 1,ntarg
       write(nu3,'(a,T30,f16.8,f10.3)') BFT(it),etarg(it),0.0
      End do 
      write(nu3,'(/a,i5)') 'units = eV    -  units for corrections, column 3'
      write(nu3,'(/a,2i5)') 'nperts = ', ntarg, nlsp1
      Do i=1,nlsp1
       write(nu3,'(10i7)') index_kpert(i,:)
      End do

      open(nu4,file=AF_thresholds) 
      write(nu4,'(f16.8)') etarg(1)
      
      End ! program
