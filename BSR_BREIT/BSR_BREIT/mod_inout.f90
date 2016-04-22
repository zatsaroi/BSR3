!======================================================================
      MODULE inout_br
!======================================================================
!
!     Containes the names and units for input/output files
!
!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number
!
!----------------------------------------------------------------------

      Implicit none
      Save

! ... runing information:

      Integer :: pri=66;  Character(40) :: AF_pri = 'bsr_breit.log'

! ... c-file:

      Integer :: nuc=1;  Character(40) :: AF_c = 'cfg.inp'
                            Character(40) :: BF_c = 'cfg.###'

! ... data bank:

      Integer :: nub=2;  Character(40) :: AF_b = 'int_bnk'
                            Character(40) :: BF_b = 'int_bnk.###'

      Logical(1) :: new     ! pointer on the previous calculation  

! ... new results if any:

      Integer :: nur=3;  Character(40) :: AF_r = 'int_res'
                            Character(40) :: BF_r = 'int_res.###'

      Logical(1) :: icalc   ! pointer for need of new calculations

! ... scratch files:

      Integer :: nui= 4  ! intermediate results

      Integer :: nus=11  ! for reallocations   
      Integer :: nud=12  ! for det. expansions
      Integer :: nua=13  ! for accumulation of data

! ... maximum record length:  

      Integer, parameter :: mrecl = 100000

!     we introduce this parameter because different platforms have 
!     different default record lengths, and it may cause troubles when 
!     change computers or compilers

! ... range of partial waves:

      Integer :: klsp1 = 1 
      Integer :: klsp2 = 1 

      End MODULE inout_br



!======================================================================
      Subroutine Open_br(nu,klsp)
!======================================================================
!
!     open (closed) files in the breit_bsr program
!
!----------------------------------------------------------------------

      USE inout_br

      Implicit none
      Integer, Intent(in) :: nu,klsp
      Integer :: i,iarg
      Character(40) :: AF,BF
      Character(3) :: ALSP
      Logical :: EX
      Integer, External :: IARGC

      Call Anumber(klsp,3,ALSP)

       Select case(nu)

       Case(1)               ! c-file

        AF = AF_c
        if(klsp.gt.0) then
         i=Index(BF_c,'.'); AF=BF_c(1:i)//ALSP; BF_c=AF
        else
         iarg = IARGC()
         if(iarg.gt.0) then
          Call GETARG(1,BF); i=INDEX(BF,'=')
          if(i.eq.0) then 
           i=INDEX(BF,'.')-1; if(i.lt.0) i=LEN_TRIM(BF)
           AF=BF(1:i)//'.c'; AF_b = BF(1:i)//'.bnk'
          end if
         end if
        end if  

        Call Check_file(AF)
        Open(nu,file=AF,status='OLD')

       Case(2)               ! bnk-file

        AF = AF_b
        if(klsp.gt.0) then
         i=Index(BF_b,'.'); AF=BF_b(1:i)//ALSP; BF_b=AF
        end if  
        Inquire (FILE=AF, exist=EX)

        new=.TRUE.
        if(EX) then
         new = .FALSE.
         Open(nub,file=AF,form='UNFORMATTED',STATUS='OLD')
        end if   

        Case(3)

         AF = AF_r
         if(klsp.gt.0) then
          i=Index(BF_r,'.'); AF=BF_r(1:i)//ALSP; BF_r=AF
         end if  
         Open(nur,file=AF,form='UNFORMATTED')

        Case(66)

         AF = AF_pri;  Open(nu,file=AF)

        Case(4,11,12,13)

         Open(nu,form='UNFORMATTED',status='SCRATCH')

        Case default
        
         Stop ' open_br: nu is out of list '

       End select

       End Subroutine Open_br

