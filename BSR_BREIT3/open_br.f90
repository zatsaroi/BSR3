!======================================================================
      Subroutine Open_br(nu,klsp)
!======================================================================
!     open (close) files in the breit_bsr program
!----------------------------------------------------------------------
      Use bsr_breit

      Implicit none
      Integer, Intent(in) :: nu,klsp
      Integer :: i,iarg
      Character(ma) :: AF,BF
      Character(3) :: ALSP
      Integer, external :: IARGC, Icheck_file

      write(ALSP,'(i3.3)') klsp 

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
           i=INDEX(BF,'.',BACK=.TRUE.)-1; if(i.lt.0) i=LEN_TRIM(BF)
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

        new=1
        if(Icheck_file(AF).eq.1) then
         new = 0
         Open(nub,file=AF,form='UNFORMATTED',STATUS='OLD')
        end if   

        Case(3)

         AF = AF_r
         if(klsp.gt.0) then
          i=Index(BF_r,'.'); AF=BF_r(1:i)//ALSP; BF_r=AF
         end if  
         Open(nur,file=AF,form='UNFORMATTED')

        Case(66)

         Open(nu,file=AF_p)

        Case(4,11,12,13)

         Open(nu,form='UNFORMATTED',status='SCRATCH')

        Case default
        
         write(*,*) 'nu,klsp =', nu,klsp
         Stop ' open_br: nu is out of list '

       End select

       End Subroutine Open_br

