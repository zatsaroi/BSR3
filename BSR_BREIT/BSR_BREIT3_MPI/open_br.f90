!======================================================================
      Subroutine Open_br(nu)
!======================================================================
!     open (close) files in the breit_bsr program
!----------------------------------------------------------------------

      USE bsr_breit

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i,iarg
      Character(ma) :: AF,BF
      Character(3) :: ALSP
      Integer, External :: IARGC, Icheck_file

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
           i=INDEX(BF,'.',BACK=.TRUE.)-1; if(i.lt.0) i=LEN_TRIM(BF)
           AF=BF(1:i)//'.c'; AF_b = BF(1:i)//'.bnk'
          end if
         end if
        end if  

        Open(nu,file=AF)

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

        Case(3)              ! results

         AF = AF_r
         if(klsp.gt.0) then
          i=Index(BF_r,'.'); AF=BF_r(1:i)//ALSP; BF_r=AF
         end if  
         Open(nur,file=AF,form='UNFORMATTED')

        Case(6)             ! log-file

         Open(nu,file=AF_p)

        Case(4,11,12,13)    !  scratch files

         Open(nu,form='UNFORMATTED',status='SCRATCH')

        Case default
        
         write(*,*) 'nu,klsp =', nu,klsp

       End select

       End Subroutine Open_br

