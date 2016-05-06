!======================================================================
      Subroutine R_arg
!======================================================================
!     read arguments, first from file unit 'nu, then from comand line
!----------------------------------------------------------------------
      Use bsr_pol
    
      Implicit none
      Integer, external :: Icheck_file

      if(Icheck_file(AF_par).ne.0) then
       Open(nup,file=AF_par)
       Call Read_ipar(nup,'klsp' ,klsp )
       Call Read_ipar(nup,'ilzero',ilzero)
       Call Read_ipar(nup,'ibzero',ibzero)
       Call Read_ipar(nup,'nortb' ,nortb )
       if(nortb.gt.0) then
        Allocate(iortb(nortb)); iortb=0
        Call Read_iarray(nup,'iortb',nortb,iortb)
       end if
      end if

      Call Read_iarg('klsp'  ,klsp )
      Call Read_iarg('ilzero',ilzero)
      Call Read_iarg('ibzero',ibzero)
      Call Read_iarg('nortb' ,nortb)
      if(nortb.gt.0) then
       if(allocated(iortb)) Deallocate(iortb)
       Allocate(iortb(nortb)); iortb=0
       Call Read_iarr('iortb',nortb,iortb)
      end if

      write(ALSP,'(i3.3)') klsp

      if(nortb.gt.0) &
       open(nuq,form='UNFORMATTED',status='SCRATCH')

      End Subroutine R_arg

