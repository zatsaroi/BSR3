!======================================================================
      Subroutine W_OUT
!======================================================================
!     define the contributions of channels and pertuber configurations
!     for each solution
!----------------------------------------------------------------------
!     Original egenproblem:
!
!     H C = S C E,   E - diagonal egenvalue matrix
!                    S - overlap matrix
!
!----------------------------------------------------------------------
!     Using squire-root matrix:
!
!     C' H C = (C' S C) E  ->  C' S C = 1  =>  (C'S^1/2) (S^1/2 C) = 1
!
!     S A = s A;  S^1/2 = A s(1/2) A'
!
!-----------------------------------------------------------------------
!     Using Cholesky decomposition:
!
!     S = L L' ->     H C = L L' C E
!                    [inv(L) H inv(L')] L'C = L'C E
!                    (L'C)' (L'C) = 1
!
!     supposed solutions are given as (L'C) in matrix A
!-----------------------------------------------------------------------
      Use bsr_hd
      Use blacs
      Use channel, only: iptar
      Use conf_ls
      Use target
      Use spline_param, only: ns

      Implicit none
      Integer :: i,i1,i2,ich,jch,is,it,ms
      Real(8) :: S,S_ch,S_pt,E_Ry
      Character(64) ::  Labl
      Character(10) ::  AS
      Real(8), allocatable :: WT(:),cval(:)
      Integer, allocatable :: iprm(:)
      Integer, external :: Icheck_file

      Call CPU_time(t0)

! ... local allocations:

      nwt = kch+kcp
      if(allocated(cval)) deallocate(cval); allocate(cval(khm))
      if(allocated(isol)) deallocate(isol); allocate(isol(khm))
      if(allocated(wt  )) deallocate(wt  ); allocate(wt  (nwt))
      if(allocated(iprm)) deallocate(iprm); allocate(iprm(nwt))

!----------------------------------------------------------------------

      if(io_processor) then      
       i = INDEX(AF_w,'.'); AF = AF_w(1:i)//ALSP
       Open(nuw,file=AF,form='UNFORMATTED')
       rewind(nuw)
       write(nuw) kch, kcp, khm
      end if 

! ... define and record weights:

      Do is = 1,khm
       
       call pdgeadd ('notrans', khm, 1, one, z, 1,is,descz, &
                                       zero, v, 1,1, descv)
       call BLACS_BARRIER (ctxt, 'all')

       if(.not.io_processor) Cycle

       cval(1:khm) = v(1:khm) * v(1:khm)

       ! ... weights of channels:

       Do ich = 1,kch
        i1=ipsol(ich-1)+1; i2=ipsol(ich); WT(ich)=SUM(cval(i1:i2))
       End do

       ! ... weights of pseudostates:

       if(kcp.gt.0) WT(kch+1:kch+kcp)=cval(ipsol(kch)+1:khm)

       ! ... find channel with maximum contribution:

       jch=0; S=0.d0
       Do ich=1,kch+kcp
        if(S.gt.WT(ich)) Cycle; jch=ich; S=WT(ich)
       End do
       isol(is) = jch

       ! ... record weights for given solution:

       write(nuw) WT,jch

      End do  ! over solutions 'is'

!----------------------------------------------------------------------
! ... define the configurations for labels

      if(io_processor) then      
       i = LEN_TRIM(AF_cfg); AF = AF_cfg(1:i-3)//ALSP
       i = Icheck_file(AF)
       if(i.eq.1) then 
        Open(nuc,file=AF)
        ncfg=0; lcfg=0; Call Add_conf_LS(nuc,0)
        Close(nuc)
       else
        write(*,'(/a,a)') 'Can not find file ',AF
        fail = 1
       end if
      end if

      Call br_ipar(fail);   if(fail.eq.1) Return

!----------------------------------------------------------------------
! ... print results:

      if(CWT.gt.0.d0.and.io_processor) then

       write(pri,'(/86(''-'')/)')

       rewind(nuw)
       read(nuw) kch, kcp, khm

       ms=msol; if(msol.le.0.or.msol.gt.khm) ms=khm

       Do is = 1,ms
        read(nuw) WT,jch
        if(eval(is).gt.Emax) Cycle

        S=SUM(WT(:)); S_ch=SUM(WT(1:kch)); S_pt=S-S_ch
        E_Ry = (eval(is)-etarg(1))*2

        S = zero            ! closed-channel contribution 
        Do ich=1,nwt
         if(ich.gt.kch) Cycle
         it = iptar(ich)
         if(Etarg(it).lt.eval(is)) Cycle
         S = S + WT(ich)
        End do

        if(S_pt.gt.0.00001) then
         write( pri,'(/i5,a,a,f16.8,a,f16.8,5x,a,f10.5,5x,a,f10.5)' ) &
           is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry, &
           '  C_cl =',S,'  C_pt =',S_pt
        elseif(S.gt.0.00001) then
         write( pri,'(/i5,a,a,f16.8,a,f16.8,5x,a,f10.5,5x,a,f10.5)' ) &
           is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry, &
           '  C_cl =',S
        else
         write( pri,'(/i5,a,a,f16.8,a,f16.8,5x,a,f10.5,5x,a,f10.5)' ) &
           is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry
        end if

        Call SORTA(nwt,WT,iprm)

        Do jch=1,nwt; ich=iprm(jch); if(abs(WT(ich)).lt.CWT) Exit

         Call Find_channel_label(ich,jch,is,eval(is),Labl)

        if(ich.gt.kch) then
         AS = 'perturber:'
         write(pri,'(a,i6,f11.5,5x,a)') AS,ich-kch,WT(ich),TRIM(Labl)
        else
         AS = 'continium:'
         it = iptar(ich)
         if(Etarg(it).gt.eval(is))  AS = 'closed ch:'
         write(pri,'(a,i6,f11.5,5x,a)') AS,ich,WT(ich),TRIM(Labl)
        end if 

       End do
      End do

      end if ! over cwt

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'W_out:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'W_out:,', (t1-t0)/60, ' min.'
      end if

      End Subroutine W_OUT

