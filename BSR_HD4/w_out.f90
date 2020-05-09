!======================================================================
      Subroutine W_OUT
!======================================================================
!     define the contributions of channels and pertuber configurations
!     for each solution
!----------------------------------------------------------------------
!     Original egeinproblem:
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
      Use conf_LS
      Use target
      Use channel,      only: iptar
      Use spline_param, only: ns

      Implicit none
      Integer :: i,i1,i2,ich,jch,is,ii,it,ms
      Real(8) :: S,S_ch,S_pt,E_Ry
      Character(64) ::  Labl(kch+kcp), Labl1
      Character(10) ::  AS
      Real(8) :: WT(kch+kcp),cval(khm)
      Integer :: iprm(kch+kcp),n_eff(kch)

      write(pri,'(/86(''-'')/)')

!----------------------------------------------------------------------
! ... define the configurations for the labels

      ii = INDEX(AF_cfg,'.',BACK=.TRUE.); AF = AF_cfg(1:ii)//ALSP
      Open(nuc,file=AF,status='OLD')
      ncfg=0; lcfg=0; Call Add_conf_LS(nuc,0) 
      Close(nuc)

!----------------------------------------------------------------------
! ... open w.nnn file:

      i = INDEX(AF_w,'.'); AF = AF_w(1:i)//ALSP
      Open(nuw,file=AF,form='UNFORMATTED')
      rewind(nuw)
      write(nuw) kch, kcp, khm

! ... define and record weights:

      Do is = 1,khm
       
       cval(1:khm) = a(1:khm,is) * a(1:khm,is)

       ! ... weights of channels:

       Do ich = 1,kch
        i1=ipsol(ich-1)+1; i2=ipsol(ich); WT(ich)=SUM(cval(i1:i2))
       End do

       ! ... weights of pseudostates:

       if(kcp.gt.0) WT(kch+1:kch+kcp)=cval(ipsol(kch)+1:khm)

       ! ... find channel with maximum contribution:

       nwt = kch+kcp 
       Call SORTA(nwt,WT,iprm);  isol(is) = iprm(1)

       ! ... find label:

       Labl = ' '
       Do jch=1,kch; ich=iprm(jch)
        Call Find_channel_label(ich,jch,is,eval(is),Labl(ich))
        if(jch.eq.1) Labl1=Labl(ich)
        if(ich.le.kch) n_eff(ich) = nn(no)
       End do

       ! ... record weights for given solution:

       write(nuw) WT,isol(is),eval(is),n_eff,Labl1,iprm

       ! ... print results:

       if(cwt.gt.0.d0) then
        if(msol.gt.0.and.is.gt.msol) Cycle

        S=SUM(WT(:)); S_ch=SUM(WT(1:kch)); S_pt=S-S_ch
        E_Ry = (eval(is)-etarg(1))*2

        write( pri,'(/i5,a,a,f16.8,a,f16.8,a,f8.5,a,f8.5/)' ) &
          is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry, &
           '  C_ch =',S_ch, '  C_pt =',S_pt

        Do jch=1,nwt; ich=iprm(jch); if(abs(WT(ich)).lt.CWT) Exit

        if(ich.gt.kch) then
         AS = 'perturber:'
         write(pri,'(a,i6,f11.5,5x,a)') AS,ich-kch,WT(ich),TRIM(Labl(ich))
        else
         AS = 'continium:'
         it = iptar(ich)
         if(Etarg(it).gt.eval(is))  AS = 'closed ch:'
         write(pri,'(a,i6,f11.5,5x,a)') AS,ich,WT(ich),TRIM(Labl(ich))
        end if 

       End do
      end if

      End do  ! over solutions 'is'

      End Subroutine W_OUT


