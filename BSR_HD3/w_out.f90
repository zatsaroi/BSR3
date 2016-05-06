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
      Character(64) ::  Labl
      Character(10) ::  AS
      Real(8), allocatable :: WT(:),cval(:)
      Integer, allocatable :: iprm(:)

! ... local allocations:
 
      Allocate(wt(kch+kcp),cval(khm),iprm(kch+kcp))
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

       jch=0; S=0.d0
       Do ich=1,kch+kcp
        if(S.gt.WT(ich)) Cycle; jch=ich; S=WT(ich)
       End do
       isol(is) = jch

       ! ... record weights for given solution:

       write(nuw) WT,jch,eval(is)

      End do  ! over solutions 'is'
!----------------------------------------------------------------------
! ... define the configurations for the labels

      ii = INDEX(AF_cfg,'.',BACK=.TRUE.); AF = AF_cfg(1:ii)//ALSP
      Open(nuc,file=AF,status='OLD')
      ncfg=0; lcfg=0; Call Add_conf_LS(nuc,0) 
      Close(nuc)

!----------------------------------------------------------------------
! ... print results:

      if(CWT.gt.0.d0) then

       write(pri,'(/86(''-'')/)')

       rewind(nuw)
       read(nuw) kch, kcp, khm
       nwt = kch+kcp

       ms=msol; if(msol.le.0.or.msol.gt.khm) ms=khm

       Do is = 1,ms
        read(nuw) WT,jch
        S=SUM(WT(:)); S_ch=SUM(WT(1:kch)); S_pt=S-S_ch
        E_Ry = (eval(is)-etarg(1))*2

        write( pri,'(/i5,a,a,f16.8,a,f16.8,a,f8.5,a,f8.5/)' ) &
          is,'.','  E_au =',eval(is),'  E_Ry =',E_Ry, &
           '  C_ch =',S_ch, '  C_pt =',S_pt

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

      Deallocate(wt,cval,iprm)

      End Subroutine W_OUT


