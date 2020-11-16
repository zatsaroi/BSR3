!======================================================================
      Subroutine  matrix_J(jot)
!======================================================================
!     generation of interaction matrix for given J (jot = 2J+1 )
!----------------------------------------------------------------------
      Use bsr_ci
      USE term_LS, IL => ILterm, IS => ISterm

      Implicit none

      Integer, intent(in) :: jot
      Integer :: i,j, ip, icase
      Real(8) :: C, phase
      Real(8) :: F_SO(nterms,nterms), F_SS(nterms,nterms)
      Real(8), external :: Z_6j

!----------------------------------------------------------------------
! ... J-dependent factors:

      Do i=1,NTERMS
       Do j=1,NTERMS
       PHASE = (-1)**((IL(j)+IS(i)+jot-3)/2)
       F_SO(i,j) = PHASE * Z_6j(IL(i),IS(i),jot,IS(j),IL(j),3)
       F_SS(i,j) = PHASE * Z_6j(IL(i),IS(i),jot,IS(j),IL(j),5)
       End do
      End do

!----------------------------------------------------------------------
! ... read LS and overlap matrix:

      rewind(nuh); read(nuh) HM
      rewind(nuo); read(nuo) HS

!----------------------------------------------------------------------
! ... Include only those interactions for which |L-S| <= jot <= L+S

      Do j=1,NCFG
       ip=LSP(j)
       if(jot.lt.iabs(IL(ip)-IS(ip))+1.or. &
          jot.gt.     IL(ip)+IS(ip)-1) then
         Do i=j,NCFG; HM(i,j) = 0.d0; End do
         Do i=j,NCFG; HS(i,j) = 0.d0; End do
         HM(j,j) =  j*1.d8
         HS(j,j) = 1.d0
       end if
      End do

!----------------------------------------------------------------------
! ... add the J-dependent contribution:

      rewind(nur)
    1 read(nur,end=2) C,i,j,icase

      if(icase.eq.10) then
        phase = F_SS(LSP(i),LSP(j))
      else
        phase = F_SO(LSP(i),LSP(j))
      end if

      C = C * phase

      if(j.gt.i) then; ip=i;i=j;j=ip; end if

      HM(i,j) = HM(i,j) + C

      go to 1
    2 Continue

      End Subroutine  matrix_J

