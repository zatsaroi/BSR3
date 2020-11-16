!======================================================================
      Subroutine g_factor(jot,k,kk,gvJ,gvLS) 
!======================================================================
!     jot -  2J+1
!     k   -  index of solution
!     kk  -  leading configuration
!     gvLS - pure LS g-factor from the leading configuration
!     gvJ  - calculated g-factor
!----------------------------------------------------------------------
      Use bsr_ci
      Use term_LS
 
      Implicit none

      Integer :: i,j,it,k,kk,jot, ILT,IST
      Real(8) :: c,sj,sl,ss,sg,gvLS,gvj 
      Real(8) :: Scorr = 1.00232
      
! ... gvLS-factor:

      i = LSP(kk); ILT = ILterm(i); IST = ISterm(i)
      sj = jot - 1; sj=sj/2
      sl = ILT - 1; sl=sl/2
      ss = IST - 1; ss=ss/2
      gvLS = 1.d0 
      if(sj.ne.0.d0) gvLS = gvLS + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 

! ... gvJ-factor:

      if(nort.lt.0) then

       gvJ = 0.d0
       Do i=1,ncfg; it = LSP(i); ILT=ILterm(it); IST=ISterm(it)
        sj = jot - 1; sj=sj/2
        sl = ILT - 1; sl=sl/2
        ss = IST - 1; ss=ss/2
        sg = 1.d0
        if(sj.ne.0.d0) sg = sg + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 
        gvJ = gvJ + HM(i,k)**2 * sg
       End do

      else

       gvJ = 0.d0
       Do i=1,ncfg; it = LSP(i); ILT=ILterm(it); IST=ISterm(it)
        sj = jot - 1; sj=sj/2
        sl = ILT - 1; sl=sl/2
        ss = IST - 1; ss=ss/2
        sg = 1.d0
        if(sj.ne.0.d0) sg = sg + Scorr*(sj*(sj+1)-sl*(sl+1)+ss*(ss+1))/(2*sj*(sj+1)) 
        Do j=1,i; if(HS(i,j).eq.0.d0) Cycle; C=1.d0; if(i.ne.j) C=2.d0
         gvJ = gvJ + HS(i,j)*HM(i,k)*HM(j,k) * sg * C
        End do
       End do

      end if

      End Subroutine g_factor

