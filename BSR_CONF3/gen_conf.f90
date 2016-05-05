!====================================================================
      Subroutine GEN_CONF
!====================================================================
!     generates all configurations from given configuration 'ic'
!     in module CONFIGS by adding the orbital 'ie' 
!     so to obtain total term IS,IL
!--------------------------------------------------------------------
      Use bsr_conf; Use conf_LS; Use orb_LS

      Call make_coupling

      ie=ie_comp; n=NEF(ie); l=LEF(ie); k=KEF(ie)
      ILT=Ltotal; IST=Stotal

! ... find if the same orbital is in the configuration ?

      no = no - 1         
      ii=0                                
      Do i=1,no
       if(n.ne.nn(i)) Cycle
       if(l.ne.ln(i)) Cycle
       if(k.ne.kn(i)) Cycle        
       ii=i; Exit
      End do

      if(ii.gt.0) then  ! ... case of the orbital trap on the existing shell   

       if(iq(ii).le.4*ln(ii)+1) then
        iq(ii)=iq(ii)+1
        IAP=LS(ii,1); ILP=LS(ii,2); ISP=LS(ii,3)
        IL_trap=ILP; IS_trap=ISP
        insert = -ii; S_cfp = 1.d0; S_recup = 1.d0
        Call SUM_TERM(ii)
       end if 
        
      else              ! ... find the position for adding orbital        

       Do i=1,no
        in = Ifind_nlk(nn(i),ln(i),kn(i),1)

        ip = (ie-in)/iabs(ie-in)
        i1 = min(in,ie)
        i2 = max(in,ie)
        if(IORT(i1,i2).eq.0) IORT(i1,i2) = 1
        if(IORT(i1,i2)*ip.lt.0) then; ii=i; Exit; end if

         if(ie.gt.in) Cycle
         ii = i; Exit

       End do

       if(ii.eq.0) then       ! ... add the new top shell     

        no=no+1
        nn(no)=n; ln(no)=l; iq(no)=1; kn(no)=k
        LS(no,1)=1; LS(no,2)=l+l+1; LS(no,3)=2
        LS(no,4)=ILT; LS(no,5)=IST
        insert = 0; S_cfp = 1.d0; S_recup = 1.d0 
        Call Check_cfg

       else                   ! ... insert new shell      

       Do i=no,ii,-1; i1=i+1
        nn(i1)=nn(i); ln(i1)=ln(i); iq(i1)=iq(i); kn(i1)=kn(i)
        LS(i1,:)=LS(i,:)
       End do
       no=no+1; nn(ii)=n; ln(ii)=l; iq(ii)=1; kn(ii)=k
       insert = ii; S_cfp = 1.d0; S_recup = 1.d0
       IAP=0; ILP=1; ISP=1
       Call SUM_TERM(ii)

       end if

      end if
      
      End Subroutine GEN_CONF


!=======================================================================
      Subroutine Sum_Term(i)
!=======================================================================
!     generates all terms for given shell 'i' 
!-----------------------------------------------------------------------
      Use bsr_conf
      Use conf_LS
      
      Implicit none
      Integer :: i,j,m,ii,jp,i1,i2,i3
      Integer, External :: Iterm_LS
      Real(8), External :: Zgen

      m = Iterm_LS(ln(i),iq(i),-1,i1,i2,i3)
      jp = Iterm_LS(ln(i),iq(i)-1,0,IAP,ILP,ISP)
      Do j=1,m
       ii = Iterm_LS(ln(i),iq(i),j,LS(i,1),LS(i,2),LS(i,3))
       if(insert.lt.0) S_cfp = zgen(iq(i)-1,ln(i),jp,j)
       CALL Sum_Iterm(i)
      End do

      End Subroutine Sum_Term


!=======================================================================
      Subroutine Sum_Iterm(istart)
!=======================================================================
!     generates all intermediate terms begining from shell 'istart'
!     which are concistent with the total TERM LS
!-----------------------------------------------------------------------
      Use bsr_conf
      USE conf_LS

      Integer :: LL_min(msh),LL_max(msh),SS_min(msh),SS_max(msh)

      i1 = istart                  ! i1 - low  limit
      i2 = no                      ! i2 - high limit in array LS(...)
      if(istart.eq.1) then
       i1=2
       LS(1,4)=LS(1,2)
       LS(1,5)=LS(1,3)
      end if

      if(no.eq.1) then
       if(LS(no,4).eq.ILT.and.LS(no,5).eq.IST)  Call Check_cfg
        Return
      end if

      i=i1
    1 j1=i-1
      j2=i

      LL_min(i)=IABS(LS(j1,4)-LS(j2,2))+1
      LL_max(i)=     LS(j1,4)+LS(j2,2) -1
      SS_min(i)=IABS(LS(j1,5)-LS(j2,3))+1
      SS_max(i)=     LS(j1,5)+LS(j2,3) -1
      LS(i,4)=LL_min(i)
      LS(i,5)=SS_min(i)

    2 if(i.lt.i2) then
         i=i+1
         go to 1
      else
         if(LS(no,4).eq.ILT.and.LS(no,5).eq.IST)  Call Check_cfg
      end if

    3 if(LS(i,5).lt.SS_max(i)) then
         LS(i,5) = LS(i,5) + 2
         go to 2
      elseif(LS(i,4).lt.LL_max(i)) then
         LS(i,4) = LS(i,4) + 2
         LS(i,5) = SS_min(i)
         go to 2
      else
         if(i.eq.i1) go to 4
         i=i-1
         go to 3
      end if

    4 Return

      End Subroutine Sum_Iterm


