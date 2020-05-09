!======================================================================
      Subroutine Find_channel_label(ich,jch,is,E,Lab)
!======================================================================

      Use bsr_hd
      Use target
      Use channel
      Use conf_LS

      Implicit none
      Integer :: is,ich,jch,i,ic,ic1,ic2,nodes
      Real(8) :: S,E
      Character(*) ::  Lab

      if(ich.le.nch) then
       i = ich
       ic1=1; if(i.gt.1) ic1=ipconf(i-1)+1
       ic2=ipconf(i)
      else
       i = ich-nch
       ic1=ippert(i-1)+1
       ic2=ippert(i)
      end if

      i=0; S=0.d0
      Do ic=ic1,ic2
       if(abs(WC(ic)).lt.S) Cycle
       S=abs(WC(ic)); i=ic
      End do
      Call Get_cfg_LS(i)

      if(ich.le.nch) then
       nn(no) = ICHAR('k')-ICHAR('1')+1
       if(i.gt.0.and.Etarg(iptar(ich)).gt.E) then
        if(jch.gt.1) then
         nn(no) = ICHAR('n')-ICHAR('1')+1
        else
         Call Find_channel_nodes(ich,is,nodes)
         nn(no) = nodes + lch(ich) + 1
        end if
       end if
      end if

      Call Label_c(Lab,1,0)

      End Subroutine Find_channel_label


