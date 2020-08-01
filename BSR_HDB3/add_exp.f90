!======================================================================
      Subroutine Add_exp
!======================================================================
!     introduction of experimental energies:
!----------------------------------------------------------------------
      Use bsr_hd
      Use blacs
      Use channel, only: iptar
      Use target,  only: Etarg

      Implicit none
      Real(8) :: S
      Integer :: i, i1, it, idim, ich

      Call CPU_time(t0)

      if(io_processor.and.debug.gt.0) & 
       write(pri,'(/a/)') 'Experimental energies:'  

      Do ich = 1,kch
        
      add = zero
      idim=ipsol(ich)-ipsol(ich-1)

      if(io_processor) then

       it=iptar(ich)
       S=E_exp(it)-Etarg(it)
       Do i=1,idim; add(i,i)=S; End do
       if(debug.gt.0) &
       write(pri,'(2i5,2F16.8,f12.6)') ich,it, &
             E_exp(it),Etarg(it),S*27.2112
      end if

      Call BLACS_BARRIER (ctxt, 'All')

      i1=ipsol(ich-1)+1
      Call pdgeadd ('notrans', idim, idim,        &
                      one,  add,  1,  1, descadd, &
                      one,    a, i1, i1, desca    )
      End do

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'Add_exp:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'Add_exp:,', (t1-t0)/60, ' min.'
      end if

      End Subroutine Add_exp

