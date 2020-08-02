!======================================================================
      Subroutine rsol_out
!======================================================================
! ... output the R-matrix solutions and 
! ... find the surface amplitudes
!----------------------------------------------------------------------
      Use blacs
      Use bsr_hd     
      Use spline_param, only: ns

      Implicit none
      Real(8) :: vb(nhm)
      Integer :: i,j,i1,i2,j1,j2, is,ich 

      Call CPU_time(t0)

      if(itype.eq.1.and.io_processor) then
       i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALSP
       Open(nur,file=AF,form='UNFORMATTED')
       write(nur) nhm,khm,kch,kcp,ns
       write(nur) (eval(i),i=1,khm)
      end if

      if(itype.ge.0.and.io_processor) then
       if(allocated(WMAT)) deallocate(WMAT)
       Allocate(WMAT(kch,khm))
       write(*,'(a,f10.2,a)') 'WMAT: ',kch*khm*8.0/(1024*1024),'  Mb'
      end if

      Do is=1,khm

       Call pdgeadd ('notrans', khm, 1, one, z, 1,is,descz, &
                                       zero, v, 1,1, descv)
       Call BLACS_BARRIER (ctxt, 'all')

       if(.not.io_processor) Cycle      

       vb = 0.d0
       Do ich = 1,kch; i1=(ich-1)*ns+1; i2=ich*ns
        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
        Do j=j1,j2
         vb(i1:i2) = vb(i1:i2) + v(j)*bb(1:ns,j)     
        End do
        WMAT(ich,is) = vb(i2)
       End do
       if(kcp.gt.0) vb(kch*ns+1:nhm)=v(ksol+1:khm)
       if(itype.eq.1) write(nur) (vb(i),i=1,nhm)

      End do
      
      if(itype.eq.1.and.io_processor) close(nur)

      if(io_processor) write(*,'(/a)') 'rsol_out: R-matrix solutions done'
      call BLACS_BARRIER (ctxt, 'all')

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'RSOL_out:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'RSOL_out:,', (t1-t0)/60, ' min.'
      end if

      End Subroutine rsol_out


