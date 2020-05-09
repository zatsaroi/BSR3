!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------
      Use MPI

      Use bsr_breit
      Use spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      Use term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, MLT,MST, ic_case
      Use conf_LS,      only: ne
      Use coef_list,    only: ntrm
      Use symc_list_LS, only: IC_need, JC_need, nsymc

      Implicit none 
      Integer :: i,j,k1,k2,iis,jjs,it,jt,MLT2,MST2, m, met,ic,jc, is,js
      Integer, external :: IDEF_cme
      Integer, allocatable :: is_ic(:), kt_is(:), kdt_is(:)
      Character(80), allocatable :: conf_is(:)
      Integer(8) :: ij
      Integer(8), external :: DEF_ij8

      t1=MPI_WTIME()
      m = ic_case
      Allocate(is_ic(m), conf_is(m), kt_is(m), kdt_is(m))
!----------------------------------------------------------------------
!                                          cycle 1 over configurations:
      rewind(nud)
      Do iis=1,ic_case

       Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT,MST
       if(Allocated(IP_kt1)) Deallocate(IP_kt1)
       Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

       if(Allocated(C_det1)) Deallocate(C_det1)
       Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

       if(Allocated(IM_det1)) Deallocate(IM_det1)
       Allocate(IM_det1(ne,kdt1)); Read(nud) IM_det1

       if(Allocated(IS_det1)) Deallocate(IS_det1)
       Allocate(IS_det1(ne,kdt1)); Read(nud) IS_det1

       read(nud) NNsym1(1:ne)
       read(nud) Lsym1(1:ne)
       
       is_ic(iis) = ic
       kt_is(iis) = kt1
       kdt_is(iis) = kdt1
       
       Call Symc_conf(ic,conf_is(iis))

       if(IC_need(ic).eq.0) Cycle

       if(mkdt.gt.0.and.kdt1.gt.mkdt) Cycle
!----------------------------------------------------------------------
!                                          cycle 2 over configurations:
      t3=MPI_WTIME()

      rewind(nud)
      Do jjs=1,iis

       Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2
       if(Allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

       if(Allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

       if(Allocated(IM_det2)) Deallocate(IM_det2)
       Allocate(IM_det2(ne,kdt2)); Read(nud) IM_det2

       if(Allocated(IS_det2)) Deallocate(IS_det2)
       Allocate(IS_det2(ne,kdt2)); Read(nud) IS_det2

       read(nud) NNsym2(1:ne)
       read(nud) Lsym2(1:ne)

       if(MLT2.ne.MLT.or.MST2.ne.MST) Cycle      ! ???

       if(MLT.ne.min(ILT1,ILT2).or.MST.ne.min(IST1,IST2)) Cycle ! ???

       ij=DEF_ij8(ic,jc);  if(JC_need(ij).eq.0) Cycle      
              
!----------------------------------------------------------------------
! ...  define number of terms:

       ntrm = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(iis.eq.jjs.and.it.gt.jt) Cycle;  ntrm = ntrm + 1
       End do; End do 
  
!----------------------------------------------------------------------
! ...  joper and JT_oper:

       if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
       Allocate(JT_oper(ntrm,noper),CT_oper(ntrm,noper))

       if(IDEF_cme(iis,jjs).eq.0) Cycle 

!----------------------------------------------------------------------
! ...  calculations:

       met = 0
       Do i=1,nprocs-1
        if(ip_proc(i).ne.0) Cycle 
        Call Send_det_exp(i,iis,jjs)
        met = i 
        ip_proc(i) = 1
        Exit
       End do

       if(met.eq.0) then
        Call Get_res(i,is,js)        

        if(is.gt.0) then
         Call Add_res(nur,is,js)
         Call Add_it_oper(is,js)          
        end if

        Call Send_det_exp(i,iis,jjs)

       end if

!----------------------------------------------------------------------

      End do    ! over jc

      t2=MPI_WTIME()                

      write(*,'(a,4i8,2f10.2,a,5x,a)') 'ic,ic_total,kt,kdt', iis,ic_case,kt1,kdt1, &
        (t2-t3)/60, (t2-t0)/60, ' min.', conf_is(iis)

      if(time_limit.gt.0.d0.and.(t2-t1)/60.gt.time_limit) Exit

      End do    ! over ic, iis

!----------------------------------------------------------------------
! ... finish the calculations:

      Do 
       if(sum(ip_proc).eq.0) Exit 
       Call Get_res(j,is,js)        

       if(is.gt.0) then
        Call Add_res(nur,is,js)
        Call Add_it_oper(is,js)          
       end if
 
       t2=MPI_WTIME()                
       write(*,'(3i5,f10.2,a)') sum(ip_proc), is, js, (t2-t0)/60, ' min.'
       if(is.lt.0) then
        is = iabs(is); js = iabs(is)
        write(*,'(a,2i10,5x,a)') 'kt,kdt =',kt_is(is),kdt_is(is),conf_is(is)
        write(*,'(a,2i10,5x,a)') 'kt,kdt =',kt_is(js),kdt_is(js),conf_is(js)
       end if

       ip_proc(j) = 0
      End do

! ... release the nodes:

       Do i=1,nprocs-1
        Call Send_det_exp(i,-1,-1)
       End do
       t2=MPI_WTIME()                
       write(*,'(a,f10.2,a)') 'conf_loop is done', (t2-t1)/60, ' min.'

      End Subroutine Conf_loop


