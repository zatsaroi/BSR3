!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------

      USE bsr_breit
      USE spin_orbitals
      USE term_exp

      USE conf_LS,      only: ne
      USE coef_list,    only: ntrm
      USE symc_list_LS, only: IC_need, JC_need

      Implicit none 
      Integer :: i,j,k1,k2,iis,jjs,it,jt,ij,MLT2,MST2, met,ic,jc, is,js
      Integer, External :: IDEF_cme, DEF_ij

      Real(8) :: t1,t2,t3 

      Call CPU_time(t1)

      nc_new = 0
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

       if(IC_need(ic).eq.0) Cycle

!----------------------------------------------------------------------
!                                          cycle 2 over configurations:
      Call CPU_time(t3)

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

       ij=DEF_ij(ic,jc);  if(JC_need(ij).eq.0) Cycle      
              
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
        Call Add_res(nui,is,js)
        Call Add_it_oper(is,js)          
        Call Send_det_exp(i,iis,jjs)
       end if

!----------------------------------------------------------------------

      End do    ! over jc

      Call CPU_time(t2)                

      write(*,'(a,4i8,2f10.2,a)') 'ic,ic_total,kt,kdt', iis,ic_case,kt1,kdt1, &
        (t2-t3)/60, (t2-t1)/60, ' min.'

      End do    ! over ic, iis

!----------------------------------------------------------------------
! ... finish the calculations:

       Do 
        if(sum(ip_proc).eq.0) Exit 
        Call Get_res(j,is,js)        
        Call Add_res(nui,is,js)
        Call Add_it_oper(is,js)          
        ip_proc(j) = 0
        Call CPU_time(t2)                
        write(*,'(a,2i5,f10.2,a)') 'proc', j, sum(ip_proc), &
          (t2-t1)/60, ' min.'
       End do

! ... release the nodes:

       Do i=1,nprocs-1
        Call Send_det_exp(i,-1,-1)
       End do
       Call CPU_time(t2)                
       write(*,'(a,f10.2,a)') 'conf_loop is done', (t2-t1)/60, ' min.'

      End Subroutine Conf_loop


