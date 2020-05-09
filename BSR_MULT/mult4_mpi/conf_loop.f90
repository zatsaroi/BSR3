!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------
      Use MPI

      Use mult_par

      Use spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      USE term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, ic_case, &
                               MLT1, MST1, MLT2, MST2

      Use conf_LS,      only: ne
      Use symc_list_LS, only: JC_need, IC_need, nsymc
      Use symt_list_LS, only: IT_done, ij
      Use coef_list,    only: ntrm,ctrm
      Use zoef_list,    only: nzoef

      Implicit none 

      Integer :: k1,k2,ic,jc,is,js,iis,jjs, it,jt, m,k, i,j
      Integer(8), external :: DEF_ij8

      Character(80) :: conf

      t1=MPI_WTIME()

!----------------------------------------------------------------------
!                                          cycle 1 over configurations:
      rewind(nud)
      Do iis=1,ic_case

       Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT1,MST1

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

      t2=MPI_WTIME()

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

       if(JC_need(DEF_ij8(ic,jc)).eq.0) Cycle      

!----------------------------------------------------------------------
! ...  define number of terms:

       ntrm = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle;  ntrm = ntrm + 1
       End do; End do 

!----------------------------------------------------------------------
! ...  JT_oper:

       if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
       Allocate(JT_oper(ntrm),CT_oper(ntrm))

       k = 0; m = 0; JT_oper=0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle;  k=k+1
        ij=DEF_ij8(it,jt) 
        if(IT_done(ij).ne.0) Cycle 
        JT_oper(k) = 1
        m = m + 1
       End do; End do 

       if(m.eq.0) Cycle

!----------------------------------------------------------------------
! ...  calculations:

       m = 0
       Do i=1,nprocs-1
        if(ip_proc(i).ne.0) Cycle 
        Call Send_det_exp(i,iis,jjs)
        m = i 
        ip_proc(i) = 1
        Exit
       End do

       if(m.eq.0) then
        Call Get_res(i,is,js)        
        Call Add_res(nui,is,js)
        Call DEF_ic(is,js)          
        Call Send_det_exp(i,iis,jjs)
       end if

      End do    ! over jc

      t3=MPI_WTIME()                

      Call Symc_conf(ic,conf)
      write(*  ,'(a,i6,a,i6,a,i6,a,i6,F10.2,a,3x,a)') &
        ' ic=',ic,'/',nsymc,'  nterm=',kt1,'  ndet=', kdt1, &
          t3-t2,' sec.',trim(conf)
      End do    ! over ic

!----------------------------------------------------------------------
! ... finish the calculations:

      Do 
       if(sum(ip_proc).eq.0) Exit 
       Call Get_res(j,is,js)        
       Call Add_res(nui,is,js)
       Call DEF_ic(is,js)          
       ip_proc(j) = 0
       t2=MPI_WTIME()                
       write(*,'(a,2i5,f10.2,a)') 'proc', j, sum(ip_proc), &
         (t2-t3)/60, ' min.'
      End do

! ... release the nodes:

       Do i=1,nprocs-1
        Call Send_det_exp(i,-1,-1)
       End do

       t2=MPI_WTIME()                

       write(*,'(a,f10.2,a)') 'conf_loop is done', (t2-t1)/60, ' min.'

      End Subroutine Conf_loop


