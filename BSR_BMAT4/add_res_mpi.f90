!======================================================================
      Subroutine Add_res(is_conf,js_conf)
!======================================================================                                              
!     extracts the data from int_bnk for specific case
!----------------------------------------------------------------------
      Use mpi

      Use bsr_breit
      Use bsr_mat
      Use conf_LS 
      Use new_dets;  Use new_defs
      Use term_exp,  only: kt1,kt2, IP_kt1, IP_kt2 
      Use coef_list, only: ncoef,coef,idfc,intc
      Use c_blocks,  only: ntype

      Implicit none
      Integer :: is_conf,js_conf
      Integer :: i,j,k,it,jt,is,js,is1,is2,js1,js2,ik,jk,ich,jch,ip1,ip2,ik1,ik2
      Integer :: jcase,kpol,itype,jtype, irecord,icase,kterm,kcoef,int,ii      
      Integer :: i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,ic,jc,io,jo,idf, m1,m2,m3,m4
      Integer :: ILT1,ILT2,IST1,IST2
      Real(8) :: C,CC,CCC
      Integer, external :: Idef_itype, no_ic_LS, Check_idef,Ifind_type
      Real(8), external :: Z_6j,c_blocks_occupation

!----------------------------------------------------------------------
! ... loop over terms:

      kterm = 0
      Do ik1=1,kt1; it=IP_kt1(ik1) 
      Do ik2=1,kt2; jt=IP_kt2(ik2)  
       if(is_conf.eq.js_conf.and.it.gt.jt) Cycle
       kterm = kterm + 1

! ... determine the range of states for terms:. 

      is1 = IT_state1(it); js1 = IT_state1(jt)
      if(is1.eq.0.or.js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

!----------------------------------------------------------------------
! ... loop over all relevant states: 

      Do ik=is1,is2; is=IP_stat(ik); ich=IP_channel(is)
!       if(my_channel(ich).eq.0) Cycle
      Do jk=js1,js2; js=IP_stat(jk); jch=IP_channel(js)
!       if(my_channel(ich).eq.0) Cycle

       if(ipert_ch.eq.0) then
        if(ich.le.nch.and.jch.gt.nch) Cycle
        if(jch.le.nch.and.ich.gt.nch) Cycle
       end if

! ... consider only low-half of interaction matrix

       if(it.eq.jt.and.js.gt.is) Cycle                   

! ... loop for angular coefficients

      Do kcoef = 1,ncoef

       C = coef(kterm,kcoef)
       if(abs(C).lt.Eps_C) Cycle

       idf = idfc(kcoef)
       int = intc(kcoef)

       Call Decode_INT (icase,k,i1,i2,i3,i4,intc(kcoef))        

       kpol = k 
       Select case(icase)
        Case(4,8,9,10); kpol=k-1
        Case(6,7,11);   kpol=0
       End select
       if(kpol.gt.mk) Cycle

! ... restriction of two-electron rel. matrix elements:

       Select case(icase)
        Case(3,4,8,9,10); if(abs(WC(is)*WC(js)).lt.Eps_soo) Cycle
       End Select

! ... include the expansion coefficients 

       C = coef(kterm,kcoef)
       if(ich.eq.jch.and.is.ne.js) C = C + C
       CC = C * WC(is) * WC(js)
       if(abs(CC).lt.Eps_C) Cycle

! ... define integral for specific orbitals

       ip1=IP_state(is)
       no1=no_ic_LS (is); np1(1:no1)=IP_orb(ip1+1:ip1+no1)

       ip2=IP_state(js)
       no2=no_ic_LS (js); np2(1:no2)=IP_orb(ip2+1:ip2+no2)

! ... define integral for specific orbitals

       if(icase.eq.11) then
        j1 = 1; j2 = 1; j3 = 1; j4 = 1
       else
        j1=np1(i1); j2=np1(i2); j3=np2(i3); j4=np2(i4)
       end if

! ... one-electron integrals:

       if(icase.eq.6.or.icase.eq.7) kpol = lbs(j1)
       if(icase.eq.7.and.kpol.gt.mlso) Cycle

! ... J-dependence for relativistic ccorrections 
       
       if(icase.gt.6.and.icase.lt.11) then
        Call Term_ic (is,ILT1,IST1)
        Call Term_ic (js,ILT2,IST2)
        k = 3; if(icase.eq.10) k = 5
        CC = CC * (-1)**((ILT2+IST1+jpar-3)/2)*   &
             Z_6j(ILT1,IST1,jpar,IST2,ILT2,k)
        if(abs(CC).lt.Eps_C) Cycle
       end if

! ... we do not need anymore the configuration index, 
! ... only pertuber index if any:   

       i=0; if(ich.gt.nch) i=ich-nch
       j=0; if(jch.gt.nch) j=jch-nch
       ic = max(i,j)
       jc = min(i,j)

! ... find canonical form:

       Call Jsym_int(icase,j1,j2,j3,j4)

! ... find overlap factors with extracted continuum:  

       Call Ndet_fact(idf,np1,np2); if(nndef.eq.0) Cycle

! ... send the final coefficients to archive:

       Do i = 1,nndef
        CCC = CC * Adef(i); io=iof(i); jo=jof(i)
        itype = Idef_itype(j1,j2,j3,j4,io,jo,ic,jc,k1,k2,k3,k4,ibi,jtype)
        if(check_target.eq.0) then
         if(check_idef(icase,jtype,ich,jch).eq.0) Cycle
        end if
        ii = Ifind_type(icase,kpol,itype)
        Call Add_coef_cblock(CCC,k1,k2,k3,k4,ii)
       End do            

      End do   ! over coefficients (kcoef)

      End do; End do ! over states (is,js)

      t2 = MPI_WTIME()    
      if( time_limit.gt.0.d0 .and. (t2-t0)/60.gt.time_limit*1.1) then
       write(*,*) 'failed at kterm,kt1,kt2', kterm,kt1,kt2
       go to 10
      end if 

      End do; End do ! over terms  (it,jt)


!-----------------------------------------------------------------------
! ... record coeficients 

      Do itype = 1,ntype
        Call Collect_coef(itype)
      End do

      Return

   10 Continue

      is_conf = -is_conf
      js_conf = -js_conf

      End Subroutine Add_res

             