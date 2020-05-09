!======================================================================
      Subroutine State_res
!======================================================================
!     extracts the data from int_bnk for specific case
!----------------------------------------------------------------------
      Use bsr_mat
      Use conf_LS 
      Use def_list; Use def_list    
      Use new_dets; Use new_defs
      Use c_data, only: ntype,kpol1,kpol2

      Implicit none
      Integer :: i,j,k,it,jt,is,js,is1,is2,js1,js2,ik,jk,ich,jch,ip1,ip2
      Integer :: jcase,kpol,itype,jtype
      Integer :: i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,ic,jc,io,jo,idf
      Integer :: ILT1,ILT2,IST1,IST2
      Integer :: ibuf, nbuf, nnbuf
      Real(8) :: C,CC,CCC, t1,t2

      Integer, external :: Idef_itype, no_ic_LS, Check_idef
      Real(8), external :: Z_6j

      if(intercase.gt.0.and.intercase.ne.icase) Return

!----------------------------------------------------------------------
! ... if we have coefficients after BSR_BMAT:

      if(bp_mode.eq.1) then;  Call State_res_bp; Return; End if

!----------------------------------------------------------------------
! ... nulify data in c_data module

      Call Initilize_c_data

!----------------------------------------------------------------------
! ... read coef.s from int.bnk to buffer:

      Call CPU_time(t1)
      if(myid.eq.0) then; rewind(nub); read(nub) nbuf; end if

      nnbuf = 0
    1 if(myid.eq.0) Call Read_buffer(nbuf)  

      nnbuf=nnbuf+1

      if(interrupt.gt.0.and.nnbuf.le.interrupt) go to 1

      interrupt = 0

!----------------------------------------------------------------------
! ... look through the buffer:

      Do ibuf=1,ncbuf

      Call Decode_INT (jcase,k,i1,i2,i3,i4,intb(ibuf)) 
      
      if(icase.ne.jcase) Cycle

      kpol = k 
      Select case(icase)
       Case(4,8,9,10); kpol=k-1
       Case(6,7,11);   kpol=0
      End select
      if(kpol.gt.mk) Cycle

! ... determine the range of states for given coeff. 

      it = itb(ibuf);   jt = jtb(ibuf)
      is1 = IT_state1(it); js1 = IT_state1(jt)
      if(is1.eq.0.or.js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

      idf = idfb(ibuf)
      C = CBUF(ibuf)
!----------------------------------------------------------------------
! ... loop over all relevant states: 

      Do ik=is1,is2; is=IP_stat(ik); ich=IP_channel(is)
       if(my_channel(ich).eq.0) Cycle
      Do jk=js1,js2; js=IP_stat(jk); jch=IP_channel(js)
       if(my_channel(ich).eq.0) Cycle

       if(ipert_ch.eq.0) then
        if(ich.le.nch.and.jch.gt.nch) Cycle
        if(jch.le.nch.and.ich.gt.nch) Cycle
       end if

! ... consider only low-half of interaction matrix

       if(it.eq.jt.and.js.gt.is) Cycle                   

! ... restriction of two-electron rel. matrix elements:

       Select case(icase)
        Case(3,4,8,9,10); if(abs(WC(is)*WC(js)).lt.Eps_soo) Cycle
       End Select

! ... include the expansion coefficients 

       C=CBUF(ibuf)
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

        if(icase.eq.7) CC = CC *zcorr      ! ???

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

       Call Det_fact(idf,np1,np2); if(nndef.eq.0) Cycle

! ... send the final coefficients to archive:

       Do i = 1,nndef
        CCC = CC * Adef(i); io=iof(i); jo=jof(i)
        itype = Idef_itype(j1,j2,j3,j4,io,jo,ic,jc,k1,k2,k3,k4,ibi,jtype)

        ! exclude the target only integrals:
        if(check_target.eq.0) then
         if(check_idef(jtype,ich,jch).eq.0) Cycle
        end if

        Call Add_coef(CCC,kpol,k1,k2,k3,k4,itype)
       End do            

      End do   !  over js
      End do   !  over is

      End do   !  over ibuf
!-----------------------------------------------------------------------
! ... generation of interaction matrix:

      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      Call CPU_time(t2)
      if(pri.gt.0.and.debug.gt.0)  write(pri,'(a,i4,a,f8.2,a)') &
       'Bufer:',nnbuf,' time:',(t2-t1)/60,' min'

      if( time_delay.gt.0.and.(t2-time0)/60 .gt. time_delay) then 
       interrupt = nnbuf; intercase=icase
       if(myid.eq.0.and.pri.gt.0)  write(pri,'(/a,3i10)') &
         'Interrupt, intercase, nnbuf:',interrupt,intercase,nnbuf
       Call Record_matrix
       if(icase.eq.11) write(nuj) 0,0
       write(nuj) mk
       write(nuj) ACF
       write(nuj) htarg
       write(nuj) otarg
       Return
      end if

      if(nbuf.gt.0) go to 1  ! go for new data from data bank

      End Subroutine State_res



!======================================================================
      Subroutine Read_buffer(nbuf)  
!======================================================================
      Use bsr_mat

      Implicit none
 
      Integer :: nbuf, i1,i2

      if(nbuf.gt.maxnc) then
       write(pri,*) 'nbuf  =',nbuf
       write(pri,*) 'maxnc =',maxnc
       Call Stop_mpi(0,0,'nbuf.gt.maxnc')
      end if
      if(nbuf.le.0) Return
      ncbuf =0
       
    1 i1 = ncbuf+1; i2=ncbuf+nbuf
      read(nub) cbuf (i1:i2)
      read(nub) itb  (i1:i2)
      read(nub) jtb  (i1:i2)
      read(nub) intb (i1:i2)
      read(nub) idfb (i1:i2)
      ncbuf = ncbuf + nbuf
      read(nub,end=2) nbuf
      if(ncbuf+nbuf.gt.maxnc) Return
      go to 1

    2 nbuf=0
      
      End Subroutine Read_buffer  
             

!======================================================================
      Subroutine State_res_bp
!======================================================================
!     extracts the data from int_bnk for specific case
!----------------------------------------------------------------------
      Use bsr_mat, k1 => itb, k2 => jtb, k3 => intb, k4 => idfb
      Use c_data, only: ntype,kpol1,kpol2
      Implicit none
      Integer :: i,jcase,kpol,itype, nc
      Real(8), external :: Cdata_occupation

! ... nulify data in c_data module

      Call Initilize_c_data

! ... read coef.s from int.bnk to buffer:

      rewind(nub)
    1 read(nub,end=2) jcase,kpol,itype,nc
      if(jcase.le.0) go to 1

      if(kpol.lt.kpol1.or.kpol.gt.kpol2)  Call Stop_mpi(0,0,'kpol in buffer - ???')
      if(nc.le.0) go to 1

      if(nc.gt.maxnc)  then
        maxnc = nc
        Deallocate(cbuf,k1,k2,k3,k4)
        Allocate(cbuf(maxnc),k1(maxnc),k2(maxnc),k3(maxnc),k4(maxnc)) 
        write(*,*) 'Warning: maxnc = ',maxnc
      end if

      read(nub)  (cbuf(i),i=1,nc)
      read(nub)  (k1(i),i=1,nc)
      read(nub)  (k2(i),i=1,nc)
      read(nub)  (k3(i),i=1,nc)
      read(nub)  (k4(i),i=1,nc)

      if(jcase.ne.icase) go to 1

      Do i = 1,nc
       Call Add_coef(Cbuf(i),kpol,k1(i),k2(i),k3(i),k4(i),itype)
      End do

      if(Cdata_occupation().lt.0.80) go to 1

      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      go to 1

    2 Continue

! ... generation of interaction matrix:

      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      End Subroutine State_res_bp

     








