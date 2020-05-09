!======================================================================
      Subroutine State_res
!======================================================================
!     extract and proceed the data from INT.BNK for specific integrals
!----------------------------------------------------------------------
      Use mpi

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
      Integer :: ibuf, nbuf, nnbuf, nccoef 
      Real(8) :: C,CC,CCC, t1,t2
      Integer, external :: Idef_itype, no_ic_LS, Check_idef
      Real(8), external :: Z_6j

      if(intercase.gt.0.and.intercase.ne.icase) Return

      if(myid.eq.0.and.pri.gt.0)  write(pri,*)

!----------------------------------------------------------------------
      if(bp_mode.eq.1) then
       Call State_res_bp
       Return
      End if

!----------------------------------------------------------------------
! ... nulify data in c_data module

      Call Initilize_c_data
!----------------------------------------------------------------------
! ... read coef.s from int.bnk to buffer:

      if(myid.eq.0) then; rewind(nub); read(nub) nbuf; end if

      nnbuf = 0
    1 if(myid.eq.0) Call Read_buffer(nbuf)  
      nnbuf=nnbuf+1

      if(interrupt.gt.0.and.nnbuf.le.interrupt) go to 1

      Call br_buffer_mpi(nbuf)

      interrupt = 0

!----------------------------------------------------------------------
! ... look through the buffer:

      nccoef = 0
      t1 =  MPI_WTIME()

      Do ibuf=1,ncbuf

      if(myid.eq.0.and.nprocs.gt.1) Exit

      Call Decode_INT (jcase,k,i1,i2,i3,i4,intb(ibuf)) 

      if(icase.ne.jcase) Cycle

      kpol = k 
      Select case(icase)
       Case(4,8,9,10); kpol=k-1
       Case(6,7,11);   kpol=0
      End select
      if(kpol.gt.mk) Cycle

! ... determine the range of states for given coeff. 

      it  = itb(ibuf);   jt = jtb(ibuf);
      is1 = IT_state1(it); js1 = IT_state1(jt)
      if(is1.eq.0.or.js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

      idf = idfb(ibuf)
      C = CBUF(ibuf)
!----------------------------------------------------------------------
! ... loop over all relevant states: 

      Do ik=is1,is2; is=IP_stat(ik); ich=IP_channel(is)
      Do jk=js1,js2; js=IP_stat(jk); jch=IP_channel(js)

! ... if this case belongs to given process:

       if(imycase(ich,jch).eq.0) Cycle

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
       end if

! ... we do not need anymore the configuration index 
! ... except pertuber (N+1)-electron configurations   

       i=0; if(ich.gt.nch) i=ich-nch
       j=0; if(jch.gt.nch) j=jch-nch
       ic = max(i,j)
       jc = min(i,j)

! ... find canonical form:

       Call Jsym_int(icase,j1,j2,j3,j4)

! ... find overlap factors with extracted continuum:  

       Call Det_fact(idf,np1,np2);  if(nndef.eq.0) Cycle

! ... send the final coefficients to archive:

       Do i = 1,nndef
        CCC = CC * Adef(i); io=iof(i); jo=jof(i)
        itype = Idef_itype(j1,j2,j3,j4,io,jo,ic,jc,k1,k2,k3,k4,ibi,jtype)
        if(check_target.eq.0) then
         if(check_idef(jtype,ich,jch).eq.0) Cycle
        end if
        Call Add_coef(CCC,kpol,k1,k2,k3,k4,itype)
        nccoef=nccoef + 1
       End do

      End do !  over js
      End do !  over is

      End do !  over ibuf

! ... generation of interaction matrix:

      if(myid.gt.0.and.pri.gt.0) write(pri,*) 'nccoef =',nccoef
      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      t2 =  MPI_WTIME()
      if(myid.gt.0.and.pri.gt.0)  write(pri,'(a,i4,T20,f10.1,a)') &
                  'Bufer:',nnbuf,(t2-t1)/60,' min'

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t2 =  MPI_WTIME()

      if(myid.eq.0.and.pri.gt.0)  write(pri,'(a,i4,T20,f10.1,a)') &
                  'Bufer:',nnbuf,(t2-t1)/60,' min'

      if( time_delay.gt.0.and.(t2-time0)/60 .gt. time_delay) then 
       interrupt = nnbuf; intercase=icase
       if(pri.gt.0)  write(pri,'(/a,3i10)') &
         'Interrupt, intercase, nnbuf:',interrupt,intercase,nnbuf
       Call Record_matrix
       if(myid.eq.0) then
        if(icase.eq.11) write(nuj) 0,0
        write(nuj) mk
        write(nuj) ACF
        write(nuj) htarg
        write(nuj) otarg
       end if
       Return
      end if

      if(nbuf.gt.0) go to 1  ! go for new data from data bank

      if(icase.eq.11) max_nbuf = nnbuf

      End Subroutine State_res


!======================================================================
      Subroutine Read_buffer(nbuf)  
!======================================================================
      Use bsr_mat

      Implicit none
 
      Integer :: nbuf,i1,i2

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
      Subroutine br_buffer_mpi(nbuf)  
!======================================================================
      Use mpi
      Use bsr_mat

      Implicit none
      Integer :: nbuf

      Call MPI_BCAST(nbuf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncbuf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(CBUF,ncbuf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(itb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jtb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(intb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(idfb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_buffer_mpi  



!======================================================================
      Subroutine State_res_bp
!======================================================================
!     extracts the data from int_list (after BSR_BMAT) for specific case
!     and update the Overlap/Hamiltoian matrix
!     (it is a replacement of State_res in case of usinf BSR_BMAT)
!----------------------------------------------------------------------
      Use mpi
      Use bsr_mat, k1 => itb, k2 => jtb, k3 => intb, k4 => idfb
      Use c_data, only: ntype,kpol1,kpol2

      Implicit none
      Integer :: i,jcase,kpol,itype, nc, ich,jch
      Real(8), external :: Cdata_occupation

!----------------------------------------------------------------------
! ... nulify data in c_data module

      Call Initilize_c_data

!----------------------------------------------------------------------
! ... read coef.s from int.bnk to buffer:

      if(myid.eq.0)  rewind(nub)

  10  Continue   

      if(myid.eq.0) then

   1  nc = 0
      read(nub,end=2,err=1) jcase,kpol,itype,nc             ! ???
      if(jcase.le.0) go to 1

      if(kpol.lt.kpol1.or.kpol.gt.kpol2)  Call Stop_mpi(pri,0,'kpol in buffer - ???')
      if(itype.lt.1.or.itype.gt.ntype)  Call Stop_mpi(pri,0,'itype in buffer - ???')

      if(nc.le.0) go to 1
      if(nc.gt.maxnc) then
       write(pri,*) 'nc  =',nc
       write(pri,*) 'maxnc =',maxnc
       Call Stop_mpi(pri,0,'nc.gt.maxnc')
      end if

      read(nub)  (cbuf(i),i=1,nc)
      read(nub)  (k1(i),i=1,nc)
      read(nub)  (k2(i),i=1,nc)
      read(nub)  (k3(i),i=1,nc)
      read(nub)  (k4(i),i=1,nc)

      if(jcase.ne.icase) go to 1

    2 Continue

      end if

      Call MPI_BCAST(nc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(nc.eq.0) then
        if(myid.eq.0) Return
        go to 20
      end if

      Call MPI_BCAST(itype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kpol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(CBUF,nc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(k1,nc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(k2,nc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(k3,nc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(k4,nc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.eq.0) go to 10
!----------------------------------------------------------------------

      Do i = 1,nc
       Call Def_ich_jch(itype,k1(i),k2(i),k3(i),k4(i),ich,jch)
       if(imycase(ich,jch).eq.0) Cycle
       Call Add_coef(Cbuf(i),kpol,k1(i),k2(i),k3(i),k4(i),itype)
      End do

      if(Cdata_occupation().lt.0.80) go to 10

      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      go to 10

!----------------------------------------------------------------------
  20  Continue

! ... generation of interaction matrix:

      Do kpol = kpol1,kpol2
       Do itype = 1,ntype
        Call Add_matrix(itype,kpol)
       End do
      End do

      End Subroutine State_res_bp


!========================================================================
      Subroutine Def_ich_jch(itype,k1,k2,k3,k4,ich,jch)
!========================================================================
!     used in State_res_bp to find channels for the given coefficient
!-----------------------------------------------------------------------

      Use bsr_mat
      Use conf_LS 

      Implicit none
      Integer, intent(in) :: itype,k1,k2,k3,k4
      Integer, intent(out) :: ich,jch
      Integer :: ic,jc,io,jo,i1,i2,j1,j2

      Select case(itype)

      Case(1)
       ic=-k3; jc=-k4
       if(ic.gt.0.and.jc.gt.0) then                !   (ic,jc)
        ich = ic + nch
        jch = ic + nch
       elseif(ic.lt.0.and.jc.gt.0) then            !   <i|.> ic
        io = -ic;  i1=io/ibo
        ich=iech(i1); jch = jc + nch
       elseif(ic.lt.0.and.jc.lt.0) then            !   <i|.> <.|j>
        io = -ic;  i1=io/ibo
        jo = -jc;  j1=jo/ibo
        ich = iech(i1); jch = iech(j1)
       elseif(ic.lt.0.and.jc.eq.0) then            !   <i|j>         
        io = -ic;  i1=io/ibo; i2=mod(io,ibo)
        ich = iech(i1); jch = iech(i2)
       else
        Call Stop_mpi(0,0,'problems to extract channel for itype=1')
       end if

      Case(2,3,4,5)                               
        ich=k3; io=k4
        if(io.gt.0) then 
          j1 = io/ibo; jch=iech(j1)
        elseif(io.lt.0) then
          jch=-io+nch
        else
          Call Stop_mpi(0,0,'problems to extract channel for itype=2,3,4,5')
        end if

      Case(6,7,8,9)                         
        ich=k1/ibi; jch=mod(k1,ibi)

      Case Default

       Call Stop_mpi(0,0,'unknown itype to extract channel')

      End Select

      End  Subroutine Def_ich_jch

