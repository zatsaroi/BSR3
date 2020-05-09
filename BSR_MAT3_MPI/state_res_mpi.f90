!======================================================================
      Subroutine State_res
!======================================================================
!     extracts the data from INT.BNK for specific case
!----------------------------------------------------------------------

      Use mpi

      USE bsr_mat
      USE conf_LS
      USE orb_LS
      USE channel
      USE target
      USE cmdata
      USE dets; USE new_dets; Use new_defs
      USE z_core, only: mlz 
      USE bsr_matrix, only: imycase

      Implicit none

      Integer :: i,j,k,it,jt,is,js,is1,is2,js1,js2,ik,jk,ich,jch, nnbuf
      Integer :: ncfg_ch,nccoef,ibuf,jcase,int,ijt,jpol,jtype, ip1,ip2
      Integer :: i1,i2,i3,i4,j1,j2,j3,j4,kk1,kk2,kk3,kk4,ic,jc,io,jo,idf
      Integer :: ILT1,ILT2,IST1,IST2

      Real(8) :: C,CC,CCC, t1,t2,t3

      Integer, External :: Idef_itype, no_ic_LS
      Real(8), External :: Z_6j

!----------------------------------------------------------------------
!                                               prepare data in cmdata:
! ... initilize all blocks:

      Do i=1,nb
       ipblk(i) = (i-1)*mb + 1; jpblk(i) = -1
      End do

! ... assign one block to each type:

      i = 0
      Do k = -1,npol
       Do j = 1,ntype
        i = i + 1; iblk(k,j)=i; nblk(k,j)=1; kblk(k,j,1)=i; jpblk(i)=0
       End do
      End do

      if(i.gt.nb/2) Call Stop_mpi(0,2*i,'increase nb parameter')

      ncfg_ch=ipconf(nch);  ncp=ncfg-ncfg_ch

!----------------------------------------------------------------------
! ... read coef.s from int.bnk to buffer:

      t1 =  MPI_WTIME()

      if(myid.eq.0) Call Read_dets
      Call br_dets
      if(icase.eq.11.and.pri.gt.0) then
       write(pri,'(/a,2i8 )') 'ndet   = ',ndet,kdet
       write(pri,'( a,2i8 )') 'ndef   = ',ndef,kdef
       C = 2*ndet + kdet + 2*ndef + kdef
       C = C /(256*1024)
       write(pri,'(/a,F8.2,a/)') 'Det. memory: ',C,'  Mb'
      end if

      nnbuf = 0
    1 Call Read_buffer_mpi  

      t2 =  MPI_WTIME()

      nnbuf=nnbuf+1

!----------------------------------------------------------------------
! ... look through the buffer:

      nccoef = 0
      Do ibuf=1,ncbuf

      if(myid.eq.0) Exit

      Call Decode_INT (jcase,k,i1,i2,i3,i4,intb(ibuf)) 

      if(icase.ne.jcase) Cycle

      Select case(icase)
       Case(3,5);      kpol=k
       Case(4,8,9,10); kpol=k-1
       Case(6,7,11);   kpol=0
      End select
      if(kpol.gt.mk) Cycle

! ... determine the range of states for given coeff. 

      ijt = ijtb(ibuf)
      it  = ijt/ibc;   jt = mod(ijt,ibc)
      is1 = IT_state1(it); js1 = IT_state1(jt)
      if(is1.eq.0.or.js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

      idf = idfb(ibuf)
      C = CBUF(ibuf)
!----------------------------------------------------------------------
!                                        loop over all relevant states: 

      Do ik=is1,is2; is=IP_stat(ik); ich=IP_channel(is)
      Do jk=js1,js2; js=IP_stat(jk); jch=IP_channel(js)

! ... if this case belongs to given process:

       if(imycase(ich,jch).eq.0) Cycle

       nccoef=nccoef+1

! ... consider only low-half of interaction matrix

       if(it.eq.jt.and.js.gt.is) Cycle                   

! ... restriction of two-electron rel. matrix elements:

       Select case(icase)
        Case(3,4,8,9,10); if(abs(WC(is)*WC(js)).lt.Eps_soo) Cycle
! Case(7); ! if(iptar(ich).gt.2.or.iptar(jch).gt.2) Cycle  ???
       End Select

! ... Include the expansion coefficients 

       C=CBUF(ibuf)
       if(ich.eq.jch.and.is.ne.js) C = C + C
       CC = C * WC(is) * WC(js)
       if(abs(CC).lt.Eps_C) Cycle

! ... Define integral for specific orbitals

       ip1=IP_state(is)
       no1=no_ic_LS (is); np1(1:no1)=IP_orb(ip1+1:ip1+no1)

       ip2=IP_state(js)
       no2=no_ic_LS (js); np2(1:no2)=IP_orb(ip2+1:ip2+no2)

       if(icase.eq.11) then
        j1 = 1; j2 = 1; j3 = 1; j4 = 1
       else
        j1=np1(i1); j2=np1(i2); j3=np2(i3); j4=np2(i4)
       end if

if(debug.gt.0) then
 if(j1.lt.1.or.j1.gt.nwf) write(pri,*) 'warning, j1:',j1
 if(j2.lt.1.or.j2.gt.nwf) write(pri,*) 'warning, j2:',j2
 if(j3.lt.1.or.j3.gt.nwf) write(pri,*) 'warning, j3:',j3
 if(j4.lt.1.or.j4.gt.nwf) write(pri,*) 'warning, j4:',j4
end if

! ... one-electron integrals:

       if(icase.eq.6.or.icase.eq.7) kpol = LEF(j1)
       if(icase.eq.7.and.kpol.gt.mlz) Cycle

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
        itype = Idef_itype(j1,j2,j3,j4,io,jo,ic,jc,kk1,kk2,kk3,kk4,ibi)
        Call Add_coef(CCC,kk1,kk2,kk3,kk4)
       End do

      End do !  over js
      End do	!  over is

      End do !  over ibuf

      jgen = 0
      Call MPI_REDUCE(igen,jgen,1,MPI_integer,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      igen = 0
      Call MPI_BCAST(jgen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(jgen.gt.-1) then
      Do jpol = -1,npol
       Do jtype = 1,ntype
        Call Gen_matrix(jtype,jpol)
       End do
      End do
      end if

      t3 =  MPI_WTIME()
      
      if(pri.gt.0)  write(pri,'(a,2i5,2f10.2,a)') &
       'jgen, bufer:',jgen,nnbuf, (t2-t1)/60,(t3-t2)/60,'  min'

      i=1; if(ncbuf.eq.maxnc) go to 1  ! go for new data from data bank

!-----------------------------------------------------------------------
! ... final generation of interaction matrix:

      Do jpol = -1,npol
       Do jtype = 1,ntype
        Call Gen_matrix(jtype,jpol)
       End do
      End do

      t2 =  MPI_WTIME()
      if(pri.gt.0)  write(pri,'(a,f10.2,a)') &
       'final:',(t2-t1)/60,'  min'


      End Subroutine State_res


!======================================================================
      Subroutine Read_buffer_mpi  
!======================================================================

      Use mpi
      Use bsr_mat
      USE cmdata

      Implicit none
 
      Integer :: i
      Real(8) :: t1,t2,t3,t4

      if(myid.eq.0) then
      i=1
    1 read(nub,end=2) Cbuf(i),ijtb(i),intb(i),idfb(i)
      i=i+1; if(i.le.maxnc) go to 1
    2 ncbuf=i-1
      end if

      Call MPI_BCAST(ncbuf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(CBUF,ncbuf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ijtb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(intb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(idfb,ncbuf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine Read_buffer_mpi  
