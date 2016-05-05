!======================================================================
      Subroutine D_matr
!======================================================================
!     This routine computes the matrix of dipole matrix elements
!     in length and velocity forms according to MULT_BNK (unit nub)
!----------------------------------------------------------------------
      Use bsr_dmat
      Use conf_LS; Use term_LS
      Use dets;    Use new_dets; Use new_defs
      Use cmdata,          only: itype,ntype,ncdata
      Use spline_orbitals, only: iech

      Implicit none
      Integer :: i,j,i1,i2,j1,j2,m,ii,kk1,kk2,kk3,io,jo
      Integer :: ijt,int,idf,int_type,kz,ip1,ip2,ich1,ich2
      Integer :: it,jt,ic,jc,ik,jk,is,js,is1,is2,js1,js2,ih,jh
      Real(8) :: C,CC, CCL,CCV, CCCL,CCCV, CL,CV
      Integer, external :: Idef_dtype, no_ic_LS, Ifind_pert

!----------------------------------------------------------------------
! ... prepare main arrays:

      Call D_prep

!----------------------------------------------------------------------
! ... prepare module cmdata:

      Call Allocate_cmdata(1)

!----------------------------------------------------------------------
! ... read determinant factors:

      Call Read_dets(nub)

!----------------------------------------------------------------------
! ... processing the data:

    1 read(nub,end=2) C, ijt,int,idf

      it=ijt/ibc;  jt=mod(ijt,ibc)

      is1 = IT_state1(it); js1 = IT_state1(jt)
      is2 = IT_state2(it); js2 = IT_state2(jt)

      if(is1.eq.0.or.js1.eq.0) go to 1   ! no relevant states

      Call Decode_mult(int_type,i1,i2,int)

!----------------------------------------------------------------------
! ... cycle over states:

      Do ik=is1,is2;  is=IP_stat(ik);  ip1=IP_state(is)
       no1=no_ic_LS (is); np1(1:no1)=IP_orb(ip1+1:ip1+no1)
      Do jk=js1,js2;  js=IP_stat(jk);  ip2=IP_state(js)
       no2=no_ic_LS (js); np2(1:no2)=IP_orb(ip2+1:ip2+no2)

      if(is.gt.ncfg1.and.js.gt.ncfg1) Cycle
      if(is.le.ncfg1.and.js.le.ncfg1) Cycle

      if(it.eq.jt.and.is.gt.js) Cycle 

! ... check if conf.s are in the right order:

      if(is.le.ncfg1) then
       ih = is;  jh = js-ncfg1
       i=LSP1(ih); j=LSP2(jh); m=0        
       kz=1
      else
       ih = js;  jh = is-ncfg1
       i=LSP1(ih); j=LSP2(jh); m=1        
       kz = (ILterm1(i)-ILterm2(j)+ISterm1(i)-ISterm2(j))/ 2        
       kz = (-1)**kz        
      end if

! ... include the expansion coefficients: 

      CC = C*C1(ih)*C2(jh)*kz;  if(abs(CC).lt.Eps_C) Cycle

! ... define integral:

      if(m.eq.0) then
       j1=IP_orb(ip1+i1); j2=IP_orb(ip2+i2)
      else
       j1=IP_orb(ip2+i2); j2=IP_orb(ip1+i1)
      end if       
      ich1=iech(j1); ich2=iech(j2)

! ... J-dependence factor:

      CCL=CC; CCV=CC
      if(jmode.eq.2) then
       CCL = CC * DJ(i,j); CCV = CC * DJM(i,j)
	      if(abs(CCL)+abs(CCV).lt.Eps_C) Cycle
      end if 
      if(int_type.eq.2) CCV = 0.d0
      if(int_type.eq.3) CCL = 0.d0

! ... find index of the pertuber if any: 

      ic = Ifind_pert(ilsp1,ih)
      jc = Ifind_pert(ilsp2,jh)                

! ... find overlap factors with extracted continuum:  

      Call Det_fact(idf,np1,np2); if(nndef.eq.0) Cycle

! ... send the final coefficients to archive:

      Do i = 1,nndef
       CCCL = CCL * Adef(i)
       CCCV = CCV * Adef(i)

       io=iof(i); jo=jof(i)

       if(jo.gt.0.and.iech(io/ibo).eq.0) then
        ii = io; io=jo; jo=ii
       end if

       itype=Idef_dtype(j1,j2,ich1,ich2,io,jo,ic,jc,kk1,kk2,kk3)
       CL=CCCL; CV=CCCV
       if(ich1+ich2.eq.0) then 
        CL=CCCL*dipL(j1,j2); CV=CCCV*dipV(j1,j2)
       end if
       Call Add_coef(CL,CV,kk1,kk2,kk3)
      End do
!-----------------------------------------------------------------------
      End do    !  over js
      End do	   !  over is
!-----------------------------------------------------------------------
      go to 1   !  new data from data bank
    2 Continue
!-----------------------------------------------------------------------
! ... final generation of interaction matrix:

      Do itype = 1,ntype; Call Gen_matrix;  End do

      End Subroutine D_MATR


