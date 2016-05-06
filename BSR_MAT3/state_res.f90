!======================================================================
      Subroutine State_res
!======================================================================
!     extracts the data from int_bnk for specific case
!----------------------------------------------------------------------
      Use bsr_mat
      Use conf_LS; Use orb_LS 
      Use channel; Use target
      Use cmdata
      Use dets; Use new_dets; Use new_defs
      Use z_core, only: mlz 

      Implicit none
      Integer :: i,j,k,it,jt,is,js,is1,is2,js1,js2,ik,jk,ich,jch,ip1,ip2
      Integer :: nccoef,ncbuf,ibuf,jcase,ijt,jpol,jtype
      Integer :: i1,i2,i3,i4,j1,j2,j3,j4,kk1,kk2,kk3,kk4,ic,jc,io,jo,idf
      Integer :: ILT1,ILT2,IST1,IST2
      Real(8) :: C,CC,CCC
      Integer, external :: Ifind_channel, Idef_itype, no_ic_LS
      Real(8), external :: Z_6j

! ... initilize blocks in cmdata module:

      Call Initilize_cmdata

! ... read det. factors (to get access to coefficients):

      Call Read_dets
      if(icase.eq.11) then
       write(pri,'(/a,2i8)') 'ndet  = ',ndet,kdet
       write(pri,'( a,2i8)') 'ndef  = ',ndef,kdef
      end if

      nccoef = 0
!----------------------------------------------------------------------
! ... read coef.s from int.bnk to buffer:

      i=1
    1 read(nub,end=2) Cbuf(i),ijtb(i),intb(i),idfb(i)
      i=i+1; if(i.le.maxnc) go to 1
    2 ncbuf=i-1

!----------------------------------------------------------------------
! ... look through the buffer:

      Do ibuf=1,ncbuf

      nccoef=nccoef+1
      Call Decode_INT (jcase,k,i1,i2,i3,i4,intb(ibuf)) 
      
      if(icase.ne.jcase) Cycle

      kpol = k 
      Select case(icase)
       Case(4,8,9,10); kpol=k-1
       Case(6,7,11);   kpol=0
      End select
      if(kpol.gt.mk) Cycle

! ... determine the range of states for given coeff. 

      ijt=ijtb(ibuf)
      it = ijt/ibc;   jt = mod(ijt,ibc)
      is1 = IT_state1(it); js1 = IT_state1(jt)
      if(is1.eq.0.or.js1.eq.0) Cycle
      is2 = IT_state2(it); js2 = IT_state2(jt)

      idf = idfb(ibuf)
      C = CBUF(ibuf)
!----------------------------------------------------------------------
!                                        loop over all relevant states: 

      Do ik=is1,is2; is=IP_stat(ik); ip1=IP_state(is)
       no1=no_ic_LS (is); np1(1:no1)=IP_orb(ip1+1:ip1+no1)
       ich=Ifind_channel(is)
       Call Term_ic (is,ILT1,IST1)
      Do jk=js1,js2; js=IP_stat(jk); ip2=IP_state(js)
       no2=no_ic_LS (js); np2(1:no2)=IP_orb(ip2+1:ip2+no2)
       jch=Ifind_channel(js)
       Call Term_ic (js,ILT2,IST2)

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

! ... we do not need anymore the configuration index, 
! ... only pertuber index if any:   

       i=0; if(ich.gt.nch) i=ich-nch
       j=0; if(jch.gt.nch) j=jch-nch
       ic = max(i,j)
       jc = min(i,j)

! ... define integral for specific orbitals

       if(icase.eq.11) then
        j1 = 1; j2 = 1; j3 = 1; j4 = 1
       else
        j1=np1(i1); j2=np1(i2); j3=np2(i3); j4=np2(i4)
       end if

! ... one-electron integrals:

       if(icase.eq.6.or.icase.eq.7) kpol = LEF(j1)
       if(icase.eq.7.and.kpol.gt.mlz) Cycle

! ... J-dependence for relativistic ccorrections: 
       
       if(icase.gt.6.and.icase.lt.11) then
        k = 3; if(icase.eq.10) k = 5
        CC = CC * (-1)**((ILT2+IST1+jpar-3)/2)*   &
             Z_6j(ILT1,IST1,jpar,IST2,ILT2,k)
        if(abs(CC).lt.Eps_C) Cycle
       end if

! ... find canonical form:

       Call Jsym_int(icase,j1,j2,j3,j4)

! ... find overlap factors with extracted continuum:  

       Call Det_fact(idf,np1,np2); if(nndef.eq.0) Cycle

! ... send the final coefficients to archive:

       Do i = 1,nndef
        CCC = CC * Adef(i); io=iof(i); jo=jof(i)
        itype = Idef_itype(j1,j2,j3,j4,io,jo,ic,jc,kk1,kk2,kk3,kk4,ibi)
        Call Add_coef(CCC,kk1,kk2,kk3,kk4)
       End do

      End do   !  over js
      End do   !  over is

      End do    !  over ibuf

      i=1; if(ncbuf.eq.maxnc) go to 1  ! go for new data from data bank

!-----------------------------------------------------------------------
! ... final generation of interaction matrix:

      Do jpol = -1,npol
       Do jtype = 1,ntype
        Call Gen_matrix(jtype,jpol)
       End do
      End do

      End Subroutine State_res
