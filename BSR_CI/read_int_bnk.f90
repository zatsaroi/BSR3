!======================================================================
      Subroutine Read_int_bnk(icase)
!======================================================================
!     read angular coefficients form int_bnk file
!     and put them in rk4_data module in ordered and packed form;
!     create list of integrals (in module dbsr_mchf).
!----------------------------------------------------------------------
      Use bsr_ci;  Use c_data;  Use dets

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: Cbuf(:)
      Integer, allocatable :: ijtb(:),intb(:),idfb(:)

      Allocate(Cbuf(maxnc),ijtb(maxnc),intb(maxnc),idfb(maxnc))
!----------------------------------------------------------------------
! ... read overlap factors:

      Call Read_dets(nub)

! ... processing the data:

    1 ncbuf = 0
      Do i = 1,maxnc
       read(nub,end=2) Cbuf(i),ijtb(i),intb(i),idfb(i) 
       ncbuf = ncbuf + 1
      End do
    2 Continue

      Do ibuf=1,ncbuf; int=intb(ibuf)

       Call Decode_int (jcase,kpol,I1,I2,I3,I4,int)

       if(jcase.ne.icase) Cycle
       if(kpol.gt.km) Cycle

! ... determine the range of states for given coeff. 

       ijt= ijtb(ibuf)
       it = ijt/ibc;   jt = mod(ijt,ibc)
       is1 = IT_state1(it); js1 = IT_state1(jt)
       if(is1.eq.0.or.js1.eq.0) Cycle
       is2 = IT_state2(it); js2 = IT_state2(jt)
       if(is2.eq.0.or.js2.eq.0) Cycle
    
       idf = idfb(ibuf); nd=0; ip=0; NP=0
       if(idf.gt.0) then
        nd=KPF(idf); ip=IPF(idf); NP(1:nd)=NPF(ip+1:ip+nd)
       end if
    
!----------------------------------------------------------------------
! ... cycle over CSF:

       Do ik=is1,is2; ic=IP_stat(ik); ip1=IP_state(ic)
     
       Do jk=js1,js2; jc=IP_stat(jk); ip2=IP_state(jc)

! ... consider only low-half of interaction matrix:

       if(it.eq.jt.and.ic.lt.jc) Cycle  

       ih = max0(ic,jc)
       jh = min0(ic,jc)

       C = Cbuf(ibuf)

!----------------------------------------------------------------------
! ... check the overlap factor assuming all orbitals are orthogonal:

       DF = 1.d0
       if(idf.gt.0) then
        Do ii=1,nd
         id=NP(ii)/idef; kd=KPD(id); ip=IPD(id)
         Do jj=1,kd
          k=NPD(jj+ip)
          jp1=k/idet +    ip1; NP1(jj)=IP_orb(jp1)
          jp2=mod(k,idet)+ip2; NP2(jj)=IP_orb(jp2)
         End do
         DF = DF * VDET(kd,NP1,NP2)**mod(NP(ii),idef)
         if(abs(DF).lt.eps_det) Exit
        End do
        if(abs(DF).lt.eps_det) Cycle
       end if 

       CC = C * DF; if(abs(CC).lt.eps_c) Cycle

!----------------------------------------------------------------------
! ... find integral:

       j1=IP_orb(i1+ip1); j2=IP_orb(i2+ip1)
       j3=IP_orb(i3+ip2); j4=IP_orb(i4+ip2)

       Select case(icase)
       Case(11) 
        HS(ih,jh) = HS(ih,jh) + CC
       Case(6,7)
        ii=min(j1,j3); jj=max(j1,j3)
        Call Add_coef(CC,0,ii,jj,ih,jh,1)
       Case(5,3,4,8,9,10)
        Call Jsym_int(5,j1,j2,j3,j4)
        ii = j1*ibi+j3; jj = j2*ibi+j4
        Call Add_coef(CC,kpol,ii,jj,ih,jh,1)
       End Select

      End do   ! over jc
      End do   ! over ic

      End do   ! over icoef

      if(ncbuf.eq.maxnc) go to 1  ! go for new data from data bank

      Deallocate(Cbuf,ijtb,intb,idfb)

      End Subroutine Read_int_bnk



!======================================================================
      Subroutine Jsym_int(met,i1,i2,i3,i4)
!======================================================================
!     use integral symmetry to obtaine the 'canonical' form
!----------------------------------------------------------------------

      m=1; ii=min(i1,i2,i3,i4)

      Select case(met)
       Case(4,5)                              !   M-,R-integrals
        if(ii.eq.i4) m = 4
        if(ii.eq.i3) m = 3
        if(ii.eq.i2) m = 2
        if(ii.eq.i1) m = 1
       Case(3)                                !   T - integrals
        if(ii.eq.i2) m = 2
        if(ii.eq.i1) m = 1
       Case(6,7,8,10)                         !   L, Z, N integrals
        if(ii.eq.i3) m = 3
        if(ii.eq.i1) m = 1
      End select

      if(m.eq.1) return

      j1 = i1; j2 = i2; j3 = i3; j4 = i4

      if(m.eq.2) then
        i1 = j2; i2 = j1; i3 = j4; i4 = j3
      elseif(m.eq.3) then
        i1 = j3; i2 = j4; i3 = j1; i4 = j2
      elseif(m.eq.4) then
        i1 = j4; i2 = j3; i3 = j2; i4 = j1
      end if

      End Subroutine Jsym_int


!======================================================================
      Real(8) Function VDET (kd,N1,N2)
!======================================================================
!     calculate the value of overlap determinant for given orbitals
!     Calls:  DET
!----------------------------------------------------------------------
      Use bsr_ci,    only: OBS
      Use conf_ls,   only: msh
      
      Implicit none
      Integer, intent(in) :: kd,N1(kd),N2(kd)
      Real(8) :: ADET(msh*msh)
      Real(8), external :: DET
      Integer :: i,j

      if(kd.eq.0) then                
       VDET = 1.d0
      elseif(kd.eq.1) then
       VDET = OBS(N1(1),N2(1))
      elseif(kd.eq.2) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))
      elseif(kd.eq.3) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2))*OBS(N1(3),N2(3)) +  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(3))*OBS(N1(3),N2(1)) +  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(1))*OBS(N1(3),N2(2)) -  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(2))*OBS(N1(3),N2(1)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))*OBS(N1(3),N2(3)) -  &
              OBS(N1(1),N2(1))*OBS(N1(2),N2(3))*OBS(N1(3),N2(2)) 
      else                
       Do i=1,kd;  Do j=1,kd
        adet((i-1)*kd+j)=OBS(N1(i),N2(j))
       End do; End do
       VDET = DET(kd,adet)      
      end if

      End Function VDET
