!======================================================================
      Subroutine Coef_ee_1conf(no,ln,iq,LS,kmax,coefs)
!======================================================================
!     compute the angular coefficients for 1 atomic states;
!----------------------------------------------------------------------
      Implicit none 

! ... input-output:

      Integer, intent(in) :: no,ln(no),iq(no),LS(5,no),kmax
      Real(8), intent(out):: coefs(no,no,0:kmax)

! ... determinant expansion:

      Integer :: kdt
      Real(8), allocatable :: Cdet(:)
      Integer, allocatable :: MLdet(:,:), MSdet(:,:) 
      Real(8) :: CC_det

! ... local variables:

      Integer :: ne, kd1, kd2, i, k, mdt, ip1(no), ip2(no)
      Integer :: in(no),md(no),nd(no),it(no)
      Integer, external :: Iglq, Iterm_LS

! ... initialize arrays:

      ne = SUM(iq(1:no))

      ip1(1)=1; ip2(1)=iq(1)
      Do i=2,no
       ip1(i)=ip2(i-1)+1
       ip2(i)=ip2(i-1)+iq(i)
      End do 

! ... determinant expansion:

      k=1; mdt=1 
      Do i=1,no
       in(i)=k; k=k+iq(i); md(i)=Iglq(ln(i),iq(i));  mdt=mdt*md(i) 
       it(i)=Iterm_LS(ln(i),iq(i),0,LS(1,i),LS(2,i),LS(3,i)) 
      End do
      if(allocated(Cdet) ) Deallocate(Cdet );  Allocate(Cdet(mdt)    )
      if(allocated(MLdet)) Deallocate(MLdet);  Allocate(MLdet(ne,mdt))
      if(allocated(MSdet)) Deallocate(MSdet);  Allocate(MSdet(ne,mdt))

      Call Det_expn_1conf  

! ... calculations:     

      coefs = 0.d0
      Do kd1 = 1,kdt
      Do kd2 = kd1,kdt
        CC_det = Cdet(kd1) * Cdet(kd2)
        if(kd1.ne.kd2) CC_det = CC_det + CC_det 
        Call me_det
      End do;  End do

CONTAINS

!======================================================================
      Subroutine Det_expn_1conf  
!======================================================================
!     determined all possible determinants and their coefficients
!     for configuration given in calling routine
!---------------------------------------------------------------------
      Implicit none
      Integer :: idet(ne)
      Integer :: ML(no),MS(no),MLp(no),MSp(no)  
      Integer :: i,j,m,ii, MLT, MST
      Real(8) :: C
      Integer, external :: ML_id, MS_id
      Real(8), external :: Clebsh, DETC_sh

! ... auxiliary arrays:
! ... in(i)  -  pointer to given shell
! ... md(i)  -  number of determinant for given shell
! ... it(i)  -  pointer to shell term
! ... mdt    -  max. possible number of determinants


      MLT = LS(4,no)
      MST = LS(5,no)

! ... exhasting all shell determinants:

      kdt=0; i=1; nd(i)=1              
    1 ii = in(i)
      Call DET_sh(ln(i),iq(i),nd(i),ML(i),MS(i),Idet(ii)) 

      m = iabs(ML(i)); if(ML(i).lt.0) m=m+2
      if(m.gt.LS(2,i)) go to 2
      m = iabs(MS(i)); if(MS(i).lt.0) m=m+2
      if(m.gt.LS(3,i)) go to 2

      if(i.eq.1) then
       MLp(1)=ML(1); MSp(1)=MS(1)
      else
       MLp(i) = MLp(i-1)+ML(i)-1
       MSp(i) = MSp(i-1)+MS(i)-1
       m = iabs(MLp(i)); if(MLp(i).lt.0) m = m + 2
       if(m.gt.LS(4,i)) go to 2
       m = iabs(MSp(i)); if(MSp(i).lt.0) m = m + 2
       if(m.gt.LS(5,i)) go to 2
      end if

      if(i.lt.no) then; i = i + 1; nd(i) = 1;  go to 1; end if

      if(MLp(no).ne.MLT) go to 2
      if(MSp(no).ne.MST) go to 2

! ... coefficient calculation:

      C = 1.d0
      Do j=1,no
       C=C*DETC_sh(ln(j),iq(j),it(j),nd(j))
       if(C.eq.0.d0) Exit
      End do
      Do j=2,no
        C=C*Clebsh(LS(4,j-1),MLp(j-1), &
                   LS(2,j  ),ML (j  ), &
                   LS(4,j  ),MLp(j  ))
        if(C.eq.0.d0) Exit
        C=C*Clebsh(LS(5,j-1),MSp(j-1), &
                   LS(3,j  ),MS (j  ), &
                   LS(5,j  ),MSp(j  ))
        if(C.eq.0.d0) Exit
      End do
      if(C.eq.0.d0) go to 2

      kdt=kdt+1;  if(kdt.gt.mdt) Stop 'kdt > mdt'
      Cdet(kdt) = C
      Do j = 1,ne;  
       MLdet(j,kdt)=(ML_id(idet(j))-1)/2
       MSdet(j,kdt)= MS_id(idet(j))
      End do

! ... selecting the next case

    2 nd(i)=nd(i)+1                
      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          
       i=i-1; go to 2
      end if
      go to 1
    3 Continue

      End Subroutine Det_expn_1conf


!======================================================================
      Subroutine me_det
!======================================================================
!     find the possible interaction orbitals in two determinants and
!     call the subroutines to calculate m.e. between possible 
!     combinations of nj-orbitals (me_ee)
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,i1,i2, j, m, idif, jdif, is,js, &
                 ii(2),jj(2), io(2), jo(2)

!----------------------------------------------------------------------
! ... the same determinants:

      if(kd1.eq.kd2) then
       Do i=1,no;  Do is = ip1(i),ip2(i)
       Do j=i,no;  Do js = ip1(j),ip2(j)
        if(js.gt.is) Call me_ee(i,j,is,js,is,js)
       End do; End do
       End do; End do
       Return
      end if

!----------------------------------------------------------------------
! ... find interaction orbitals:
 
       idif=0
       Do i=1,no
        Do i1=ip1(i),ip2(i)
         m = 0
         Do i2=ip1(i),ip2(i)
          if(MLdet(i1,kd1).ne.MLdet(i2,kd2)) Cycle
          if(MSdet(i1,kd1).ne.MSdet(i2,kd2)) Cycle
          m=1; Exit 
         End do
         if(m.eq.1) Cycle  
         idif=idif+1
         if(idif.gt.2) Return
         ii(idif)=i1; io(idif)=i
        End do
       End do
       if(idif.ne.2) Stop 'me_det: idif <> 2'
 
       jdif=0
       Do i=1,no
        Do i2=ip1(i),ip2(i)
         m = 0
         Do i1=ip1(i),ip2(i)
          if(MLdet(i1,kd1).ne.MLdet(i2,kd2)) Cycle
          if(MSdet(i1,kd1).ne.MSdet(i2,kd2)) Cycle
          m=1; Exit 
         End do
         if(m.eq.1) Cycle  
         jdif=jdif+1
         if(jdif.gt.2) Return
         jj(jdif)=i2; jo(jdif)=i
        End do
       End do
       if(jdif.ne.2) Stop 'me_det: jdif <> 2'

       if(io(1).ne.jo(1))  Stop 'me_det: io(1) <> jo(1)'
       if(io(2).ne.jo(2))  Stop 'me_det: io(2) <> jo(2)'
       Call me_ee(io(1),io(2),ii(1),ii(2),jj(1),jj(2))

       End Subroutine me_det

!======================================================================
      SUBROUTINE me_ee (i,j,i1,j1,i2,j2)
!======================================================================
!     angular part of matrix elements in nlms-representation
!     for two-electron operator
!     Calls: Check_boef
!----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, intent(in) :: i,j,i1,i2,j1,j2
      Integer :: k,kz,ib,ii,jj

      ii = max(i,j);  jj = min(i,j)

      if(MLdet(i1,kd1) + MLdet(j1,kd1).ne. &
         MLdet(i2,kd2) + MLdet(j2,kd2)) Return
      if(MSdet(i1,kd1) + MSdet(j1,kd1).ne. &
         MSdet(i2,kd2) + MSdet(j2,kd2)) Return

      kz = (-1)**(i1+i2+j1+j2)

      Call Check_boef(ln(i),MLdet(i1,kd1),MSdet(i1,kd1), & 
                      ln(j),MLdet(j1,kd1),MSdet(j1,kd1), &
                      ln(i),MLdet(i2,kd2),MSdet(i2,kd2), &
                      ln(j),MLdet(j2,kd2),MSdet(j2,kd2))

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib); if(k.gt.kmax) Cycle
       coefs(ii,jj,k) = coefs(ii,jj,k) + Boef(ib)*CC_det*kz
      End do

      kz = -kz

      Call Check_boef(ln(i),MLdet(i1,kd1),MSdet(i1,kd1), & 
                      ln(j),MLdet(j1,kd1),MSdet(j1,kd1), &
                      ln(j),MLdet(j2,kd2),MSdet(j2,kd2), &
                      ln(i),MLdet(i2,kd2),MSdet(i2,kd2))

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib); if(k.gt.kmax) Cycle
       coefs(jj,ii,k) = coefs(jj,ii,k) + Boef(ib)*CC_det*kz
      End do

      End Subroutine me_ee

      End Subroutine Coef_ee_1conf

