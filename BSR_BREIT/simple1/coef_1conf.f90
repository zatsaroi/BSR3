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
      Integer :: kd,i,j,k,m,ii, MLT, MST
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
      Integer :: k,kz,ib

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
       coefs(i,j,k) = coefs(i,j,k) + Boef(ib)*CC_det*kz
      End do

      kz = -kz

      Call Check_boef(ln(i),MLdet(i1,kd1),MSdet(i1,kd1), & 
                      ln(j),MLdet(j1,kd1),MSdet(j1,kd1), &
                      ln(j),MLdet(j2,kd2),MSdet(j2,kd2), &
                      ln(i),MLdet(i2,kd2),MSdet(i2,kd2))

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib); if(k.gt.kmax) Cycle
       coefs(j,i,k) = coefs(j,i,k) + Boef(ib)*CC_det*kz
      End do

      End Subroutine me_ee

      End Subroutine Coef_ee_1conf



!======================================================================
      Integer Function Iterm_LS (l,iq,k,IA,IL,IS)
!======================================================================
!
!     provides information about shell terms:
!
!     k>0  --> IA,IL,IS = k-th term of lq-subshell
!     k=0  --> Iterm_LS = position of the IA,IL,IS term in subshell list
!     k=-1 --> Iterm_LS = number of terms in lq-subshell
!
!----------------------------------------------------------------------

      Implicit none

      Integer :: l,iq,k,IA,IL,IS,   kl,jq,i,j,n

      Integer, parameter :: nsh=10, nterms=419
      Integer :: lqnp(4,nsh)
      Data lqnp / 0, 0,   1,   0,    &    ! full 
		                1, 3,   3,   1,    &    ! p3
                  2, 3,   8,   4,    &    ! d3
			            	  2, 4,  16,  12,    &    ! d4
				              2, 5,  16,  28,    &    ! d5
              		  3, 3,  17,  44,    &    ! f3
			            	  3, 4,  47,  61,    &    ! f4
				              3, 5,  73, 108,    &    ! f5
				              3, 6, 119, 181,    &    ! f6
                  3, 7, 119, 300/         ! f7

      Integer :: ILS(3,nterms)
      Data ILS/                                              &
         0,1,1,                                              &! full
         3,1,4,1,3,2,3,5,2,                                  &! p3
         3,3,2,3,3,4,1,5,2,3,5,2,3,7,2,3,7,4,3, 9,2,3,11,2,  &! d3
         0,1,1,4,1,1,2,3,3,4,3,3,2,5,1,4,5,1,4, 5,3,4, 5,5,  &! d4
         4,7,1,2,7,3,4,7,3,2,9,1,4,9,1,4,9,3,4,11,3,4,13,1,  &
         5,1,2,5,1,6,3,3,2,3,3,4,1,5,2,3,5,2,5, 5,2,5, 5,4,  &! d5
         3,7,2,5,7,2,3,7,4,3,9,2,5,9,2,5,9,4,3,11,2,5,13,2,  &
         3, 1, 4, 3, 5, 4, 3, 7, 4, 3, 9, 4, 3,13, 4, &       ! f3 
         3, 3, 2, 1, 5, 2, 2, 5, 2, 1, 7, 2, 2, 7, 2, &       
         1, 9, 2, 2, 9, 2, 1,11, 2, 2,11, 2, 3,13, 2, &       
         3,15, 2, 3,17, 2,                            &       
         4, 1, 5, 4, 5, 5, 4, 7, 5, 4, 9, 5, 4,13, 5, &       ! f4 
         1, 3, 3, 2, 3, 3, 3, 3, 3, 1, 5, 3, 2, 5, 3, &       
         1, 7, 3, 2, 7, 3, 3, 7, 3, 4, 7, 3, 1, 9, 3, &       
         2, 9, 3, 3, 9, 3, 1,11, 3, 2,11, 3, 3,11, 3, &       
         4,11, 3, 1,13, 3, 2,13, 3, 1,15, 3, 2,15, 3, &       
         4,17, 3, 4,19, 3, 1, 1, 1, 2, 1, 1, 1, 5, 1, &       
         2, 5, 1, 3, 5, 1, 4, 5, 1, 4, 7, 1, 1, 9, 1, &       
         2, 9, 1, 3, 9, 1, 4, 9, 1, 1,11, 1, 2,11, 1, &       
         1,13, 1, 2,13, 1, 3,13, 1, 4,15, 1, 1,17, 1, &       
         2,17, 1, 4,21, 1,                            &       
         5, 3, 6, 5, 7, 6, 5,11, 6, 3, 1, 4, 1, 3, 4, &       ! f5 
         2, 3, 4, 1, 5, 4, 2, 5, 4, 3, 5, 4, 1, 7, 4, &       
         2, 7, 4, 3, 7, 4, 4, 7, 4, 1, 9, 4, 2, 9, 4, &       
         3, 9, 4, 4, 9, 4, 1,11, 4, 2,11, 4, 3,11, 4, &       
         1,13, 4, 2,13, 4, 3,13, 4, 1,15, 4, 2,15, 4, &       
         5,17, 4, 5,19, 4, 1, 3, 2, 2, 3, 2, 3, 3, 2, &       
         4, 3, 2, 1, 5, 2, 2, 5, 2, 3, 5, 2, 4, 5, 2, &       
         5, 5, 2, 1, 7, 2, 2, 7, 2, 3, 7, 2, 4, 7, 2, &       
         5, 7, 2, 6, 7, 2, 7, 7, 2, 1, 9, 2, 2, 9, 2, &       
         3, 9, 2, 4, 9, 2, 5, 9, 2, 6, 9, 2, 1,11, 2, &       
         2,11, 2, 3,11, 2, 4,11, 2, 5,11, 2, 6,11, 2, &       
         7,11, 2, 1,13, 2, 2,13, 2, 3,13, 2, 4,13, 2, &       
         5,13, 2, 1,15, 2, 2,15, 2, 3,15, 2, 4,15, 2, &       
         5,15, 2, 1,17, 2, 2,17, 2, 3,17, 2, 1,19, 2, &       
         2,19, 2, 5,21, 2, 5,23, 2,                   &       
         6, 7, 7, 4, 1, 5, 6, 3, 5, 1, 5, 5, 2, 5, 5, &       ! f6 
         3, 5, 5, 1, 7, 5, 2, 7, 5, 1, 9, 5, 2, 9, 5, &       
         3, 9, 5, 1,11, 5, 2,11, 5, 1,13, 5, 2,13, 5, &       
         6,15, 5, 6,17, 5, 1, 3, 3, 2, 3, 3, 3, 3, 3, &       
         4, 3, 3, 5, 3, 3, 6, 3, 3, 1, 5, 3, 2, 5, 3, &       
         3, 5, 3, 4, 5, 3, 5, 5, 3, 1, 7, 3, 2, 7, 3, &       
         3, 7, 3, 4, 7, 3, 5, 7, 3, 6, 7, 3, 7, 7, 3, &       
         8, 7, 3, 9, 7, 3, 1, 9, 3, 2, 9, 3, 3, 9, 3, &       
         4, 9, 3, 5, 9, 3, 6, 9, 3, 7, 9, 3, 1,11, 3, &       
         2,11, 3, 3,11, 3, 4,11, 3, 5,11, 3, 6,11, 3, &       
         7,11, 3, 8,11, 3, 9,11, 3, 1,13, 3, 2,13, 3, &       
         3,13, 3, 4,13, 3, 5,13, 3, 6,13, 3, 1,15, 3, &       
         2,15, 3, 3,15, 3, 4,15, 3, 5,15, 3, 6,15, 3, &       
         1,17, 3, 2,17, 3, 3,17, 3, 1,19, 3, 2,19, 3, &       
         3,19, 3, 6,21, 3, 6,23, 3, 1, 1, 1, 2, 1, 1, &       
         3, 1, 1, 4, 1, 1, 6, 3, 1, 1, 5, 1, 2, 5, 1, &       
         3, 5, 1, 4, 5, 1, 5, 5, 1, 6, 5, 1, 1, 7, 1, &       
         2, 7, 1, 3, 7, 1, 4, 7, 1, 1, 9, 1, 2, 9, 1, &       
         3, 9, 1, 4, 9, 1, 5, 9, 1, 6, 9, 1, 7, 9, 1, &       
         8, 9, 1, 1,11, 1, 2,11, 1, 3,11, 1, 4,11, 1, &       
         1,13, 1, 2,13, 1, 3,13, 1, 4,13, 1, 5,13, 1, &       
         6,13, 1, 7,13, 1, 1,15, 1, 2,15, 1, 3,15, 1, &       
         1,17, 1, 2,17, 1, 3,17, 1, 4,17, 1, 1,19, 1, &       
         2,19, 1, 1,21, 1, 2,21, 1, 6,25, 1,          &       
         7, 1, 8, 5, 3, 6, 7, 5, 6, 5, 7, 6, 7, 9, 6, &       ! f7 
         5,11, 6, 7,13, 6, 1, 1, 4, 2, 1, 4, 1, 3, 4, &       
         2, 3, 4, 1, 5, 4, 2, 5, 4, 3, 5, 4, 4, 5, 4, &       
         5, 5, 4, 6, 5, 4, 1, 7, 4, 2, 7, 4, 3, 7, 4, &       
         4, 7, 4, 5, 7, 4, 1, 9, 4, 2, 9, 4, 3, 9, 4, &       
         4, 9, 4, 5, 9, 4, 6, 9, 4, 7, 9, 4, 1,11, 4, &       
         2,11, 4, 3,11, 4, 4,11, 4, 5,11, 4, 1,13, 4, &       
         2,13, 4, 3,13, 4, 4,13, 4, 5,13, 4, 1,15, 4, &       
         2,15, 4, 3,15, 4, 1,17, 4, 2,17, 4, 3,17, 4, &       
         5,19, 4, 7,21, 4, 1, 1, 2, 2, 1, 2, 1, 3, 2, &       
         2, 3, 2, 3, 3, 2, 4, 3, 2, 5, 3, 2, 1, 5, 2, &       
         2, 5, 2, 3, 5, 2, 4, 5, 2, 5, 5, 2, 6, 5, 2, &       
         7, 5, 2, 1, 7, 2, 2, 7, 2, 3, 7, 2, 4, 7, 2, &       
         5, 7, 2, 6, 7, 2, 7, 7, 2, 8, 7, 2, 9, 7, 2, &       
         0, 7, 2, 1, 9, 2, 2, 9, 2, 3, 9, 2, 4, 9, 2, &       
         5, 9, 2, 6, 9, 2, 7, 9, 2, 8, 9, 2, 9, 9, 2, &       
         0, 9, 2, 1,11, 2, 2,11, 2, 3,11, 2, 4,11, 2, &       
         5,11, 2, 6,11, 2, 7,11, 2, 8,11, 2, 9,11, 2, &       
         1,13, 2, 2,13, 2, 3,13, 2, 4,13, 2, 5,13, 2, &       
         6,13, 2, 7,13, 2, 8,13, 2, 9,13, 2, 1,15, 2, &       
         2,15, 2, 3,15, 2, 4,15, 2, 5,15, 2, 6,15, 2, &       
         7,15, 2, 1,17, 2, 2,17, 2, 3,17, 2, 4,17, 2, &       
         5,17, 2, 1,19, 2, 2,19, 2, 3,19, 2, 4,19, 2, &       
         1,21, 2, 2,21, 2, 5,23, 2, 7,25, 2/                  

!    empty or full shell - 1S0
!    l1 or l[4l+1] - 2L1
!    l2 or l[4l]   - 1S0, 3P2, 1D2, 3F2, 1G2, ...
!    p3 - 4S3, 2P1, 2D3
!    d3 - 2P3, 4P3, 2D1, 2D3, 2F3, 4F3, 2G3, 2H3
!    d4 - 1S0, 1S4, 3P2, 3P4, 1D2, 1D4, 3D4, 5D4,
!         1F4, 3F2, 3F4, 1G2, 1G4, 3G4, 3H4, 1I4
!    d5 - 2S5, 6S5, 2P3, 4P3, 2D1, 2D3, 2D5, 4D5,
!       - 2F3, 2F5, 4F3, 2G3, 2G5, 4G5, 2H3, 2I5
!    f3 - 4S 3, 4D 3, 4F 3, 4G 3, 4I 3, 2P 3, 2D13, 2D23, 2F11, 2F23, 
!         2G13, 2G23, 2H13, 2H23, 2I 3, 2K 3, 2L 3 
!    f4 - 5S 4, 5D 4, 5F 4, 5G 4, 5I 4, 3P12, 3P24, 3P34, 3D14, 3D24, 
!         3F12, 3F24, 3F34, 3F44, 3G14, 3G24, 3G34, 3H12, 3H24, 3H34, 
!         3H44, 3I14, 3I24, 3K14, 3K24, 3L 4, 3M 4, 1S10, 1S24, 1D12, 
!         1D24, 1D34, 1D44, 1F 4, 1G12, 1G24, 1G34, 1G44, 1H14, 1H24, 
!         1I12, 1I24, 1I34, 1K 4, 1L14, 1L24, 1N 4 
!    f5 - 6P 5, 6F 5, 6H 5, 4S 3, 4P15, 4P25, 4D13, 4D25, 4D35, 4F13, 
!         4F25, 4F35, 4F45, 4G13, 4G25, 4G35, 4G45, 4H15, 4H25, 4H35, 
!         4I13, 4I25, 4I35, 4K15, 4K25, 4L 5, 4M 5, 2P13, 2P25, 2P35, 
!         2P45, 2D13, 2D23, 2D35, 2D45, 2D55, 2F11, 2F23, 2F35, 2F45, 
!         2F55, 2F65, 2F75, 2G13, 2G23, 2G35, 2G45, 2G55, 2G65, 2H13, 
!         2H23, 2H35, 2H45, 2H55, 2H65, 2H75, 2I13, 2I25, 2I35, 2I45, 
!         2I55, 2K13, 2K25, 2K35, 2K45, 2K55, 2L13, 2L25, 2L35, 2M15, 
!         2M25, 2N 5, 2O 5 
!    f6 - 7F 6, 5S 4, 5P 6, 5D14, 5D26, 5D36, 5F14, 5F26, 5G14, 5G26, 
!         5G36, 5H16, 5H26, 5I14, 5I26, 5K 6, 5L 6, 3P12, 3P24, 3P34, 
!         3P46, 3P56, 3P66, 3D14, 3D24, 3D36, 3D46, 3D56, 3F12, 3F24, 
!         3F34, 3F44, 3F56, 3F66, 3F76, 3F86, 3F96, 3G14, 3G24, 3G34, 
!         3G46, 3G56, 3G66, 3G76, 3H12, 3H24, 3H34, 3H44, 3H56, 3H66, 
!         3H76, 3H86, 3H96, 3I14, 3I24, 3I36, 3I46, 3I56, 3I66, 3K14, 
!         3K24, 3K36, 3K46, 3K56, 3K66, 3L14, 3L26, 3L36, 3M14, 3M26, 
!         3M36, 3N 6, 3O 6, 1S10, 1S24, 1S36, 1S46, 1P 6, 1D12, 1D24, 
!         1D34, 1D44, 1D56, 1D66, 1F14, 1F26, 1F36, 1F46, 1G12, 1G24, 
!         1G34, 1G44, 1G56, 1G66, 1G76, 1G86, 1H14, 1H24, 1H36, 1H46, 
!         1I12, 1I24, 1I34, 1I46, 1I56, 1I66, 1I76, 1K14, 1K26, 1K36, 
!         1L14, 1L24, 1L36, 1L46, 1M16, 1M26, 1N14, 1N26, 1Q 6 
!    f7 - 8S 7, 6P 5, 6D 7, 6F 5, 6G 7, 6H 5, 6I 7, 4S13, 4S27, 4P15, 
!         4P25, 4D13, 4D25, 4D35, 4D47, 4D57, 4D67, 4F13, 4F25, 4F35, 
!         4F45, 4F57, 4G13, 4G25, 4G35, 4G45, 4G57, 4G67, 4G77, 4H15, 
!         4H25, 4H35, 4H47, 4H57, 4I13, 4I25, 4I35, 4I47, 4I57, 4K15, 
!         4K25, 4K37, 4L15, 4L27, 4L37, 4M 5, 4N 7, 2S17, 2S27, 2P13, 
!         2P25, 2P35, 2P45, 2P57, 2D13, 2D23, 2D35, 2D45, 2D55, 2D67, 
!         2D77, 2F11, 2F23, 2F35, 2F45, 2F55, 2F65, 2F75, 2F87, 2F97, 
!         2F07, 2G13, 2G23, 2G35, 2G45, 2G55, 2G65, 2G77, 2G87, 2G97, 
!         2G07, 2H13, 2H23, 2H35, 2H45, 2H55, 2H65, 2H75, 2H87, 2H97, 
!         2I13, 2I25, 2I35, 2I45, 2I55, 2I67, 2I77, 2I87, 2I97, 2K13, 
!         2K25, 2K35, 2K45, 2K55, 2K67, 2K77, 2L13, 2L25, 2L35, 2L47, 
!         2L57, 2M15, 2M25, 2M37, 2M47, 2N15, 2N27, 2O 5, 2Q 7  
 
      Iterm_LS = 0
      kl=4*l+2
      jq=iq; if(iq.gt.kl/2) jq=kl-iq
!----------------------------------------------------------------------
! ... the number of terms in lq-subshell:

      if(k.eq.-1) then            
      if(l.lt.0.or.iq.lt.0.or.iq.gt.kl) then
       write(*,*) 'l,q = ', l,iq
       Stop ' Iterm_LS: not possible lq combination'
      end if
      if(jq.eq.0.or.jq.eq.1) then;   Iterm_LS=1              ! l^1
      elseif(jq.eq.2) then;          Iterm_LS=l+l+1          ! l^2
      else
       Do i=1,nsh
        if(l .ne.lqnp(1,i)) Cycle
        if(jq.ne.lqnp(2,i)) Cycle
        Iterm_LS=lqnp(3,i)
      		Exit
	      End do
      end if
      if(Iterm_LS.eq.0) then
       write(*,*)'Iterm_LS: not included case of lq'
       Stop ' '
      end if
      Return
      end if

!----------------------------------------------------------------------
! ... position of the IA,IL,IS term in list

      if(k.eq.0) then          

      if(jq.eq.0.and.IL.eq.1.and.IS.eq.1.and.IA.eq.0)     Iterm_LS=1
      if(jq.eq.1.and.IS.eq.2.and.IL.eq.l+l+1.and.IA.eq.1) Iterm_LS=1
      if(jq.eq.2.and.IL.ge.1.and.IL.le.4*l+1.and.  &
         IS.eq.2*mod((IL-1)/2,2)+1.and.&
         ((IL.eq.1.and.IA.eq.0).or.(IL.gt.1.and.IA.eq.2))) &
         Iterm_LS=(IL+1)/2
      if(Iterm_LS.gt.0) Return

      n = 0
      Do i=1,nsh
       if(l.ne.lqnp(1,i)) Cycle
       if(jq.ne.lqnp(2,i)) Cycle
       n = lqnp(3,i)
       j = lqnp(4,i)
       Exit
	     End do

      Do i=1,n
       if(IA.eq.ILS(1,i+j).and.IL.eq.ILS(2,i+j).and. &
          IS.eq.ILS(3,i+j)) Iterm_LS=i
      End do

      if(Iterm_LS.eq.0) then
       write(*,'(a,5i5)') 'lq, ALS=',l,iq,IA,IL,IS 
       Stop ' Iterm_LS: incorrect term'
      end if

      Return
      end if

!----------------------------------------------------------------------
! .. k-th term of lq-subshell:

     if(k.gt.0) then                       

      IL=0
      if(jq.eq.0.and.k.eq.1) then
       IA=0; IL=1; IS=1
      elseif(jq.eq.1.and.k.eq.1) then
       IA=1; IL=l+l+1; IS=2
      elseif(jq.eq.2.and.k.le.l+l+1) then
       IA=0; if(k.gt.1) IA=2; IL=k+k-1; IS=2*mod(k-1,2)+1
      else

       n = 0
       Do i=1,nsh
        if(l.ne.lqnp(1,i)) Cycle
        if(jq.ne.lqnp(2,i)) Cycle
        n = lqnp(3,i)
        j = lqnp(4,i)
        Exit
	      End do

       if(n.eq.0.or.k.gt.n) then
        write(*,'(a,5i5)') 'lq, ALS=',l,iq,IA,IL,IS 
        Stop ' Iterm_LS: incorrect term request'
       end if

       IA=ILS(1,k+j)
       IL=ILS(2,k+j)
       IS=ILS(3,k+j)

      end if

      Iterm_LS=k
      Return

     end if  ! over k > 0

     End Function Iterm_LS

!--------------------------------------------------------------------
      Real(8) FUNCTION Z_3j (j1,m1,j2,m2,j3,m3) 
!--------------------------------------------------------------------
!
!     determines the value of the 3j-symbols without direct using of
!     factorials. The following expression for the 3j-symbols is used:
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977)
!
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} *
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ]
!                         SUM(z) {   (-1)^z  /
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! *
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] }
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values a(i) and b(i)
!     (see below the text of program) then
!
!     3j =         (-1) ^ Sum[a(i)]
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] }
!                  Sum(z) { (-1)^z  /
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] }
!
!     (below the moments are used in (2J+1)-representation)
!
!--------------------------------------------------------------------

      Implicit none
 
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3

      Integer(4) :: i,i_max,k,kk,m,iz,iz_min,iz_max
      Real(8) :: x,y,z

      Integer(4) a(3),b(3),J(16) 

      Z_3j=0.0

      IF(M1+M2+M3-3.ne.0) RETURN    ! check of conservation rules
      J(1)= J1+J2-J3-1
      J(2)= J1-J2+J3-1
      J(3)= J2-J1+J3-1
      J(4)= J1+M1-2
      J(5)= J1-M1
      J(6)= J2-M2
      J(7)= J2+M2-2
      J(8)= J3+M3-2
      J(9)= J3-M3
      Do I=1,9
       IF(J(i).lt.0.or.mod(J(i),2).eq.1) RETURN
      End do

      a(1) = 0                         ! auxiliary values
      a(2) = (j2-j3-m1+1)/2
      a(3) = (j1-j3+m2-1)/2
      b(1) = (j1+j2-j3-1)/2
      b(2) = (j1-m1)/2
      b(3) = (j2+m2-2)/2

      IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum
      IZ_max=MIN0(b(1),b(2),b(3))
      IF(IZ_max.LT.IZ_min) Return

      Do I=1,3                         ! constant factorial parameters
      Do K=1,3
       J(I+3*K-3)=b(i)-a(k)
      End do
      End do
      J(10)=(j1+j2+j3-3)/2+1

      Do I=1,3
       J(I+10)=IZ_min-a(i)               ! initial factorial parameters
       J(I+13)=b(i)-IZ_min               ! in the sum
      End do

      Z=0.0
      DO IZ=IZ_min,IZ_max                 ! summation

       I_max=0                            ! max. factorial
       Do I=1,16
        if(J(i).gt.I_max) I_max=J(i)
       End do

       Y=1.0
       DO I=2,I_max         ! estimation of one term in sum
        K=0                 ! K - the extent of the integer I in term
        DO M=1,9
         IF(J(M).GE.I) K=K+1
        End do
        IF(J(10).GE.I) K=K-1
        DO M=11,16
         IF(J(M).GE.I) K=K-2
        End do
        IF(K.EQ.0) Cycle

        X=DBLE(I)                   ! Y = Y * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) Y=Y*X
          IF(K.LT.0) Y=Y/X
         END DO
        END IF
        IF(mod(K,2).EQ.+1) Y=Y*SQRT(X)
        IF(mod(K,2).EQ.-1) Y=Y/SQRT(X)
       End do

       IF(mod(IZ,2).eq.1) Y=-Y
       Z=Z+Y

       Do I=11,13                  ! new factorial parameters in sum
        J(I)=J(I)+1
       End do
       DO I=14,16
        J(I)=J(I)-1
       End do

      End do                       ! end of summation

      K=a(1)+a(2)+a(3)
      if(mod(k,2).ne.0) Z=-Z
      Z_3j=Z

      END FUNCTION Z_3j


!--------------------------------------------------------------------
      Real(8) FUNCTION Z_3jj(j1,m1,j2,m2,j3,m3)
!--------------------------------------------------------------------

      IMPLICIT NONE
  
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), External :: Z_3j
       
      Z_3jj=Z_3j(j1+j1+1,m1+m1+1,j2+j2+1,m2+m2+1,j3+j3+1,m3+m3+1)

      End FUNCTION Z_3jj

!--------------------------------------------------------------------
      Real(8) FUNCTION Z_3j2(j1,m1,j2,m2,j3,m3)
!--------------------------------------------------------------------

      IMPLICIT NONE
  
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), External :: Z_3j
       
      Z_3j2=Z_3j(j1+1,m1+1,j2+1,m2+1,j3+1,m3+1)

      End FUNCTION Z_3j2


!====================================================================
      Real(8) FUNCTION CLEBSH(J1,M1,J2,M2,J,M)
!====================================================================
!
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are used in (2J+1)-representation)
!
!     Call:  Z_3j
!--------------------------------------------------------------------

      Implicit None
      
      Integer, intent(in) :: J, M, J1, M1, J2, M2
      Real(8), External :: Z_3j

      Clebsh=(-1)**((j1-j2+m-1)/2)*sqrt(DBLE(J))*   &
             Z_3j(j1,m1,j2,m2,J,-m+2)

      END FUNCTION CLEBSH


!======================================================================
      Real(8) FUNCTION CLEBCH(L1,M1,L2,M2,L,M)
!======================================================================

      Implicit None
      
      Integer, intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebch = Clebsh(l1+l1+1,m1+m1+1,l2+l2+1,m2+m2+1,l+l+1,m+m+1)

      End FUNCTION CLEBCH


!======================================================================
      Real(8) FUNCTION CLEBSH2(L1,M1,L2,M2,L,M)
!======================================================================

      Implicit None
      
      Integer(4), intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebsh2 = Clebsh(l1+1,m1+1,l2+1,m2+1,l+1,m+1)

      End FUNCTION CLEBSH2

!====================================================================
      Integer(4) Function Iglq(l,iq)
!====================================================================
!     total stat.weight of lq-subshell (Newton's binom)
!--------------------------------------------------------------------
      kl=4*l+2
      if(iq.gt.kl) Stop 'Iglq:  iq > max.'
      S=1.0
      Do i=iq+1,kl
       S=S*i/(i-iq)
      End do
      Iglq=S+0.5

      End Function Iglq

