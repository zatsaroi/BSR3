!=======================================================================
      Subroutine Mult_2conf(no1,nn1,ln1,iq1,LS1,no2,nn2,ln2,iq2,LS2, &
                            kpol, ck,ik,jk)
!=======================================================================
!     compute the angular coefficients for multipole operator "kpol"
!     between 2 atomic states in case of orthogonal orbitals.
!     In this case, we are expecting only 1 radial integral:
!          ck  * INT[ P_ik * r^k * P_jk, r=0,inf]
!-----------------------------------------------------------------------
      Implicit none 

! ... input-output:

      Integer, intent(in)  :: no1,nn1(no1),ln1(no1),iq1(no1),LS1(5,no1)
      Integer, intent(in)  :: no2,nn2(no2),ln2(no2),iq2(no2),LS2(5,no2)
      Integer, intent(in)  :: kpol
      Integer, intent(out) :: ik, jk
      Real(8), intent(out) :: ck

! ... determinant expansion:

      Integer :: kdt1, kdt2
      Real(8), allocatable :: Cdet1(:), Cdet2(:)
      Integer, allocatable :: MLdet1(:,:), MSdet1(:,:), &
                              MLdet2(:,:), MSdet2(:,:) 
      Real(8) :: CC_det, eps_C = 1.d-6

! ... local variables:

      Integer :: ne, kd1, kd2, i, k, mdt1, mdt2, i1,i2,i3, j1,j2, qpol, &
                 LT1, LT2, MT1, MT2
      Integer, allocatable :: nl1(:), nl2(:), ip1(:), ip2(:)
      Integer, external :: Iglq
      Real(8) :: CN
      Real(8) , external :: Z_3j   


! ... check the selection rules for electric multipole transition:

      ck = 0.d0; ik=0; jk=0; if(kpol.le.0) Return
      i1=(-1) ** SUM(ln1(1:no1)*iq1(1:no1)) 
      i2=(-1) ** SUM(ln2(1:no2)*iq2(1:no2)) 
      if(i1.eq.i2.and.mod(kpol,2).ne.0) Return
      if(i1.ne.i2.and.mod(kpol,2).ne.1) Return
      if(LS1(5,no1).ne.LS2(5,no2)) Return
      i1 = (LS1(4,no1)-1)/2;  i2 = (LS2(4,no2)-1)/2; i3=kpol
      if(i1.gt.i2+i3.or.i2.gt.i1+i3.or.i3.gt.i1+i2) Return
      if(i1.lt.iabs(i2-i3).or.i2.lt.iabs(i1-i3).or. &
         i3.lt.iabs(i1-i2)) Return
      LT1 = LS1(4,no1);  MT1 = LT1
      LT2 = LS2(4,no2);  MT2 = LT2

! ... initialize arrays:

      ne = SUM(iq1(1:no1))
      if(ne.ne.SUM(iq2(1:no2))) Stop 'Coef_ee_2conf: ne1 <> ne2' 
      if(allocated(nl1)) Deallocate(nl1,ip1); Allocate(nl1(ne),ip1(ne))
      if(allocated(nl2)) Deallocate(nl2,ip2); Allocate(nl2(ne),ip2(ne))
      
      k=1; Do i=1,no1; nl1(k:k+iq1(i)-1)=nn1(i)*1000+ln1(i); k=k+iq1(i); End do
      k=1; Do i=1,no1; ip1(k:k+iq1(i)-1)=i; k=k+iq1(i); End do
      k=1; Do i=1,no2; nl2(k:k+iq2(i)-1)=nn2(i)*1000+ln2(i); k=k+iq2(i); End do
      k=1; Do i=1,no2; ip2(k:k+iq2(i)-1)=i; k=k+iq2(i); End do

! ... determinant expansion 1:

      mdt1=1;  Do i=1,no1;   mdt1=mdt1*Iglq(ln1(i),iq1(i)); End do
      if(allocated(Cdet1) ) Deallocate(Cdet1 );  Allocate(Cdet1(mdt1)    )
      if(allocated(MLdet1)) Deallocate(MLdet1);  Allocate(MLdet1(ne,mdt1))
      if(allocated(MSdet1)) Deallocate(MSdet1);  Allocate(MSdet1(ne,mdt1))

      Call Det_expn_1conf(no1,ln1,iq1,LS1,ne,mdt1,kdt1,Cdet1,MLdet1,MSdet1)  

! ... determinant expansion 2:

      mdt2=1;  Do i=1,no2;   mdt2=mdt2*Iglq(ln2(i),iq2(i)); End do
      if(allocated(Cdet2) ) Deallocate(Cdet2 );  Allocate(Cdet2(mdt2)    )
      if(allocated(MLdet2)) Deallocate(MLdet2);  Allocate(MLdet2(ne,mdt2))
      if(allocated(MSdet2)) Deallocate(MSdet2);  Allocate(MSdet2(ne,mdt2))

      Call Det_expn_1conf(no2,ln2,iq2,LS2,ne,mdt2,kdt2,Cdet2,MLdet2,MSdet2)  

! ... calculations:     

      qpol = (LS1(4,no1)-LS2(4,no2))/2

      CN = Z_3j(LT1,-MT1+2,2*kpol+1,MT1-MT2+1,LT2,MT2) * (-1)**((LT1-MT1)/2)

      if(CN.eq.0.d0) Return

      Do kd1 = 1,kdt1
      Do kd2 = 1,kdt2
       CC_det = Cdet1(kd1) * Cdet2(kd2)
       Call me_det

write(*,*) kd1,kd2, CC_det

      End do;  End do

write(*,*) 'CN,CK',CN, CK

      CK = CK / CN


CONTAINS

!======================================================================
      Subroutine Det_expn_1conf(no,ln,iq,LS,ne,mdt,kdt,Cdet,MLdet,MSdet)    
!======================================================================
!     determined all possible determinants and their coefficients
!     for configuration given in calling routine
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: no,ln(no),iq(no),LS(5,no), ne,mdt
      Integer, intent(out) :: MLdet(ne,mdt), MSdet(ne,mdt)
      Real(8), intent(out) :: Cdet(mdt)
     
      Integer :: idet(ne)
      Integer :: ML(no),MS(no),MLp(no),MSp(no)  
      Integer :: kdt,i,j,k,m,ii, MLT, MST
      Real(8) :: C, eps_det = 1.d-8
      Integer :: in(no),md(no),nd(no),it(no)
      Integer, external :: Iglq, Iterm_LS, ML_id, MS_id
      Real(8), external :: Clebsh, DETC_sh

! ... auxiliary arrays:
! ... in(i)  -  pointer to given shell
! ... md(i)  -  number of determinant for given shell
! ... it(i)  -  pointer to shell term
! ... mdt    -  max. possible number of determinants

      k=1 
      Do i=1,no
       in(i)=k; k=k+iq(i); md(i)=Iglq(ln(i),iq(i)) 
	      it(i)=Iterm_LS(ln(i),iq(i),0,LS(1,i),LS(2,i),LS(3,i)) 
      End do

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
      if(abs(C).lt.eps_det) go to 2

      kdt=kdt+1;  if(kdt.gt.mdt) Stop 'kdt > mdt'
      Cdet(kdt) = C
      Do j = 1,ne  
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
!     call the subroutine to calculate m.e. between possible 
!     combinations of nj-orbitals ()
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,i1,i2, j,j1,j2, m, idif, jdif, is,js, &
                 k1,k2,k3,k4
      Real(8), external :: D_matr

! ... find interaction orbitals:

       idif=0
       Do i1=1,ne
        m=0
        Do i2=1,ne 
         if(nl1(i1).ne.nl2(i2)) Cycle
         if(MLdet1(i1,kd1).ne.MLdet2(i2,kd2)) Cycle
         if(MSdet1(i1,kd1).ne.MSdet2(i2,kd2)) Cycle
         m=1; Exit 
        End do
        if(m.eq.1) Cycle  
        idif=idif+1
        if(idif.gt.1) Return
        i=i1
       End do

       jdif=0
       Do i2=1,ne
        m=0
        Do i1=1,ne 
         if(nl1(i1).ne.nl2(i2)) Cycle
         if(MLdet1(i1,kd1).ne.MLdet2(i2,kd2)) Cycle
         if(MSdet1(i1,kd1).ne.MSdet2(i2,kd2)) Cycle
         m=1; Exit 
        End do
        if(m.eq.1) Cycle  
        jdif=jdif+1
        if(jdif.gt.1) Return
        j=i2
       End do

       if(idif.ne.jdif) Stop 'me_det: problems with determinants'

       if(idif.ne.1) Return

       if(ik.eq.0) ik = ip1(i)
       if(jk.eq.0) jk = ip2(j)
       if(ip1(i).ne.ik) Stop 'me_det: problems with int. orbitals' 
       if(ip2(j).ne.jk) Stop 'me_det: problems with int. orbitals' 

       CK = CK + CC_det * D_matr(kpol,qpol,ln1(ik),MLdet1(i,kd1), &
                                           ln2(jk),MLdet2(j,kd2)) &
                        * (-1) ** (i+j)
       End Subroutine me_det


       End Subroutine Mult_2conf


!====================================================================
      Real(8) Function D_matr(kpol,qpol,l1,m1,l2,m2)
!====================================================================
!     angular part of electric transition operator 
!     between 'nlms' orbitals:
!
!            <n1,l1,m1,s1| T(kq) | n2,l2,m2,s2>
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: kpol,qpol,l1,m1,l2,m2
   	  Real(8), external :: Z_3jj, ZCLKL

      D_matr =  (-1)**(l1-m1) * Z_3jj(l1,-m1,kpol,qpol,l2,m2) &
                              * ZCLKL(l1,kpol,l2)
   	  End Function D_matr

!======================================================================
      Integer Function Iterm_LS (l,iq,k,IA,IL,IS)
!======================================================================
!     provides information about shell terms:
!
!     k>0  --> IA,IL,IS = k-th term of lq-subshell
!     k=0  --> Iterm_LS = position of the IA,IL,IS term in subshell list
!     k=-1 --> Iterm_LS = number of terms in lq-subshell
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
      
      Integer, intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebsh2 = Clebsh(l1+1,m1+1,l2+1,m2+1,l+1,m+1)

      End FUNCTION CLEBSH2



!======================================================================
      Real(8) Function ZCB(K1,K2,K3)
!======================================================================
!     CB =  3j(k1,0,k2,0,k3,0)**2
!
!     3j =  sqrt[ (2g-2k1)! (2g-2k2)! (2g-2k3)! / (2g+1)! ] *
!                   g! / [ (g-k1)! (g-k2)! (g-k3)! ],
!
!     where 2g = k1+k2+k3
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k1,k2,k3
      Real(8) :: A,B,C
      Integer :: I,J, K,L, M,N, N1,N2,N3, M1,M2,M3

      ZCB=0.0;  M=K1+K2+K3;  N=M/2; if(N+N.NE.M) Return
      N1=N-K1; N2=N-K2; N3=N-K3
      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) Return

      M1=N1+N1; M2=N2+N2; M3=N3+N3
      A = 1.d0/(M + 1)
      DO I = 1,N
       K = +1
       IF(M1.GE.I) K=K+1
       IF(M2.GE.I) K=K+1
       IF(M3.GE.I) K=K+1
       IF(N1.GE.I) K=K-2
       IF(N2.GE.I) K=K-2
       IF(N3.GE.I) K=K-2
       J = I + N
       L = -1 
       IF(M1.GE.J) L=L+1 
       IF(M2.GE.J) L=L+1 
       IF(M3.GE.J) L=L+1 
       B=I;  C=J;  A = A * B**K * C**L
      End do
      ZCB = A

      END FUNCTION ZCB



!======================================================================
      Real(8) Function ZCLKL(K1,K2,K3)
!======================================================================
!     reduced matrix elements of spherical garmonics:
!
!     <k1||C[k2]||k3> = (-1)^k1 * sqtr[ (2*k1+1) * (2*k3+1) ] *
!                     3j(k1,0,k2,0,k3,0)
!
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(In) :: k1,k2,k3
      Integer :: M, N
      Real(8), external :: ZCB

      M=K1+K2+K3
      N=M/2
      IF(N+N.NE.M) Return

      ZCLKL = ZCB(K1,K2,K3)
      ZCLKL = SQRT(ZCLKL*(K1+K1+1)*(K3+K3+1))
      IF(mod(N+K1,2).eq.1) ZCLKL=-ZCLKL

      End Function ZCLKL


!====================================================================
      Integer Function Iglq(l,iq)
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
