!----------------------------------------------------------------------
      Integer FUNCTION JDET_SIMP(kz,kn,N1,N2)
!----------------------------------------------------------------------
!     simplify the determinant kn,N1,N2
!     accoding the orthogonality conditions for radial w.f.'s
!     kz - number of needed permutations (added !!!)
!     JDET_SIMP = 0,1,2 with overlap determinant = 0,1 or some value  
!----------------------------------------------------------------------

      Use orb_LS, ONLY: nwf,IORT

      Implicit none

      Integer :: kz,kn
      Integer, Dimension(kn) :: N1(kn),N2(kn)
      Integer :: i,j,ii,i1,i2, k,kk,k1,k2, m1,m2 

      if(kn.le.0) Stop ' JDET_SIMP: kn <= 0'
      i = minval(N1); j=minval(N2)
      if(min(i,j).lt.1) Stop 'JDET_SIMPL: orbital index < 0' 
      i = maxval(N1); j=maxval(N2)
      if(max(i,j).gt.nwf) Stop 'JDET_SIMPL: orbital index > nwf' 

      JDET_SIMP=0

!----------------------------------------------------------------------
!                       Check for a row with only one non-zero element:
    1  Do i1=1,kn                
       k=0                      
       Do i2=1,kn
        m1=max(N1(i1),N2(i2)); m2=min(N1(i1),N2(i2)); ii=IORT(m1,m2)
        if(ii.ne.0) then
         k=k+1; kk=ii;  if(k.gt.1.or.kk.ne.1) Exit
         k1=i1; k2=i2
        end if
       End do
       if(k.eq.0) Return; if(k.eq.1.and.kk.eq.1) go to 2
      End do

!----------------------------------------------------------------------
!                   Check for a coulomb with only one non-zero element:

      Do i2=1,kn                
       k=0                    
       Do i1=1,kn
         m1=max(N1(i1),N2(i2)); m2=min(N1(i1),N2(i2)); ii=IORT(m1,m2)
        if(ii.ne.0) then
         k=k+1; kk=ii; if(k.gt.1.or.kk.ne.1) Exit
         k1=i1; k2=i2
        end if
       End do
       if(k.eq.0) Return; if(k.eq.1.and.kk.eq.1) go to 2
      End do

      go to 3
!-----------------------------------------------------------------------
!                                                 the case of <k1|k2>=1:
    2 kn=kn-1                       
      if(kn.eq.0) then
       JDET_SIMP=1; Return
      end if
      kz=kz+k1+k2
      Do i=k1,kn; N1(i)=N1(i+1); End do
      Do i=k2,kn; N2(i)=N2(i+1); End do
      go to 1
!-----------------------------------------------------------------------
!                                                  ordering of elements:
    3 Continue                     

      Do i1=1,kn-1
       Do i2=i1+1,kn
        if(N1(i1).gt.N1(i2)) then
         kk=N1(i1);  N1(i1)=N1(i2); N1(i2)=kk; kz=kz+1
        end if
        if(N2(i1).gt.N2(i2)) then
         kk=N2(i1);  N2(i1)=N2(i2); N2(i2)=kk; kz=kz+1
        end if
       End do
      End do

      JDET_SIMP=2

      End Function JDET_SIMP


