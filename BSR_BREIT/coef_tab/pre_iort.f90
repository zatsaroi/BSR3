!---------------------------------------------------------------------
      Subroutine Pre_iort(nu)
!---------------------------------------------------------------------
!
!     prepares the orthogonality conditions (IORT) for all orbitals
!     according to the orbital set numbers KEF, and additional
!     conditions from c-file (unit nu)
!
!     It is assumed that the orbitals of one set are orthogonal,
!     and all orbitals are orthonomal to the ones with KEF=0
!
!     JORT<0 - full orthogonality
!     JORT=0 - full non-orthogonality
!     JORT>0 - partial orthogonality
!
!---------------------------------------------------------------------

      USE configs

      IMPLICIT NONE
      
      Integer(4), Intent(in) :: nu
      Integer(4) :: i1,i2
      
      Do i1=1,nwf
      Do i2=1,i1

       IORT(i1,i2)=0
       if(LEF(i1).ne.LEF(i2)) Cycle

       if(JORT.lt.0) then             ! full orthogonality

        IORT(i1,i2)=0
        if(i1.eq.i2) IORT(i1,i2)=1

       elseif(JORT.eq.0) then         ! full non-orthogonality

        IORT(i1,i2)=2
        if(i1.eq.i2.and.NEF(i1).le.nmax) IORT(i1,i2)=1

       elseif(JORT.gt.0) then         ! partial orthogonality

        if(KEF(i1).eq.KEF(i2)) then
         IORT(i1,i2)=0
         if(NEF(i1).eq.NEF(i2).and.NEF(i1).le.nmax) IORT(i1,i2)=1         
         if(NEF(i1).eq.NEF(i2).and.NEF(i1).gt.nmax) IORT(i1,i2)=2         
        elseif(KEF(i1)*KEF(i2).eq.0) then
         IORT(i1,i2)=0
        else
         IORT(i1,i2)=2
        end if

       end if

      End do
      End do

! ... additional orthogonality conditions 

      if(nu.gt.0) Call R_orth(nu)
      
      End Subroutine Pre_iort
