!======================================================================
      Subroutine R_orth(nuc)
!======================================================================
!
!     read from file nuc and store in array IORT (module CONFIGS)
!     the imposed orthogonal conditions 
!
!----------------------------------------------------------------------

      Use configs

      Character(13) :: Aort

      rewind(nuc)
    1 read(nuc,'(a)',end=2) Aort
      if(Aort(1:1).ne.'<') go to 1 
      Call EL4_nlk(Aort(2: 5),n1,l1,k1)
      i1 = Ifind_orb(n1,l1,k1)
      if(i1.eq.0.and.k1.ne.0) go to 1
      Call EL4_nlk(Aort(7:10),n2,l2,k2)
      i2 = Ifind_orb(n2,l2,k2)
      if(i2.eq.0.and.k2.ne.0) go to 1
      if(l1.ne.l2) go to 1 
      read(Aort(13:13),'(i1)') ii

      if(k1.gt.0.and.k2.gt.0) then
       i=max(i1,i2); j=min(i1,i2); IORT(i,j)=ii
      elseif(k1.gt.0.and.k2.eq.0) then
       Do i2=1,nwf
        if(l2.ne.LEF(i2).or.n2.ne.NEF(i2)) Cycle
        i=max(i1,i2); j=min(i1,i2); IORT(i,j)=ii
       End do       
      elseif(k1.eq.0.and.k2.gt.0) then
       Do i1=1,nwf
        if(l1.ne.LEF(i1).or.n1.ne.NEF(i1)) Cycle
        i=max(i1,i2); j=min(i1,i2); IORT(i,j)=ii
       End do       
      elseif(k1.eq.0.and.k2.eq.0) then
       Do i1=1,nwf
        if(l1.ne.LEF(i1).or.n1.ne.NEF(i1)) Cycle
       Do i2=1,nwf
        if(l2.ne.LEF(i2).or.n2.ne.NEF(i2)) Cycle
        i=max(i1,i2); j=min(i1,i2); IORT(i,j)=ii
       End do; End do       
      end if

      go to 1
    2 Continue

      End Subroutine R_orth 

