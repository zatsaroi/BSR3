!======================================================================
      Integer(4) Function Idet_fact (j1,j2,j3,j4)
!======================================================================
!
!     determines the overlap factor and its position in NDEF list
!     for matrix element between two determinant wave functions
!     located in module 'spin_orbitals' 
!
!     j1,j2,j3,j4 - interacting electrons
!
!     Calls:  Nadd_det, Nadd_def, ISORT
!
!----------------------------------------------------------------------

      USE param_br, ONLY: ibd,ibf
      USE spin_orbitals

      Implicit none
      Integer(4), Intent(in) :: j1,j2,j3,j4
      Integer(4) :: i,i1,i2, k,k1,k2, is,id,kd
      Integer(4), External :: Nadd_det, Nadd_def, ISORT

      kd = 0
      Do is = 1,NSYM

       i1 = 1; if(is.gt.1) i1 = IPsym1(is-1)+1; i2 = IPsym1(is)
       k1 = 0
       Do i = i1,i2
        if(i.eq.j1.or.i.eq.j2) Cycle
        k1 = k1 + 1;  N1(k1) = nnsym1(Isym1(i))
       End do

       i1 = 1; if(is.gt.1) i1 = IPsym2(is-1)+1; i2 = IPsym2(is)
       k2 = 0
       Do i = i1,i2
        if(i.eq.j3.or.i.eq.j4) Cycle
        k2 = k2 + 1;  N2(k2) = nnsym2(Isym2(i))
       End do

       if(k1.ne.k2) Stop 'Idet_fact: k1.ne.k2'
       if(k1.eq.0.or.k2.eq.0) Cycle

       NP(1:k1) = N1(1:k1)*ibd + N2(1:k1);  id = Nadd_det(k1,NP)

       Do k=1,kd
        if(N3(k).ne.id) Cycle;  N4(k)=N4(k)+1; id=0; Exit
       End do
       if(id.eq.0) Cycle

       kd=kd+1; N3(kd)=id; N4(kd)=1

      End do

      Idet_fact=0; if(kd.eq.0) Return

      NP(1:kd)=N3(1:kd)*ibf + N4(1:kd)

      k=ISORT(kd,NP)

      Idet_fact = Nadd_def(kd,NP)

      End Function Idet_fact

