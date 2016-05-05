!======================================================================
      Integer Function Idet_fact (j1,j2,j3,j4)
!======================================================================
!     determines the overlap factor and its position in NDEF list
!     for matrix element between two determinant wave functions
!     located in module 'spin_orbitals' 
!
!     j1,j2,j3,j4 - interacting electrons
!
!     Calls:  Iadd_det, Iadd_def, ISORT
!----------------------------------------------------------------------
      Use conf_LS,       only: ne
      Use det_list,      only: ibd
      Use def_list,      only: ibf
      Use spin_orbitals, only: nsym, IPsym1,NNsym1,Isym1, &
                                     IPsym2,NNsym2,Isym2

      Implicit none
      Integer, Intent(in) :: j1,j2,j3,j4
      Integer :: i,i1,i2, k,k1,k2, is,id,kd
      Integer :: np(2*ne),n1(2*ne),n2(2*ne),n3(2*ne),n4(2*ne)     
      Integer, External :: Iadd_ndet, Iadd_ndef, ISORT

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

       NP(1:k1) = N1(1:k1)*ibd + N2(1:k1);  id = Iadd_ndet(k1,NP)

       Do k=1,kd
        if(N3(k).ne.id) Cycle;  N4(k)=N4(k)+1; id=0; Exit
       End do
       if(id.eq.0) Cycle

       kd=kd+1; N3(kd)=id; N4(kd)=1

      End do

      Idet_fact=0; if(kd.eq.0) Return

      NP(1:kd)=N3(1:kd)*ibf + N4(1:kd)

      k=ISORT(kd,NP)

      Idet_fact = Iadd_ndef(kd,NP)

      End Function Idet_fact

