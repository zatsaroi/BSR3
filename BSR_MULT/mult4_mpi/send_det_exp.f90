!======================================================================
      Subroutine Send_det_exp(id,ic,jc)
!======================================================================
      Use MPI

      Use mult_par,      only: ierr,JT_oper
      Use conf_LS,       only: ne
      Use coef_list,     only: ntrm
      Use spin_orbitals, only: NNsym1,NNsym2, Lsym1,Lsym2
      Use term_exp

      Implicit none
      
      Integer :: id,ic,jc 

      Call MPI_SEND(ic  ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      if(ic.le.0) Return

      Call MPI_SEND(kt1 ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kdt1,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ip_kt1,kt1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IM_det1,ne*kdt1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IS_det1,ne*kdt1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(jc  ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(kt2 ,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kdt2,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ip_kt2,kt2,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IM_det2,ne*kdt2,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IS_det2,ne*kdt2,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(MLT1,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(MST1,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(MLT2,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(MST2,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ILT1,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IST1,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ILT2,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IST2,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(nnsym1,ne,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(nnsym2,ne,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(lsym1,ne,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(lsym2,ne,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(ntrm,1,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(JT_oper,ntrm,MPI_INTEGER,id,0,MPI_COMM_WORLD,ierr)

      End Subroutine send_det_exp
