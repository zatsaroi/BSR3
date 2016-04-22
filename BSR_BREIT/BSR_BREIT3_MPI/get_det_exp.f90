!======================================================================
      Subroutine get_det_exp(ic,jc)
!======================================================================

      USE MPI

      USE bsr_breit, only: ierr,noper,joper,JT_oper,CT_oper
      USE conf_LS,   only: ne
      USE coef_list, only: ntrm
      USE spin_orbitals, only: NNsym1,NNsym2, Lsym1,Lsym2
      USE term_exp
          
      Implicit none

      Integer :: status(MPI_STATUS_SIZE), ic,jc

      Call MPI_RECV(ic,1,MPI_INTEGER, 0,0, MPI_COMM_WORLD, status, ierr)

      if(ic.le.0) Return

      Call MPI_RECV(kt1 ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(kdt1,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(ip_kt1)) Deallocate(ip_kt1); Allocate(ip_kt1(kt1))
      Call MPI_RECV(ip_kt1,kt1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(IM_det1)) Deallocate(IM_det1); Allocate(IM_det1(ne,kdt1))
      Call MPI_RECV(IM_det1,ne*kdt1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(IS_det1)) Deallocate(IS_det1); Allocate(IS_det1(ne,kdt1))
      Call MPI_RECV(IS_det1,ne*kdt1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(C_det1)) Deallocate(C_det1); Allocate(C_det1(kt1,kdt1))
      Call MPI_RECV(C_det1,kt1*kdt1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(jc  ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(kt2 ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(kdt2,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(ip_kt2)) Deallocate(ip_kt2); Allocate(ip_kt2(kt2))
      Call MPI_RECV(ip_kt2,kt2,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(IM_det2)) Deallocate(IM_det2); Allocate(IM_det2(ne,kdt2))
      Call MPI_RECV(IM_det2,ne*kdt2,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(IS_det2)) Deallocate(IS_det2); Allocate(IS_det2(ne,kdt2))
      Call MPI_RECV(IS_det2,ne*kdt2,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      if(allocated(C_det2)) Deallocate(C_det2); Allocate(C_det2(kt2,kdt2))
      Call MPI_RECV(C_det2,kt2*kdt2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(MLT,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(MST,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(ILT1,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(IST1,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(ILT2,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(IST2,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(nnsym1,ne,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(nnsym2,ne,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(lsym1,ne,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(lsym2,ne,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)

      Call MPI_RECV(ntrm,1,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      Call MPI_RECV(joper,noper,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      if(allocated(JT_oper)) Deallocate(JT_oper); Allocate(JT_oper(ntrm,noper))
      Call MPI_RECV(JT_oper,ntrm*noper,MPI_INTEGER,0,0,MPI_COMM_WORLD, status,ierr)
      if(allocated(CT_oper)) Deallocate(CT_oper); Allocate(CT_oper(ntrm,noper))

      End Subroutine get_det_exp
