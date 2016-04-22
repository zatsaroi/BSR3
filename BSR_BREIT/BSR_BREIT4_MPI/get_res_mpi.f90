!======================================================================
      Subroutine GET_res(id,is,js)
!======================================================================
! ... get results from processor 'id'
!----------------------------------------------------------------------
      USE MPI

      USE bsr_breit, only: ierr,noper,koper,JD_oper, pri
      USE term_exp,  only: jt1,jt2, JP_kt1,JP_kt2
      USE coef_list, only: ktrm,ncoef,idfc,intc,coef
      USE ndet_list
      USE ndef_list

      Implicit none
      
      Integer :: id, is,js, status(MPI_STATUS_SIZE)

! ... coeffs:

      Call MPI_RECV(ncoef,1,MPI_INTEGER, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      id = status(MPI_SOURCE)

      Call MPI_RECV(is,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(js,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(jt1,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(jt2,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(ktrm,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(koper,noper,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if(allocated(JD_oper)) Deallocate(JD_oper)
      Allocate(JD_oper(ktrm,noper))
      Call MPI_RECV(JD_oper,ktrm*noper,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

! ... term indexes:

      if(allocated(JP_kt1)) Deallocate(JP_kt1); Allocate(JP_kt1(jt1))
      Call MPI_RECV(JP_kt1,jt1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if(allocated(JP_kt2)) Deallocate(JP_kt2); Allocate(JP_kt2(jt2))
      Call MPI_RECV(JP_kt2,jt2,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(ncoef.eq.0) Return

! ... coefficients:

      if(allocated(coef)) Deallocate(coef); Allocate(coef(ktrm,ncoef))
      if(allocated(idfc)) Deallocate(idfc); Allocate(idfc(ncoef))
      if(allocated(intc)) Deallocate(intc); Allocate(intc(ncoef))

      Call MPI_RECV(idfc,ncoef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(intc,ncoef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(coef,ktrm*ncoef,MPI_DOUBLE_PRECISION, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

! ... determinants:

      Call MPI_RECV(ndet,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(ndet.eq.0) Return

      Call MPI_RECV(ldet,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(allocated(KPD)) Deallocate(KPD); Allocate(KPD(ndet))
      if(allocated(IPD)) Deallocate(IPD); Allocate(IPD(ndet))
      if(allocated(NPD)) Deallocate(NPD); Allocate(NPD(ldet))
      mdet = ndet; kdet=ldet

      Call MPI_RECV(KPD,ndet,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(IPD,ndet,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(NPD,ldet,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

! ... det. factors:

      Call MPI_RECV(ndef,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      if(ndef.eq.0) Return
      Call MPI_RECV(ldef,1,MPI_INTEGER, id, &
                     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if(allocated(KPf)) Deallocate(KPf); Allocate(KPf(ndef))
      if(allocated(IPf)) Deallocate(IPf); Allocate(IPf(ndef))
      if(allocated(NPf)) Deallocate(NPf); Allocate(NPf(ldef))
      mdef = ndef; kdet=ldef

      Call MPI_RECV(KPf,ndef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(IPf,ndef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(NPf,ldef,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      End Subroutine get_res
