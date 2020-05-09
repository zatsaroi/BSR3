!=======================================================================
      Subroutine DD_OUT
!=======================================================================
!     output full d-matrix
!-----------------------------------------------------------------------
      Use bsr_dmat 
      Use dmatrix

      Implicit none
      Real(8), allocatable :: rsol1(:,:), rsol2(:,:)
      Real(8), allocatable :: BBL(:,:), AAL(:,:), BBV(:,:), AAV(:,:)

      Real(8), allocatable ::  eval(:)
      Real(8) :: S,SL,SV,E1,E2
      Integer :: i,n,m,nhm,khm,kch,kcp,noterm, n1,n2, m1,m2, ik,jk
      Integer :: JLT1,JST1,JLT2,JST2,JPT1,JPT2

      write(AF,'(a,i3.3,a,i3.3)') 'd.',ilsp1,'_',ilsp2      

      Call Read_aarg('AF_dd',AF)

      open(nudd,file=AF,form='UNFORMATTED')

!-----------------------------------------------------------------------
! ... read initial states: inner region R-matrix solutions

      i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALS1
      Call Check_file(AF)
      Open(nuo,file=AF,form='UNFORMATTED',STATUS='OLD')
      read(nuo) nhm,khm,kch,kcp
      if(nhm.ne.kdm1) Stop 'DD_out, rsol: nhm <> kdm1 '
      if(kch.ne.kch1) Stop 'DD_out, rsol: kch <> kch1 '
      if(kcp.ne.kcp1) Stop 'DD_out, rsol: kcp <> kcp1 '
      if(allocated(eval)) Deallocate(eval);  Allocate(eval(khm)) 
      read(nuo) eval
      nstate1 = khm
      Allocate(rsol1(nstate1,kdm1)) 
      Do i = 1,nstate1
       read(nuo) rsol1(i,:)
      End do
      Close(nuo)

!-----------------------------------------------------------------------
! ... read final states: inner region R-matrix solutions

      i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALS2
      Call Check_file(AF)
      Open(nuo,file=AF,form='UNFORMATTED',STATUS='OLD')
      read(nuo) nhm,khm,kch,kcp
      if(nhm.ne.kdm2) Stop 'DD_out, rsol: nhm <> kdm2 '
      if(kch.ne.kch2) Stop 'DD_out, rsol: kch <> kch2 '
      if(kcp.ne.kcp2) Stop 'DD_out, rsol: kcp <> kcp2 '
      if(allocated(eval)) Deallocate(eval);  Allocate(eval(khm)) 
      read(nuo) eval
      nstate2 = khm
      Allocate(rsol2(kdm2,nstate2)) 
      Do i = 1,nstate2
       read(nuo) rsol2(:,i)
      End do
      Close(nuo)

!-----------------------------------------------------------------------
! ... output:

      JPT1 = (IPT1-1)/2
      JPT2 = (IPT2-1)/2
      if(jmode.eq.0) then
       JLT1 = (ILT1-1)/2; JST1 = IST1 
       JLT2 = (ILT2-1)/2; JST2 = IST2 
      else
       JLT1 = JOT1-1; JLT2 = JOT2-1; JST1=0; JST2=0
      end if

      noterm = kns - 1
      if(noterm.gt.nstate1) noterm=nstate1
      if(noterm.gt.nstate2) noterm=nstate2

      write(nudd) noterm,nstate2,kch2,JLT2,JPT2,JST2,nstate1,kch1,JLT1                          

!-----------------------------------------------------------------------
! ... D(kdm1,kdm2) * rsol2(kdm2,nstate2)

      Allocate (BBL(kdm1,nstate2),BBV(kdm1,nstate2))

      BBL = MATMUL (DL,rsol2)
      BBV = MATMUL (DV,rsol2)
      Deallocate (rsol2,DL,DV)

! ... rsol1(nstate1,kdm1) * BB(kdm1,nstate2)  

      Allocate (AAL(nstate1,nstate2),AAV(nstate1,nstate2))

      AAL = MATMUL (rsol1,BBL)
      AAV = MATMUL (rsol1,BBV)

      Deallocate (rsol1,BBL,BBV)
      Allocate (DL(nstate2,nstate1),DV(nstate2,nstate1))

      DL = TRANSPOSE(AAL)
      DV = TRANSPOSE(AAV)

!      DL = INVERT(AAL)
!      DV = INVERT(AAV)

      Deallocate(AAL,AAV)

!-----------------------------------------------------------------------
! ... output:

! ... d-matrix:

      n2 = 0
      Do jk = 1, 1 + (nstate1-1)/noterm
       n1=n2+1; n2 = min(n2+noterm,nstate1)
       m2 = 0
       Do ik = 1, 1 + (nstate2-1)/noterm
        m1=m2+1; m2 = min(m2+noterm,nstate2)
        write(nudd) DL(m1:m2,n1:n2), DV(m1:m2,n1:n2)
       End do
      End do

! ... Buttle corrections:

      write(nudd) 0.d0
      write(nudd) 0.d0
      write(nudd) 0.d0

! ... Clebsch-Gordan:

      m = min(JLT1,JLT2)
      write(nudd) m,(0.d0,n=1,m)

! ... outer region:

      Allocate(AAL(kch2,kch1), BBL(kch2,kch1), BBV(kch2,kch1))
      AAL = TRANSPOSE(AC)
      BBL = TRANSPOSE(BLC)
      BBV = TRANSPOSE(BVC)

      write(nudd) ((AAL(m,n),n=1,kch1),m=1,kch2), &
                  ((BBL(m,n),n=1,kch1),m=1,kch2), &
                  ((BBV(m,n),n=1,kch1),m=1,kch2)

      Close(nudd)

      End Subroutine DD_out

