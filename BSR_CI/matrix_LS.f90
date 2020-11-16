!======================================================================
      Subroutine matrix_LS
!======================================================================
!     generation of interaction and overlap matrices accoding to INT_BNK
!----------------------------------------------------------------------
      Use bsr_ci
      Use c_data
      Use L_core

      Implicit real(8) (A-H,O-Z)

!----------------------------------------------------------------------
! ... allocations of interaction and overlap matrixes:

      Allocate(HM(ncfg,ncfg), HS(ncfg,ncfg))
      HM = 0.d0; HS = 0.d0

! ... overlap matrix: 

      Call Read_int_bnk(11)
     
! ... L-ingegrals:

      Call Gen_Lval

      Call Alloc_c_data(1,0,0,mblock,nblock,kblock,eps_c)
      Call Read_int_bnk(6)
      Call Add_matrix(1,0)

      Do j=1,ncdata; i=IPT(j)
       ic = k3(i);  jc = k4(i)
       ii = k1(i);  jj = k2(i)
       HM(ic,jc) = HM(ic,jc) + Cdata(i)*L_int(ii,jj)
      End do

      Call Alloc_Lcore(0,0,0)

! ... Rk-integrals:

      kmax = km       !  ???

      Call Alloc_c_data(1,0,kmax,mblock,nblock,kblock,eps_c)

      Call Read_int_bnk(5)

      Do k=0,kmax
       Call Add_matrix(1,k)
       Call I_data(5,k) 
      End do

! ... o-o inreraction:

      if(ioo.gt.0) then
       
      Call Alloc_c_data(1,0,kmax,mblock,nblock,kblock,eps_c)

      Call Read_int_bnk(3)
      Do k=0,kmax
       Call Add_matrix(1,k)
       Call I_data(3,k) 
      End do

      Call Read_int_bnk(4)
      Do k=0,kmax
       Call Add_matrix(1,k)
       Call I_data(4,k) 
      End do

      end if  ! over ioo

! ... record the LS matrix into scratch file:

      open(nuh,form='UNFORMATTED',status='SCRATCH')
      write(nuh) HM

!----------------------------------------------------------------------
! ... record the overlap matrix:

      open(nuo,form='UNFORMATTED',status='SCRATCH')
      write(nuo) HS

!----------------------------------------------------------------------
! ... print non-trivial TOTAL overlap matrix elements:

      write(iwrite,'(/a,f10.5/)') &
                   'Non-trivial total overlaps with eps_o:',eps_o
      Do j=1,ncfg
       Do i=j+1,ncfg
        if(abs(HS(i,j)).lt.eps_o) Cycle
           write(iwrite,'(2I8,3f10.3,5x,a,5x,a)') &
                i, j, HS(i,j), WC(i),WC(j), trim(LABEL(i)),trim(LABEL(j))
       End do
      End do
      write(iwrite,'(a)') '***'
      write(iwrite,'(/a,f10.5/)') &
                   'Non-trivial normalization with eps_d:',eps_d
      Do i=1,ncfg
       if(abs(HS(i,i)-1.d0).gt.eps_d) &
          write(iwrite,'(F13.5,2x,2I8)') HS(i,i), i, i
      End do

!----------------------------------------------------------------------
! ... find NORT --> flag for the generelized eigenvalue problem

      NORT = -1

      Do j=1,NZERO;  Do i=j,NZERO
        C = HS(i,j); if(i.eq.j) C = C - 1.D0
        if(abs(C).lt.eps_ovl) Cycle
        NORT = 0;  Exit
      End do; End do

      if(NORT.eq.0) write(iwrite,'(/a/)') &
        'NORT = 0 --> generelized eigenvalue problem'
      if(NORT.eq.-1) write(iwrite,'(/a/)') &
        'NORT = 0 --> simple eigenvalue problem'

      End Subroutine matrix_LS

