!======================================================================
      Subroutine R_bsrmat
!======================================================================
!     read the data from bsr_mat.nnn file for the given partial wave
!----------------------------------------------------------------------
      Use bsr_pol
      Use spline_param, only: ns
      Use channel,      only: nch,ncp
      
      Implicit none
      Integer :: i, i1,i2,i3

! ... hamiltonian matrix file:

      i=LEN_TRIM(AF_int); AF_int(i-2:i)=ALSP
      Call Check_file(AF_int)
      Open(nui,file=AF_int,form='UNFORMATTED')

! ... check the dimensions:

      read(nui) i1,i2,i3
      if(i1.ne.ns ) Stop ' BSR_POL: different ns  in BSR_MAT file'
      if(i2.ne.nch) Stop ' BSR_POL: different kch in BSR_MAT file'
      if(i3.ne.ncp) Stop ' BSR_POL: different kcp in BSR_MAT file'
      
      nhm = nch*ns + ncp
      mhm = nhm + nort + nortb

      write(pri,'( a,i5,a)') 'nhm  = ',nhm,'  - interaction matrix'
      write(pri,'( a,i5,a)') 'mhm  = ',mhm,'  - full size of matrix'

! ... overlap matrix:

      if(allocated(c)) Deallocate(c); Allocate(c(mhm,mhm))
      Call Read_bsr_matrix(nui,mhm,ns,nch,c)

! ... interaction matrix:

      if(allocated(a)) Deallocate(a); Allocate(a(mhm,mhm))
      Call Read_bsr_matrix(nui,mhm,ns,nch,a)

      End Subroutine R_bsrmat 


!======================================================================
      Subroutine Read_bsr_matrix(nui,mhm,ns,kch,a)
!======================================================================
!     read the inreraction (overlap) matrix in BSR-format
!----------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: nui,mhm,ns,kch
      Integer :: i,i1,i2, j,j1,j2, ic,jc, k
      Real(8) :: a(mhm,mhm)
       
      a = 0.d0;  k = ns*kch-kch

! ... diagonal blocks:

      Do i=1,kch; i1=(i-1)*ns+1; i2=i*ns
       read(nui) a(i1:i2,i1:i2)
      End do

! ... other elements if any:

      Do 
       read(nui) ic,jc;  if(ic.le.0) Exit
       if(ic.gt.kch.and.jc.gt.kch) then            !  pert-pert
        read(nui) a(ic+k,jc+k)
       elseif(ic.gt.kch) then                      !  ch-pert
        i=ic+k; j1=(jc-1)*ns+1; j2=jc*ns
        read(nui) a(i,j1:j2)
       else                                        !  ch-ch
        i1=ns*(ic-1)+1; i2=ns*ic
        j1=ns*(jc-1)+1; j2=ns*jc
        read(nui) a(i1:i2,j1:j2)
       end if
      End do

! ... get upper-half from symmetry:

      Do i=1,mhm; Do j=1,i; a(j,i)=a(i,j);  End do; End do

      End Subroutine Read_bsr_matrix
