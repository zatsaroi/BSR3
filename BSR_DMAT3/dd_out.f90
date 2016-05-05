!=======================================================================
      Subroutine DD_OUT
!=======================================================================
!     output full d-matrix
!-----------------------------------------------------------------------
      Use bsr_dmat 
      Use dmatrix

      Implicit none
      Integer :: i

      i = INDEX(AF_d,'.'); AF = AF_d(1:i)//ALS1//'_'//ALS2
      Open(nud,file=AF,form='UNFORMATTED')

      write(nud) kdm1,kdm2
      Do i=1,kdm2; write(nud) DL(1:kdm1,i); End do
      Do i=1,kdm2; write(nud) DV(1:kdm1,i); End do

      Close(nud)

      End Subroutine DD_out

