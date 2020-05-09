!======================================================================
      Subroutine rw_mask
!======================================================================
!     Read interaction matrix from the disk
!----------------------------------------------------------------------
      Use bsr_mat
      Use L_core

      Implicit none
      Real(8) :: S, a(ns,ns), aaa(nch,nch,0:mk)

      Integer :: ich,jch, i1,i2,i3,ii

      read(nud) i1,i2,i3
      if(i1.ne.ns   ) Stop 'different ns  in BSR_MAT_mask file'
      if(i2.ne.nch  ) Stop 'different nch in BSR_MAT_mask file'
      if(i3.ne.npert) Stop 'different npert in BSR_MAT_mask file'
      if(i3.ne.0    ) Stop 'npert /= 0 in BSR_MAT_mask file'
      write(nui) ns,nch,npert

! ... overlap matrix:

      Do 
       read(nud) ich,jch,ii
       write(nui) ich,jch,ii
       if(ich.le.0) Exit
       read(nud) a
       write(nui) a
      End do

! ... L-integrals:

      Call Gen_Lval
      hl_full = -0.5d0 * hl_full 

! ... Hamiltonian matrix:

      Do 
       read(nud) ich,jch,ii
       write(nui) ich,jch,ii
       if(ich.le.0) Exit
       read(nud) a
       if(ich.eq.jch)  a(:,:) = a(:,:) + hl_full(:,:,lch(ich)) 
       write(nui) a
      End do

! ... asymptotic coefficients:

      read(nud) i1
      read(nud) aaa
      read(nud) S
      if(i1.ne.mk) Stop 'different mk in BSR_MAT_mask file'
      write(nui) mk;  write(nui) aaa;  write(nui) t(ns+1),ns

      ich = (nch+1)*nch/2
      if(allocated(htarg)) Deallocate(htarg); Allocate(htarg(ich))
      if(allocated(otarg)) Deallocate(otarg); Allocate(otarg(ich))

      read(nud) htarg;     write(nui) htarg
      read(nud) otarg;     write(nui) otarg
      read(nud) Etarg;     write(nui) Etarg
      read(nud) EC;        write(nui) EC

      close(nud);          close(nui)

      End Subroutine rw_mask

