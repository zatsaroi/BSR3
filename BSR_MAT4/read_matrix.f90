!======================================================================
      Subroutine Read_matrix
!======================================================================
!     Read interaction/overlap matrix from the file 'nui';
!     matrix is added to the memory data (do we need that ???)
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer ::  ic,jc, ii

      Do 
       read(nui) ic,jc
       if(ic.le.0) Exit

      if(ic.gt.nch.and.jc.gt.nch) then  !  pert-pert

       ic=ic-nch; jc=jc-nch; ii = ibb(ic,jc)
       read(nui) hbb(ii)
        
      elseif(ic.gt.nch) then            !  ch-pert

       ic=ic-nch; ii = icb(jc,ic)
       read(nui) hcb(:,ii)

      else                              !  ch-ch

       ii = icc(ic,jc)
       read(nui) hcc(:,:,ii)

      end if 

      End do

      End Subroutine read_matrix

!======================================================================
      Subroutine Skip_matrix
!======================================================================
!     skip overlap matrix from the file 'nui'
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Real(8) :: S, w(ns), ww(ns*ns)
      Integer ::  ic,jc

      Do 
       read(nui) ic,jc; if(ic.le.0) Exit
       if(ic.gt.nch.and.jc.gt.nch) then  !  pert-pert
        read(nui) S
       elseif(ic.gt.nch) then            !  ch-pert
        read(nui) w
       else                              !  ch-ch
        read(nui) ww
       end if 
      End do

      End Subroutine skip_matrix

