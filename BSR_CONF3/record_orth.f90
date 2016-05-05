!======================================================================
      Subroutine Record_orth
!======================================================================
! ... record orthogonal conditions for sct. orbitals:
!----------------------------------------------------------------------
      USE bsr_conf
      USE target; USE channel; USE conf_LS;; USE orb_LS; USE phys_orb_LS

      Implicit none
      Integer :: ich,io, i,j, ii,jj, it,irec

      write(AF,'(a,i3.3)') 'cfg.',ilsp
      Open (nuc,file=AF,position='APPEND')
      write(nuc,'(/a/)') 'Derived orth. conditions:'

! ... record the orth. conditions with indication of main compensation
! ... configuration (see pert_comp.nnn for total coefficients)      

      Do ich=1,nch;  ii = ipch(ich); it = iptar(ich)
       irec = 0 
       Do io = 1,nphys_sub; jj = jp_sub(io)
        if(LEF(ii).ne.LEF(jj)) Cycle
        i = max(ii,jj); j = min(ii,jj)
        if(IORT(i,j).eq.0) then
         write(nuc,'(a1,a4,a1,a4,a3,3x,a12,5x,a)') &
              '<',ELF(ii),'|',ELF(jj),'>=0',trim(AFT(it)),trim(BFT(it))
         irec = 1
        end if
       End do
       if(irec.gt.0) write(nuc,*)
      End do

      End Subroutine Record_orth
