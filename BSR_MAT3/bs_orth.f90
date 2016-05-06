!=====================================================================
      Subroutine BS_ORTH
!=====================================================================
!     Including the orthogonal conditions for scattering orbitals
!     by updating the interaction matrix according the procedure
!     suggested by M.Bentley J.Phys.B27 (1994) 637. 
!
!     Let i-th channel is supposed to be orthogonal to orbital p.
!     It can be achieved by following modification of Hamiltonian:
!
!                H -->  (1 - Bcc') H (1 - cc'B) 
!
!     where c - full solution vector with all zero elements except
!     the channel i is replaced on B-spline expansion of orbital p
!     primes denote the transportasion. 
!
!     It leads to transformation of separate blocks  H(i,j):
!
!     H(i,j) with j <> i -->  (1 - Bpp') H(i,j)    
!     H(j,i) with j <> i -->   H(i,j) (1 - pp'B)   
!     H(i,i) with j  = i -->  (1 - Bpp') H(i,j) (1 - pp'B)    
!     j also includes the block for interaction with perturbers
!---------------------------------------------------------------------
      Use bsr_mat, only: pri
      Use bsr_matrix
      Use channel
      Use spline_param
      Use spline_orbitals, p => pbs
      Use spline_galerkin, b => bb
      
      Implicit none
      Character(80) :: line
      Integer :: i,j, ii,jj, ij, ich,jch, nort, iline
      Real(8), allocatable :: y(:,:),xx(:,:),yy(:,:),zz(:,:),z(:)

      write(pri,'(/a/)') 'Orthogonal conditions:'
      Allocate(y(ns,ns),xx(ns,ns),yy(ns,ns),zz(ns,ns),z(ns))
      Do ich = 1,nch;  i = ipch(ich)

       zz = 0.d0
       nort = 0
       write(line,'(a,a)') ELC(ich),' -> '; iline=11 
       Do j = 1,nbf
        if(iech(j).gt.0) Cycle
        if(lbs(i).ne.lbs(j)) Cycle        
        if(IBORT(max(i,j),min(i,j)).ne.0) Cycle
        nort = nort + 1 
        write(line(iline:),'(a5)') EBS(j); iline=iline+5
        if(iline.gt.75) then
         write(pri,'(a)') trim(line); iline=11
        end if
        z = p(:,j)
        Do ii=1,ns; Do jj=1,ns
         zz(ii,jj) = zz(ii,jj) + z(ii)*z(jj)
        End do; End do        
       End do
       if(nort.eq.0) Cycle

       if(iline.gt.11) write(pri,'(a)') trim(line)

       x = matmul (B,zz)
       y = matmul (zz,B)

       Do jch = 1,nch
        if(ich.gt.jch) then
         ij=ich*(ich-1)/2+jch
         xx = matmul (x,hcc(:,:,ij))
         hcc(:,:,ij) = hcc(:,:,ij) - xx
        elseif(ich.eq.jch) then
         ij=ich*(ich+1)/2
         xx = matmul (x,hcc(:,:,ij))
         hcc(:,:,ij) = hcc(:,:,ij) - xx
         yy = matmul (hcc(:,:,ij),y)
         hcc(:,:,ij) = hcc(:,:,ij) - yy
        else
         ij=jch*(jch-1)/2+ich
         yy = matmul (hcc(:,:,ij),y)
         hcc(:,:,ij) = hcc(:,:,ij) - yy
        end if
       End do

       if(kcp.gt.0) then
        Do i = 1,kcp 
         z = matmul(hcb(:,ich,i),y)
         hcb(:,ich,i)=hcb(:,ich,i)-z
        End do
       end if
       
      End do ! over channels (ich)

      Deallocate(y,xx,yy,zz,z)

      End Subroutine BS_ORTH