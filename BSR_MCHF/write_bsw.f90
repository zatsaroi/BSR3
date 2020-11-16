!==================================================================
      Subroutine write_bsw
!==================================================================
!     Ouput the radial functions
!------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer :: i

      AF_out = trim(name)//'.bsw'
      Call Read_ipar(inp,'out',AF_out)
      Call Read_iarg('out',AF_out)
      Open(nuw,file=AF_out,form='UNFORMATTED')

      rewind(nuw)
      Do i=1,nbf
       write(nuw) ebs(i), Z, h, hmax, rmax,ks,ns,mbs(i)
       write(nuw) p(:,i)
      End do
      Close(nuw)

      End Subroutine write_bsw


