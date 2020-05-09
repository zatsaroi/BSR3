!======================================================================
      SUBROUTINE mrk_cell(k)
!======================================================================
!
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm
!
!     Calls: rk_moments
!
!----------------------------------------------------------------------
!
!     on entry      k        multipole index
!     --------
!       
!     on exit       rkb     four-dimensional array of Slater integrals 
!     -------               of power k in the B-spline basis
!                           (in module spline-integrals)
!----------------------------------------------------------------------

      USE spline_param
      USE spline_integrals
      USE spline_moments
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: k
   
      ! .. local variables
   
      INTEGER(4) :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      REAL(8) :: c
   
      ! .. check the need of calculations
   
      if(itype == 'aaa') Call allocate_integrals
      if(itype == 'rk ' .and. krk == k) Return
   
      ! .. compute the moments in the spline basis
   
      CALL rk_moments(k)
   
      ! .. generate the rkb array
   
      rkb=0.d0
   
      DO jv=1,nv
       jj = 0
       DO jh=1,ks
        j = jv + jh - 1
        DO jhp=jh,ks
         jp = jhp - jh + 1
         jj = jj + 1
   
         DO iv=1,nv
          ii = 0
          DO ih=1,ks
           i = iv + ih -1
           DO ihp=ih,ks
            ip = ihp - ih + 1
            ii = ii + 1
   
            if( iv < jv ) then
             c = rkd1(ii,iv)*rkd2(jj,jv)
            else if( iv > jv ) then
             c = rkd1(jj,jv)*rkd2(ii,iv)
            else
             c = rkd(ii,jj,iv) 
            end if
         
            rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c 
          
           END DO
          END DO
         END DO

        END DO
       END DO
      END DO
   
      itype='rk '
      krk=k
   
      END SUBROUTINE mrk_cell
   
   