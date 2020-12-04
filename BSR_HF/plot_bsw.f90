!======================================================================
      Subroutine plot_bsw
!======================================================================
! ... Computes and tabular the radial orbitals in all gausian points, 
! ... plus border values, for further plots
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals
 
      Implicit none
      Real(8) :: yp(nv*ks+2,nbf),r(nv*ks+2)
      Integer :: i,j,io,i1,i2,m
 
      if(out_plot.eq.0) Return
      AF_plt = trim(name)//BF_plt   
      Call Read_apar(inp,'plot',AF_plt)
      Call Read_aarg('plot',AF_plt)
      open(nup,file=AF_plt)
 
! ... radial points for output:

      m=1; r(1)=t(1)
      Do i=1,nv; Do j=1,ks; m=m+1; R(m)=gr(i,j); End do; End do
      m=m+1; R(m)=t(ns+1)
 
      yp = 0.d0
      Do io=1,nbf
       Call Bvalue_bm(p(1,io),yp(1,io))
      End do 

      rewind(nup)
      i2=nbf; i1=max(1,i2-nit+1)
      write(nup,'(100(6x,a6))') 'r',(ebs(i),i=i1,i2)
      Do i=1,m
       write(nup,'(100(1PE12.4))') r(i),(yp(i,io),io=i1,i2)
      End do

      End Subroutine plot_bsw


!======================================================================
      Subroutine BVALUE_bm(a,y)
!======================================================================
!     Computes the function y(r) = SUM(i)  a_i * B_i(r)
!     in all gausian points + border values
!----------------------------------------------------------------------
      Use bsr_hf

      Implicit none
      Integer :: i,j,iv,ith,m
      Real(8), intent(in)  :: a(ns)
      Real(8), intent(out) :: y(nv*ks+2)

      y = 0.d0
      Do iv = 1,nv                 ! over intervals
       Do m = 1,ks                 ! over gausian points
        Do ith = 1,ks              ! over B-splines in given interval
         i = iv+ith-1              ! B-spline index
         j = 1 + (iv-1)*ks + m     ! radial point index
         y(j) = y(j) + a(i)*bsp(iv,m,ith)
        End do
       End do
      End do

      m = nv*ks+2                  ! last point
      Do i=1,ks
       y(m) = y(m) + bsp(nv+1,1,i) * a(ns-ks+i)
      End do

      m = 1                        ! first point
      Do i=1,ks
       y(m) = y(m) + bsp(nv+1,2,i) * a(i)
      End do

      End Subroutine BVALUE_bm
