!======================================================================
     Subroutine get_spline_param 
!======================================================================
!    This routine gets information about the spline parameters and
!    and set up spline knots and all relevant arrays.
!    Program checks file AF_grid (knot.dat or name.knot)
!    for B-spline parameters  or use default ones.
!    The spline paremeters can also be modified from command line. 
!----------------------------------------------------------------------
     Use bsr_hf

     Implicit none
     Logical :: EX

! .. define the knot sequence:

     AF = ' ';  Call Read_apar(inp,'knot',AF) 
                Call Read_aarg('knot',AF)
     if(len_trim(AF).ne.0) then
      AF_grid=trim(AF)
     else
      Inquire(file=AF_grid,exist=EX)
      if(.not.EX) AF_grid=trim(name)//trim(BF_grid)
     end if

     Call Read_iarg('ks',ks)
     Call Read_rarg('h',h)                                       
     Call Read_rarg('hmax',hmax)
     Call Read_rarg('rmax',rmax)

     Call define_grid(z)

! .. initialize other B-spline arrays:

     Call define_spline

     Call allocate_RK_integrals

     write(log,'(/a/)') 'B-spline parameters:'
     write(log,'(a,i5,t25,a)')   'ns   =',ns,  '-  number of splines' 
     write(log,'(a,i5,t25,a)')   'ks   =',ks,  '-  order  of splines' 
     write(log,'(a,f9.3,t25,a)') 'h    =',h ,  '-  initial point' 
     write(log,'(a,f9.3,t25,a)') 'hmax =',hmax,'-  exponetial step size' 
     write(log,'(a,f9.3,t25,a)') 'rmax =',rmax,'-  maximum radius' 

     End Subroutine get_spline_param 

