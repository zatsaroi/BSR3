!=====================================================================
      Subroutine solve_HF
!=====================================================================
!     Solve the HF equations in turns 
!---------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals

      Implicit none
      Real(8) :: hfm(ns,ns), v(ns), et
      Real(8), external :: QUADR_hf, BVMV
      Integer :: ip, i,j, it

      Real(8) :: t1,t3,t4, S,S1,S2
    
      Call CPU_time(t1)

      et = etotal
      dpm = 0.d0
      Do it = 1,max_it
       write(log,'(//A,I6/A/)') 'Iteration ',it, &
                                '----------------'
       write(log,'(2x,a,9x,a,10x,a,9x,a/)')  'nl', 'e(nl)', 'dpm', 'ns'

! ... main iterations other orbitals:

       Do ip = 1,nwf; i=iord(ip); if(i.eq.0) Cycle

        Do j=i+1,nwf
         if(it.le.2) Cycle; if(rotate.eq.0) Cycle;  Call Rotate_ij(i,j)
        End do

        Call hf_matrix(i,hfm)

        ! .. diagonalize the hf matrix

!        if(it <= 2 .or. newton == 0  .or. dpm(i) > nr_tol ) then
         Call hf_eiv (i,hfm,v) 
!        else
!         Call hf_nr (i,hfm,v,e(i,i))
!        end if

        S = maxval(abs(p(1:ns,i)-v(1:ns)))/maxval(abs(p(:,i)))
        if(ip.eq.1) then
         if(S.lt.orb_tol) iord(ip)=0
        else
         if(S.lt.orb_tol.and.iord(ip-1).eq.0) iord(ip)=0
        end if

        if(it.gt.1.and.(S.gt.dpm(i).or.acc.eq.1)) then
         v(1:ns) = aweight * v(1:ns) + bweight * p(1:ns,i)
         S1 = BVMV (ns,ks,sb,'s',v,v); S2 = sqrt(S1)
         v = v / S2
        end if
       
        dpm(i)=maxval(abs(p(1:ns,i)-v(1:ns)))/maxval(abs(p(:,i)))

        write(log,'(A5,f16.8,1PD13.2,I8)')  ebs(i),e(i,i),dpm(i),mbs(i)

        p(1:ns,i) = v(1:ns)

        Call Check_tails(i)

        Call CPU_time(t3)
        Call Update_int(i)
        Call CPU_time(t4)
        time_update_int=time_update_int+t4-t3

        ! .. remove tail zero

        if(debug.gt.0) then
         Do j = 1,nwf 
          if(i.eq.j) Cycle 
          if(lbs(i).ne.lbs(j)) Cycle
          write(log,'(a,a,a,f16.8)') &
           'Orthogonality ',ebs(i),ebs(j),QUADR_hf(p(1,i),p(1,j),0)
         End do
        end if

       End do ! over functions 

       Call Energy
       write(log,'(/a,F16.8)') 'Etotal = ',etotal

       orb_diff =  maxval(abs(dpm(1:nwf)))
       scf_diff =  abs(et-etotal)/abs(etotal)
       et = etotal
       write(log,*)
       write(log,'(A,T25,1PD10.2)') 'Maximum orbital diff. ', orb_diff
       write(log,'(A,T25,1PD10.2)') 'Energy difference     ', scf_diff
       write(scr,'(a,F20.12,i5,1P2E10.2)') &
        'etotal,it,E_diff,orb_diff = ',etotal,it,scf_diff,orb_diff
       
       Call Boundary_conditions

       if ( orb_diff < orb_tol .and. scf_diff  < scf_tol ) Exit

      End do ! over iterations

      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'time_eiv:',time_hf_eiv,'  sec' 
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'time_mat:',time_hf_matrix,'  sec' 
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'time_int:',time_update_int,'  sec' 

      End Subroutine solve_HF


!==================================================================
      Subroutine hf_eiv(i,hfm,v)
!==================================================================
!     Find the eigenvector v of hfm for the "positive-eneregy"
!     eigenvalue m=n-l. Each orthogonality condition to the lower 
!     state reduces m - 1.  We supposed that nl orbitals are
!     ordered to their n-values.
!------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals

      Implicit none
      Integer, intent(in)    :: i
      Real(8), intent(inout) :: hfm(ns,ns)
      Real(8), intent(out)   :: v(ns)
      Real(8) :: aa(ns,ns),ss(ns,ns),w(3*ns),eval(ns),a(ns),s(ns),t1,t2
      Integer :: j, jp, info, k,kk,m, ipos(1)
    
      Call CPU_time(t1)

! ... apply orthogonality conditions for orbitals  ??? to all n'l with n'<n  ???

      m = nbs(i)-lbs(i)
      Do j = 1,nwf 
       if(i.eq.j) Cycle 
       if(lbs(i).ne.lbs(j)) Cycle
       if(j.lt.i) m = m - 1
       Call orthogonality(hfm,p(1,j))
       if(debug.gt.0) write(log,'(a,a,a,a,i3)') &
        'Orthogonality ',ebs(i),ebs(j),' is applied, m =',m
      End do

! ... apply boundary conditions (delete extra B-splines)

      kk=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle; kk=kk+1
       k=0
       Do jp=1,ns
        if(iprm(jp,i).eq.0) Cycle
        k=k+1; a(k)=hfm(jp,j); s(k)=bb(jp,j)  
       End do
       aa(1:k,kk)=a(1:k); ss(1:k,kk)=s(1:k)
      End do 

! ... evaluates the eigenvalues and eigenvectors (LAPACK routine):

      Call dsygv(1,'V','L',kk,aa,ns,ss,ns,eval,w,3*ns,INFO)
      if (info /= 0) then
       write(scr,'(a,i6)') 'Error in Eigenvalue routine, dsygv', info
       Stop ' '
      end if

! ... restore the solutions in original B-spline net:

      a(1:ns) = aa(1:ns,m);  v=0.d0; k=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
      End do 

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

!      if (v(ks) < 0.d0) v = -v

      e(i,i) = eval(m)

!      if(debug.gt.0) then
!       write(log,'(a,5E15.5)') 'eval =',eval(mm-2:mm+2)
!       write(log,'(a,2i5,E15.5)') 'we choose m,mm,e =',m,mm,eval(mm)
!      end if

      Call CPU_time(t2)
      time_hf_eiv = time_hf_eiv + t2-t1

      End Subroutine hf_eiv

