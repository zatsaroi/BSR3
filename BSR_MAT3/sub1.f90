!======================================================================
      Subroutine SUB1
!======================================================================
!     provides calculations for one partial wave
!----------------------------------------------------------------------
      Use bsr_mat
      Use cmdata
      Use bsr_matrix

      Use target,          only: nelc,Etarg, coupling
      Use channel,         only: nch,ncp, ipconf, iptar, jpar,npert

      Use spline_atomic,   only: EC,kclosd
      Use spline_orbitals, only: nbf,lbs
      Use spline_param,    only: ns,ks      
      Use spline_galerkin, only: sb
      Use spline_grid,     only: t

      Implicit none
      Real(8) :: C,t3,t4
      Integer :: i,j,k,l,ich,jch, ij, i1,i2
      Character(100) :: line

!-----------------------------------------------------------------------
! ... read configuration expansion and orbitals information: 

      Call Read_data

!-----------------------------------------------------------------------
! ... initialize arrays:

      Call Allocate_matrix(nch,ns,npert,mk,ipconf(nch))

      write(pri,'(/a/  )') 'Main dimensions in bsr_matrix module:'
      write(pri,'( a,i6,a)') 'kch    = ',kch, &
       '  -  number of channels '
      write(pri,'( a,i6,a)') 'kcp    = ',kcp, &
       '  -  number of perturbers '
      write(pri,'( a,i6,a)') 'kns    = ',kns, &
       '  -  number of splines ' 
      i = kch*kns+kcp; C = i*(i+1)/2; C = C /(128*1024)
      write(pri,'(/a,T33,i8)') 'Interaction matrix dimension:',i
      write(pri,'(/a,T33,f8.1,a)') 'Reqired memory for bsr-matrix:',C,' Mb'

      l=maxval(lbs(1:nbf));  npol=max(l,mk); k=(npol+2)*ntype*2
      if(nb.lt.k) then
       nb=k
       write(pri,'(/a,i8,a)') &
       'nb = ',nb,'  -  number of blocks re-assigned !!! '
      end if

      i = nb*(4+mb);  C = 28*i; C = C /(1024*1024)
      write(pri,'(/a,T33,f8.1,a)') 'Reqired memory for cmdata:',C,' Mb'

      Call Allocate_cmdata 

      Call Allocate_ndets(-1)  
      Call Allocate_ndefs(-1)  

      if(.not.allocated(CBUF)) &
       Allocate(CBUF(maxnc),ijtb(maxnc),intb(maxnc),idfb(maxnc))

       C = 20*maxnc; C = C /(1024*1024)
      write(pri,'(/a,T33,f8.1,a)') 'Reqired memory for bufer:',C,' Mb'

!-----------------------------------------------------------------------
!                                                        overlap matrix:
      Call CPU_time(t3)

! ... B-spline overlaps:

      Do ich=1,nch;  Call UPDATE_HL(ich,ich,ns,ks,sb,1.d0); End do

! ... read data from INT.BNK ...

      icase = 11;  kpol = 0;  Call State_res

! ... symmetrize the int. matrix ...

      Do ich = 1,nch; ij=ich*(ich+1)/2
      Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij);  hcc(i,j,ij)=C/2.d0;  hcc(j,i,ij)=C/2.d0
      End do; End do
      End do

! ... check big overlaps:

      Call Check_mat(k)

! ... redo overlap matrix:

      if(k.gt.0) then

       Call Allocate_matrix(nch,ns,npert,mk,ipconf(nch))
       Do ich=1,nch;  Call UPDATE_HL(ich,ich,ns,ks,sb,1.d0); End do
       icase = 11;  kpol = 0;  Call State_res
       Do ich = 1,nch; ij=ich*(ich+1)/2
       Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij);  hcc(i,j,ij)=C/2.d0;  hcc(j,i,ij)=C/2.d0
       End do; End do
       End do

      end if

! ... record overlap matrix: 

      rewind(nui)
      write(nui) ns,nch,kcp

      Call Record_matrix(nui)  

      Call CPU_time(t4)
      write(pri,'(/a,f10.2,a)') 'Overlaps:     ',(t4-t3)/60,' min '
      write(  *,'( a,f10.2,a)') 'Overlaps:     ',(t4-t3)/60,' min '

!----------------------------------------------------------------------
! ... Interaction matrix:
!----------------------------------------------------------------------
! ... core-energy shift

      hcc = hcc * EC       
      if(kcp.gt.0) then;  hcb = hcb * EC; hbb = hbb * EC;  end if
!----------------------------------------------------------------------
!                                                          L-integrals:
      Call CPU_time(t3)

      Call Gen_Lval

      icase=6; Call State_res   

      Call Alloc_Lcore(0,0,0)

      Call CPU_time(t4)
      write(pri,'(/a,f10.2,a)') 'L-integrals:  ',(t4-t3)/60,' min '
      write(*,'(a,f10.2,a)') 'L-integrals:  ',(t4-t3)/60,' min '


!----------------------------------------------------------------------
!                                                          Z-integrals:
      if(mso.gt.0) then    

       Call CPU_time(t3)

       Call Gen_Zval;  icase=7;  Call State_res   

       Call Alloc_zcore(0,0,0)

       Call CPU_time(t4)
       
       write(*,'(a,f10.2,a)') 'Z-integrals:  ',(t4-t3)/60,' min '

       write(pri,'(/a,f10.2,a)') 'Z-integrals:  ',(t4-t3)/60,' min '

      end if
!----------------------------------------------------------------------
!                                                          R-integrals:

      Do icase = 3,10

       Select case(icase)
        case(6,7); Cycle
        case(8,9); if(msoo.eq.-1.or.(mrel.lt.3.and.msoo.ne.1)) Cycle
        case( 10); if(mss.eq.-1.or.(mrel.lt.4.and.mss.ne.1)) Cycle
        case(3,4); if(moo.eq.-1.or.(mrel.lt.5.and.moo.ne.1)) Cycle
       End Select

        Call CPU_time(t3);   Call State_res;   Call CPU_time(t4)

        write(pri,'(/a,a,f10.2,a)') Aint(icase),'-integrals:  ',&
                                      (t4-t3)/60,' min '
        write(*,'(a,a,f10.2,a)') Aint(icase),'-integrals:  ',&
                                    (t4-t3)/60,' min '

      End do  ! over icase

!----------------------------------------------------------------------
!                                                      target energies:
      Do i = 1,nch
       j=iptar(i); C = Etarg(j)-EC
       Call UPDATE_HL(i,i,ns,ks,sb,C)
      End do

!----------------------------------------------------------------------
!                                            output interaction matrix:
! ... symmetrize the diagonal blocks:

      Do ich = 1,nch; ij=ich*(ich+1)/2
      Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij);  hcc(i,j,ij)=C/2.d0;  hcc(j,i,ij)=C/2.d0
      End do; End do
      End do

! ... orthogonal conditions:

      Call CPU_time(t3);   Call BS_ORTH;  Call CPU_time(t4)

      write(pri,'(/a,4x,f10.2,a)') 'BS_ORTH:  ',(t4-t3)/60,' min '
      write(*  ,'( a,4x,f10.2,a)') 'BS_ORTH:  ',(t4-t3)/60,' min '

! ... record interaction matrix: 

      Call Record_matrix(nui)

!----------------------------------------------------------------------
! ... asymptotic coefficients:  

! ... Symmetrize the ACF - matrix and add the correction
! ... from core screening for k=0:

      eps_acf = 1.d-5
      Do k = 0,mk
       Do i = 1,nch
        Do j = i,nch
         C = ACF(i,j,k) + ACF(j,i,k); if(i.ne.j) C=C*2.d0
         if(abs(C).lt.eps_acf) C=0.d0
         ACF(i,j,k) = C; ACF(j,i,k) = C
        End do
       End do
      End do

      k=0; Do i=1,kclosd; k=k+2*(4*lbs(i)+2); End do
      Do i = 1,nch;  ACF(i,i,0)=ACF(i,i,0)+k; End do

      write(pri,'(/a,i2)') 'Asymptotic coefficients: mk = ',mk

      k=0
      Do i = 1,nch
       if(abs(ACF(i,i,0)-2*nelc).lt.eps_acf) Cycle
       if(k.eq.0) then
       write(pri,'(/a,i3/)') &
        'For k=0, afc(i,i) should be equal to 2*nelectrons =', 2*nelc
       write(pri,'(a,1Pe10.1/)') 'k=0 deviations if > eps_acf = ', eps_acf
       end if
       k=k+1
       write(pri,'(i5,2F15.6)') i,ACF(i,i,0)-2*nelc
      End do
      
      if(pri_ac.gt.0) then
      write(pri,'(/a/)') 'Asymptotic coefficients: i,j, ACF(i,j,k)'
      Do k=0,mk
       if(SUM(acf(:,:,k)).eq.0) Cycle
       write(pri,'(a,i2)') 'k = ',k
       ij = 0
       Do i=1,nch; Do j = 1,i      
        if(abs(acf(i,j,k)).lt.eps_acf) Cycle
        i1=ij*20+1; i2=i1+19 
        write(line(i1:i2),'(2i4,E12.3)') j,i,acf(i,j,k)
        ij=ij+1
        if(ij.lt.5) Cycle
        write(pri,'(a)') line; ij=0
       End do; End do
       if(ij.eq.0) Cycle
       i1=1; i2=ij*20
       write(pri,'(a)') line(i1:i2)
      End do
      end if

      write(nui) mk;  write(nui) ACF;  write(nui) t(ns+1),ns

      if(pri_f.ne.0) Call f_values

      if(iitar.ne.0) Call Target_new

!----------------------------------------------------------------------
!                                            target interaction matrix:
      write(nui) htarg
      write(nui) otarg
      write(nui) Etarg
      write(nui) EC
      close(nui)

      Call Target_print(pri,EC,eps_tar)
      
      End Subroutine SUB1 



