!======================================================================
      Subroutine Diag(et)
!======================================================================
! ... diagonalization of Hamiltonian
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer :: il
      Real(8) :: et,t1,t2
      Real(8), external :: Ecore_hf

      Call CPU_time(t1)

      if(icore.gt.0) then 
       Ecore = Ecore_hf()
       write(log,'(/a,f16.8)')  'Ecore  = ',Ecore
      end if 

      Call Update_int

      Call Diag_block(ncfg)

      elevel = elevel + Ecore
      et = SUM(elevel*weight)

      write(log,'(/a/)') ' Level     Energy      Leading CSFs'
      Do il=1,nlevels
       Call print_level(il,ncfg)
      End do

      Call CPU_time(t2)
      time_diag = time_diag + (t2-t1)

      End Subroutine Diag


!======================================================================
      Subroutine  Diag_block(nc)
!======================================================================
! ... diagonalization of the block ib:
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer, intent(in) :: nc
      Integer :: i,j,k,ic,jc,info
      Real(8) :: HM(nc,nc),eval(nc),S

! ... set-up Hamitonian matrix for given block:

      HM = 0.d0

! ... L-intrgrals:

      Do i=1,Lint
       S = L_int(i) 
       Do j = ip_Lint(i-1)+1,ip_Lint(i) 
        ic = ic_Lcoef(j)
        jc = jc_Lcoef(j)
        HM(ic,jc) = HM(ic,jc) + S*L_coef(j)
       End do
      End do

! ... Rk-intrgrals:

      Do i=1,nint
       S = Rk_int(i) 
       Do j = ip_int(i-1)+1,ip_int(i) 
        ic = ic_coef(j)
        jc = jc_coef(j)
        HM(ic,jc) = HM(ic,jc) + S*Rk_coef(j)
       End do
      End do

      if(debug.gt.0) then
       write(log,'(/a)') 'matrix:'
       Do i=1,nc
        write(log,'(10E15.5)') HM(i,1:nc)
       End do
      end if

! ... maximum number of needed solutions

      k = 0
      Do i=1,nlevels
       if(level(i).gt.k) k=level(i)
      End do

      Call LAP_DSYEVX('V','L',nc,nc,HM,eval,k,info)
      if(info.ne.0) Stop "Diag_block: DSYEVX failed"

      Do i=1,nlevels
       elevel(i) = eval(level(i))
       coefs(ip_level(i)+1:ip_level(i)+nc) = HM(1:nc,level(i))
      End do

      End Subroutine Diag_block


!======================================================================
      Subroutine Update_int
!======================================================================
! ... print major contributors to ASF:
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer :: i,k,nu
      Real(8), external :: HLC, rk, hlc_hf 

      Do i=1,Lint
       if(if_Lint(i).eq.1 .and. L_int(i).ne.0.d0) Cycle
       L_int(i) = HLC_hf(i1_Lint(i),i2_Lint(i))

      End do

      Do k = 0,kmax

      Do i=nk_int(k-1)+1,nk_int(k)
       if(if_int(i).eq.1 .and. Rk_int(i).ne.0.d0) Cycle
       Rk_int(i) = rk (i1_int(i),i2_int(i),i3_int(i),i4_int(i),k)
      End do; End do

! ... debug printing of integrals:

      if(debug.gt.0) then
       Call Find_free_unit(nu)
       AF = trim(name)//'.int_value'
       open(nu,file=AF)
       Call Print_int_value(nu) 
       Call Print_int_value(log)                 ! ???
       Close(nu)
      end if

 
      End Subroutine Update_int


!======================================================================
      REAL(8) FUNCTION rk (i1,j1,i2,j2,k)
!======================================================================
!               k
!     Returns  R (i1, j1; i2, j2) base on the assembling the B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. If not, they are calculated by programs
!     mrk_diff or mrk_cell, depending of the parameter 'meth':
!     meth = 'd' - differential equation method
!          = 'c' - cell integration method
!----------------------------------------------------------------------
      Use bsr_mchf
  
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      INTEGER :: i,ip, j,jp
      Real(8) :: a(ns,ks),b(ns,ks)
      Real(8) :: rkj
  
      ! .. check the B-spline integrals in module spline-integrals
  
      if(k.ne.krk.or.itype.ne.'rk') then
       if(meth.eq.'d') then
         Call MRK_diff(k)
       else
         Call MRK_cell(k)
       end if
      end if
  
      ! .. form cross-products
  
      Call density (ns,ks,a,p(1,i1),p(1,i2),'s')
      Call density (ns,ks,b,p(1,j1),p(1,j2),'s')
  
      ! .. assembling the B-spline integrals

      rk = 0.d0
      do ip = 1,ks
        do i = 1,ns-ip+1
          rkj = 0.d0
          do jp = 1,ks
            do j = 1,ns-jp+1
              rkj = rkj+b(j,jp)*rkb(j,i,jp,ip)
            end do
          end do
          rk = rk + a(i,ip)*rkj
        end do
      end do
  
      End Function rk
  

!======================================================================
      Subroutine Print_int_value(nu)
!======================================================================
!     print integrals values
!----------------------------------------------------------------------
      Use bsr_mchf
 
      Implicit none
      Integer :: k,j1,j2,j3,j4, i, nu

      Do i=1,Lint
       j1 = i1_Lint(i); j2 = i2_Lint(i)
       write(nu,'(a,a,a,a,a,a,f12.8)') 'L',' (',ebs(j1),',',ebs(j2),') = ', L_int(i)
      End do

      Do k = kmin,kmax
      Do i = nk_int(k-1)+1,nk_int(k)
       j1 = i1_int(i); j2 = i2_int(i); j3 = i3_int(i); j4 = i4_int(i)
       write(nu,'(a,i2,9a,F12.8)') &
        'R',k,' (',ebs(j1),',',ebs(j2),';',ebs(j3),',',ebs(j4),') = ', Rk_int(i)
      End do; End do

      End Subroutine Print_int_value

