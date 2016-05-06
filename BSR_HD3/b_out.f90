!=======================================================================
      Subroutine b_out 
!=======================================================================
!     output the bound solutions
!-----------------------------------------------------------------------
      Use bsr_hd
      Use target, only: nelc,nz,Etarg,coupling
      Use channel
      Use conf_LS
      Use spline_param, only: ns
      
      Implicit none
      Character(64) :: Labl
      Real(8) :: Ebind(khm), eff_n(khm)
      Real(8) :: E1,E2,zion
      Integer :: n_eff(khm)
      Integer :: i,j,i1,i2,j1,j2,ich,is,js,it,nbound,ms

! ... output file:

      i = INDEX(AF_b,'.',BACK=.TRUE.); AF = AF_b(1:i)//ALSP
      Open(nub,file=AF)
      i = INDEX(AF_ub,'.',BACK=.TRUE.); AF = AF_ub(1:i)//ALSP
      Open(nuu,file=AF,form='UNFORMATTED')

! ... local allocations:

      if(allocated(v)) Deallocate(v);  Allocate(v(1:nhm))

      Ebind = 0.d0; eff_n = 0.d0; n_eff = -1

      zion = nz - nelc; if(zion.eq.0.d0) zion = 1.d0
!----------------------------------------------------------------------
! ... define number of bound states for output:

      E1=Emin; ! if(E1.eq.0.d0) E1 = Etarg(1)*1.5
      E2=Emax; ! if(E2.eq.0.d0) E2 = -0.1
      if(IT_max.gt.0) E2 = Etarg(IT_max)

      ms=msol; if(msol.le.0.or.msol.gt.khm) ms=khm

      nbound = 0
      Do is = 1,ms
       if(eval(is).gt.E2.and.E2.ne.0.d0) Cycle
       if(eval(is).lt.E1.and.E1.ne.0.d0) Cycle
       ich = isol(is); it = 0
       if(ich.le.kch) then
        it=iptar(ich); Ebind(is)=eval(is)-Etarg(it)  

        ! ... define 'n' for outer electron: 

        if(Ebind(is).gt.0.d0) then
         n_eff(is) = ICHAR('k')-ICHAR('1')+1
         eff_n(is) = sqrt(2*Ebind(is)) 
        else
         eff_n(is) = zion / sqrt(abs(2*Ebind(is))) 
         n_eff(is) = ICHAR('n')-ICHAR('1')+1
        end if
       else  ! perturber
        n_eff(is) = -2
       end if
       nbound = nbound + 1
      End do        
        
!----------------------------------------------------------------------
! ... store the solutions:

      write(nub,'(5i10,3i5,T74,a)') &
                  ns,kch,kcp,nhm,nbound,lpar,ispar,ipar, &
              '=> ns,kch,kcp,nhm,nbound,lpar,ispar,ipar'
      write(nuu) ns,kch,kcp,nhm,nbound,lpar,ispar,ipar

      js = 0
      Do is = 1,ms
       if(n_eff(is).eq.-1) Cycle
       js = js + 1

       ! ... define label:

       Call Find_channel_label(isol(is),1,is,eval(is),LabL)

       write(nub,'(i5,2x,a)') js,trim(LABL)
       write(nub,'(E20.10,f15.5,f10.2,T74,a)') &
             eval(is),Ebind(is)*au_eV,eff_n(is), &
             '=> E(au), E_bind(eV), n_effective '          
       write(nuu) js,LABL
       write(nuu) eval(is),Ebind(is)*au_eV,eff_n(is)

       ! ...  define solution in original B-spline basis:

       v = 0.d0
       Do ich = 1,kch; i1=(ich-1)*ns+1; i2=ich*ns
        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
        Do j=j1,j2
         v(i1:i2) = v(i1:i2) + a(j,is)*bb(1:ns,j)     
        End do
       End do
       if(kcp.gt.0) v(kch*ns+1:nhm)=a(ksol+1:khm,is)

       write(nub,'(5D15.8)') v(1:nhm)
       write(nuu) v(1:nhm)

      End do

! ... energies output:

      write(nub,*)
      write(nub,'(a,i8)')  'khm =',khm
      write(nub,*)
      Do i=1,khm
       write(nub,'(2D16.8)') eval(i), (eval(i)-etarg(1))*2
      End do

      Do i=1,khm
       write(nuu) eval(i), (eval(i)-etarg(1))*2
      End do

      Close(nub); Close(nuu)

      End Subroutine b_out



