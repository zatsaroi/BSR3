!=======================================================================
      Subroutine b_out 
!=======================================================================
!     output the bound solutions
!-----------------------------------------------------------------------
      USE bsr_hd
      USE blacs
      USE target, only: nelc,nz,Etarg,coupling
      USE channel
      USE conf_LS
      USE spline_param, only: ns
      
      Implicit none
      Character(64) ::  Labl
      Real(8), Allocatable :: Ebind(:), eff_n(:), vb(:)
      Real(8) :: E1,E2,zion
      Integer, Allocatable :: n_eff(:)
      Integer :: i,j,i1,i2,j1,j2,ich,is,js,it,nbound,ms

      Call CPU_time(t0)     

      if(io_processor) then 

! ... output file:

      i = INDEX(AF_bound,'.',BACK=.TRUE.); AF = AF_bound(1:i)//ALSP
      Open(nub,file=AF)

      Allocate(Ebind(khm),eff_n(khm),n_eff(khm),vb(nhm))
      Ebind = 0.d0; eff_n = 0.d0; n_eff = -1

      zion = nz - nelc; if(zion.le.0.d0) zion = 1.d0

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
         n_eff(is) = ICHAR('n')-ICHAR('1')+1
         eff_n(is) = zion/sqrt(abs(2*Ebind(is))) 
        end if
       else  ! perturber
        n_eff(is) = -2
       end if
       nbound = nbound + 1
      End do        

      end if  ! io_processor

!----------------------------------------------------------------------
! ... broadcast some parameters:

      Call br_ipar(ms)
      if(io_processor) then
       Call igebs2d (ctxt, 'all', ' ', khm, 1, n_eff, khm)
      else
       if(allocated(n_eff)) deallocate(n_eff); allocate(n_eff(khm) )
       Call igebr2d (ctxt, 'all', ' ', khm, 1, n_eff, khm, rsrc, csrc)
      end if    

!----------------------------------------------------------------------
!                                                  store the solutions:
       if(io_processor) &
       write(nub,'(5i10,3i5,T74,a)') &
                  ns,kch,kcp,nhm,nbound,lpar,ispar,ipar, &
              '=> ns,kch,kcp,nhm,nbound,lpar,ispar,ipar'

       js = 0
       Do is = 1,ms
        call BLACS_BARRIER (ctxt, 'all')
        if(n_eff(is).eq.-1) Cycle

        call pdgeadd ('notrans', khm, 1, one, z, 1,is, descz, &
                                        zero, v, 1, 1, descv)
        call BLACS_BARRIER (ctxt, 'all')

        if(.not.io_processor) Cycle      

        vb = 0.d0
        Do ich = 1,kch; i1=(ich-1)*ns+1; i2=ich*ns
         j1 = ipsol(ich-1)+1; j2=ipsol(ich)
         Do j=j1,j2
          vb(i1:i2) = vb(i1:i2) + v(j)*bb(1:ns,j)     
         End do
        End do
        if(kcp.gt.0) vb(kch*ns+1:nhm)=v(ksol+1:khm)

        js = js + 1

        ! ... define label:

        Call Find_channel_label(isol(is),1,is,eval(is),LabL)

        write(nub,'(i5,2x,a)') js,trim(LABL)
        write(nub,'(E20.10,f15.5,f10.2,T74,a)') &
             eval(is),Ebind(is)*au_eV,eff_n(is), &
             '=> E(au), E_bind(eV), n_effective '          

        write(nub,'(5D15.8)') vb(1:nhm)

       End do

       if(io_processor) then

       ! ... energies output:

        write(*,'(/a,i10)') 'khm =',khm

        write(nub,*)
        write(nub,'(a,i10)')  'khm =',khm
        write(nub,*)
        Do i=1,khm
         write(nub,'(2D16.8)') eval(i), (eval(i)-etarg(1))*2
        End do

        Close(nub)
        Deallocate(Ebind, eff_n, n_eff, vb)

        write(*,'(/a)') 'b_out: bound.nnn is recorded'

       end if

       if(io_processor) then           
        Call CPU_time(t1)
        write (pri,'(/a,T30,f10.2,a)') 'B_out:,', (t1-t0)/60, ' min.'
        write (*  ,'(/a,T30,f10.2,a)') 'B_out:,', (t1-t0)/60, ' min.'
       end if

       End Subroutine b_out



