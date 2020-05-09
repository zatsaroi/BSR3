!======================================================================
      Subroutine SUB_check_orb
!======================================================================
!     analize the orbitals in the file 'AFW' and assign new set indexis
!     if needed
!     orbitals recorded in the "orbitals" list: module orb_LS, and in
!     the "radial" list: module "spline_orbitals" 
!     pointer IEF provides connection between these two lists
!     (two spectroscopic orbitals may have the same radial functions) 
!----------------------------------------------------------------------
      Use bsr_prep

      Implicit real(8) (A-H,O-Z)
      Character(4) :: elw
      Character(4), external :: ELF4
      Integer, allocatable :: IP_occ(:)

      Open(nuw,file=AFW,form='UNFORMATTED')
      rewind(nuw)

!----------------------------------------------------------------------
! ... read radial functions< one by one:

    1 m=nbf+1; if(m.ge.mbf) CALL Allocate_bsorb(mbf+jbf)                                       

      read(nuw,end=2) elw,zw,hw,hmw,rmw,ksw,nsw,mw

      read(nuw) p(1:mw,m); if(mw.lt.ns) p(mw+1:ns,m)=0.d0

      ! check consistence with the B-spline basis:
      if(zw.ne.z)     Stop ' Read_bsw:  z <> zw'
      if(hw.ne.h)     Stop ' Read_bsw:  h <> hw'
      if(abs(hmw-hmax).gt.1.D-12) then
                      write(*,*) elw,hmw,hmax, hmw-hmax
                      Stop ' Read_bsw:  hmw <> hmax'
      end if
      if(rmw.ne.rmax) Stop ' Read_bsw:  rmw <> rmax'
      if(ksw.ne.ks)   Stop ' Read_bsw:  ksw <> ks'
      if(nsw.ne.ns)   Stop ' Read_bsw:  nsw <> ns'

! ... find orbital in orbital list:

       Call EL4_nlk(elw,n,l,k);  k=k+kshift
       ii=Ifind_nlk(n,l,k,0)

       if(ii.gt.0) then                              
        mbs(m)=mw; nbs(m)=n; lbs(m)=l; kbs(m)=k; ebs(m)=elw
       else
        write(pri,'(a4,12x,a)') elw,' excessive orbital'
        go to 1
       end if

! ...  define overlaps with existing orbitals:

       OBS(:,m) = 0.d0
       Do i = 1,m
        if(lbs(i).ne.lbs(m)) Cycle
        OBS(i,m)=QUADR(i,m,0)
        OBS(m,i)=OBS(i,m)
       End do

!----------------------------------------------------------------------
! ...  compare with the existing orbitals:
       
       SM1 = QUADR(m,m, 1); SM2 = QUADR(m,m, 2)
       Do i = 1,nbf;     if(abs(OBS(i,m)).lt.eps_ovl) Cycle

        ! ... check orthogonality to core: 

        if(i.le.nclosd.and.ii.gt.nclosd) then
         write(pri,'(a,a,a,a,a,f12.8)') ' file ', trim(AFC), '  orbital ', &
               elw, ' does not orthogonal to core:', OBS(i,m)
         Stop 'problem with orthogonality to core'
        end if

        ! ... define if it is approximately the same orbital:

        S  = abs(OBS(i,m)-1.d0);      if(S .gt.eps_ovl) Cycle
        S1 = abs(QUADR(i,i, 1)-SM1);  if(S1.gt.eps_ovl) Cycle
        S2 = abs(QUADR(i,i, 2)-SM2);  if(S2.gt.eps_ovl) Cycle

        IEF(ii) = i
        write(pri,'(a,a,a,a)')  elw,' --> ',ebs(i),'    the same'
        go to 1
       
       End do

!---------------------------------------------------------------------
! ...  core orbitals should be the same:   

       if(ii.le.nclosd) then
        write(pri,'(a,a,a)') 'file ',AFW,'  has another core orbital'
        Stop ' another core orbital? '
       end if

! ... assign set index for new radial orbital: 

       Call Assign_index(m); nbf=m; IEF(ii)=m 

       go to 1    ! go to next orbital
   2  Continue

! ... check if bsw-file contains all radial orbitals:

      Do i=ncore+1,nwf; ii=IEF(i)
       if(ii.eq.0) then
        write(pri,'(a,a,a,a)') &
        trim(AFC), '-  orbital ',ELF(i),' was not found in the w-file'
        Stop ' unknown orbitals ! '
       else
        NEF(i)=nbs(ii); LEF(i)=lbs(ii); KEF(i)=kbs(ii); ELF(i)=ebs(ii)
       end if
      End do

      End Subroutine SUB_check_orb
