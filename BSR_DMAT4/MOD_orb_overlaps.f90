!======================================================================
      Module orb_overlaps
!======================================================================
!     contains the desription of one-electron overlaps 
!----------------------------------------------------------------------
      Implicit none

! ... number of orbitals:

      Integer :: norb = 0

! ... pointer to orbitals in l-order and copy of l-values:

      Integer, allocatable :: ipl(:),jpl(:), lorb(:), chan(:) 

! ... base pointer to orbitals for given l

      Integer :: max_l = 0
      Integer, allocatable :: ip_l(:), jp_l(:)
      
! ... base pointer to overlaps for given l

      Integer, allocatable :: ip_ovl(:)

! ... curent number of overlaps

      Integer :: nobs = 0      

! ... value of the bound-bound one-electron overlap
! ... or orthogonal conditions for continuum orbitals

      Real(8), allocatable :: Cobs(:)

! ... memory in words (4 bytes):

      Integer :: mem_orb_overlaps = 0

      Integer :: ibo = 2**15 

      End Module orb_overlaps

!======================================================================
      Subroutine Alloc_orb_overlaps(nbf,lbs,iech,ncore)
!======================================================================
!     allocate, de-allocate or re-allocate arrays in given module
!----------------------------------------------------------------------
      Use orb_overlaps

      Implicit none
      Integer, intent(in) :: nbf,lbs(nbf),iech(nbf),ncore
      Real(8), external :: QUADR
      Integer :: i,j,l, il,jl,is,js, ip 

      if(allocated(ipl)) Deallocate(ipl,jpl,lorb,chan,ip_l,jp_l,ip_ovl,Cobs)
      norb = 0; mem_orb_overlaps = 0
      if(nbf.le.0) Return

! ... sort the orbitals according l-values:
  
      norb = nbf
      Allocate(lorb(norb)); lorb(1:norb) = lbs(1:norb)
      Allocate(chan(norb)); chan(1:norb) = iech(1:norb)
      Allocate(ipl(norb),jpl(norb))
      Call SORTI(norb,lorb,ipl)

      Do i=1,norb; j = ipl(i); jpl(j) = i; End do


! ... find last entry for each l:

      max_l  = maxval(lorb)
      Allocate(ip_l(0:max_l),jp_l(0:max_l));  jp_l=0
      Do j = 1, norb; i = ipl(j)
       l=lorb(i); jp_l(l) = j
      End do

! ... shift ip_l to get "base" for given l:     

      Do l=1,max_l;  if(jp_l(l).eq.0) jp_l(l)=jp_l(l-1); End do
      Do l=max_l,1,-1;  ip_l(l) = jp_l(l-1); End do;  ip_l(0)=0

! ... allocate the overlaps:

      Allocate(ip_ovl(0:max_l)); ip_ovl=0
      i = jp_l(0)-ip_l(0); nobs = i*(i+1)/2  
      Do l = 1,max_l
       ip_ovl(l) = nobs
       i = jp_l(l)-ip_l(l)
       nobs = nobs + i*(i+1)/2  
      End do
      Allocate(Cobs(nobs)); Cobs = 0.d0
      mem_orb_overlaps = 2*nobs + 4*norb + 3*max_l

! ... radial overlaps

      Do il = 1,norb;  i = ipl(il);  l = lorb(i)
      Do jl = 1,il;    j = ipl(jl);  if(l.ne.lorb(j)) Cycle              
       ip = ip_l(l)

       is = max(il,jl)-ip; js = min(il,jl)-ip
       ip = ip_ovl(l) + is*(is-1)/2 + js

       if(chan(i).eq.0.and.chan(j).eq.0) then
        Cobs(ip) = QUADR(i,j,0)
       elseif(chan(i).ne.0.and.chan(j).eq.0) then
        Cobs(ip) = min(i,j)*ibo+max(i,j) ! i*ibo+j        
        if(j.le.ncore) Cobs(ip) = 0.d0
       elseif(chan(i).eq.0.and.chan(j).ne.0) then
        Cobs(ip) = min(i,j)*ibo+max(i,j) ! j*ibo+i        
        if(i.le.ncore) Cobs(ip) = 0.d0
       else
        Cobs(ip) = min(i,j)*ibo+max(i,j)  ! i*ibo+j        
       end if

      End do
      End do

      Call Check_orb_overlaps

      End Subroutine Alloc_orb_overlaps


!======================================================================
      Real(8) Function OBS(io,jo)
!======================================================================
!     find overlaps <io|jo>
!----------------------------------------------------------------------
      Use orb_overlaps

      Implicit none
      Integer, intent(in) :: io,jo
      Integer :: i,j,l, ip, is,js

! ... check the arguments:

      obs = 0.d0; if(lorb(io).ne.lorb(jo)) Return;  l=lorb(io)

      i = jpl(io);  j = jpl(jo)
      ip = ip_l(l)
      is = max(i,j)-ip; js = min(i,j)-ip
      ip = ip_ovl(l) + is*(is-1)/2 + js

      obs = Cobs(ip)

      End Function OBS



!======================================================================
      Real(8) Function VDET (kd,N1,N2)
!======================================================================
!     calculate the value of overlap determinant for given orbitals
!
!     Calls:  DET
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: kd,N1(kd),N2(kd)
      Integer :: i,j
      Real(8) :: ADET(kd*kd)
      Real(8), external :: DET, OBS

      if(kd.eq.0) then                
       VDET = 1.d0
      elseif(kd.eq.1) then
       VDET = OBS(N1(1),N2(1))
      elseif(kd.eq.2) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))
      elseif(kd.eq.3) then
       VDET = OBS(N1(1),N2(1))*OBS(N1(2),N2(2))*OBS(N1(3),N2(3)) +  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(3))*OBS(N1(3),N2(1)) +  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(1))*OBS(N1(3),N2(2)) -  &
              OBS(N1(1),N2(3))*OBS(N1(2),N2(2))*OBS(N1(3),N2(1)) -  &
              OBS(N1(1),N2(2))*OBS(N1(2),N2(1))*OBS(N1(3),N2(3)) -  &
              OBS(N1(1),N2(1))*OBS(N1(2),N2(3))*OBS(N1(3),N2(2)) 
      else                
       Do i=1,kd;  Do j=1,kd
         adet((i-1)*kd+j)=OBS(N1(i),N2(j))
       End do; End do
       VDET = DET(kd,adet)      
      end if

      End Function VDET

!======================================================================
      Integer Function IBORT(io,jo)
!======================================================================
!     recover old definition IBORT
!----------------------------------------------------------------------
      Use orb_overlaps

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), external :: OBS

      if(chan(io).eq.0.and.chan(jo).eq.0) &
        Call Stop_mpi (0,0,'IBORT for bound orbitals?')

      IBORT = NINT(OBS(io,jo))

      End Function IBORT



!======================================================================
      Subroutine Check_orb_overlaps
!======================================================================
!     the imposed orthogonal conditions given as < n1k1| n2k2>=0  
!----------------------------------------------------------------------
      Use orb_overlaps
      Use internal_file            !  ???

      Implicit none
      Integer :: i,j,l, i1,n1,l1,k1, i2,n2,l2,k2, is,js,ip, ii  
      Integer, external :: Ifind_bsorb
      Character(13) :: Aort

      Do ii = 1,nlines
      Aort = aline(ii)
      if(Aort(1:1).ne.'<') Cycle 
      Call EL4_nlk(Aort(2: 5),n1,l1,k1)
      i1 = Ifind_bsorb(n1,l1,k1)
      if(i1.eq.0) Cycle
      Call EL4_nlk(Aort(7:10),n2,l2,k2)
      i2 = Ifind_bsorb(n2,l2,k2)
      if(i2.eq.0) Cycle

      if(l1.ne.l2) Cycle; l = l1 
      read(Aort(13:13),'(i1)') j
      if(j.ne.0) Cycle

      i = jpl(i1);  j = jpl(i2)
      ip = ip_l(l)
      is = max(i,j)-ip; js = min(i,j)-ip
      ip = ip_ovl(l) + is*(is-1)/2 + js

      Cobs(ip) = 0.d0

      End do

      End Subroutine Check_orb_overlaps 

