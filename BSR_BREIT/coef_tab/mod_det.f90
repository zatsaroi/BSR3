!----------------------------------------------------------------------
!  Here are placed the modules for processing the overlap determinant
!  factors. In principal, all four modules DET_lis, DEF_list,
!  DET_lis, DEF_list have the identical structure and differ only
!  by names for some variables. 
!  The connection between DET, DEF and NDET,NDEF lists are given
!  in subroutine NDET_IDET.
!  We introduced the additional NDET, NDEF lists in hope to reduce
!  the seeking time in large common list DET, DETF. 
!----------------------------------------------------------------------




!======================================================================
      MODULE DET_list
!======================================================================
!
!     Containes the overlap determinants for all config. symmetries.
!
!----------------------------------------------------------------------

      Use param_br, ONLY: isd,jsd

      Implicit none
      Save

      Integer(4) :: ndet = 0       ! number of determinants
      Integer(4) :: mdet = 0       ! current dimentsion of list
      Integer(4) :: idet = isd     ! supposed max. dimentsion
      Integer(4) :: jdet = jsd     ! average size of one det. 
      Integer(4) :: kdet = 0       ! dimension of all det.s 
      
      Integer(4), Allocatable, Dimension(:) :: KPD,IPD,NPD   

      End MODULE DET_list


!======================================================================
      Subroutine alloc_det(m)
!======================================================================

      Use param_br, ONLY: nus
      Use DET_list

      Implicit none
      Integer(4), Intent(in) :: m
      Integer(4) :: k

      if(m.le.0) then
       if(allocated(KPD)) Deallocate (KPD,IPD,NPD)
       mdet = 0; ndet = 0; kdet =0
      elseif(.not.allocated(KPD)) then
       mdet = m; kdet = mdet*jdet
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet))
      elseif(m.le.mdet) then
       Return
      elseif(ndet.eq.0) then
       Deallocate (KPD,IPD,NPD)
       mdet = m; kdet = mdet*jdet
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) KPD(1:ndet)
       write(nus) IPD(1:ndet)
       write(nus) NPD(1:kdet)
       Deallocate (KPD,IPD,NPD)
       k=kdet; mdet = m; kdet = mdet*jdet
       Allocate(KPD(1:mdet),IPD(1:mdet),NPD(1:kdet))
       rewind(nus)
       read(nus) KPD(1:ndet)
       read(nus) IPD(1:ndet)
       read(nus) NPD(1:k)
       Close(nus)
       write(*,*) ' Realloc_det: new dimension = ', mdet,jdet
      end if

      End Subroutine alloc_DET


!----------------------------------------------------------------------
      Integer(4) Function Iadd_det (kd,NP)
!----------------------------------------------------------------------
!
!     add the overlap determinant to DET_list
!
!----------------------------------------------------------------------

      Use det_list

      Implicit none 
      Integer(4) , Intent(in) :: kd
      Integer(4) , Intent(in), Dimension(kd) :: NP
      Integer(4) :: i,j,ip

      Iadd_det = 0
      if(kd.le.0) Return
      if(mdet.eq.0) Call Alloc_DET(idet)

! ... check: is the same det. in the list:

      Do i=1,ndet
       if(KPD(i).ne.kd) Cycle
       ip=IPD(i); Iadd_det = i
       Do j=1,kd
        if(NP(j).eq.NPD(ip+j)) Cycle; Iadd_det = 0; Exit
       End do
       if(Iadd_det.ne.0) Return
      End do

! ... Add new det.:

      ndet=ndet+1; ip=0; if(ndet.gt.1) ip=IPD(ndet-1)+KPD(ndet-1)
      if(ndet.eq.mdet.or.ip+kd.gt.kdet) then
       jdet = kdet/ndet + 1; Call Alloc_det(mdet+idet)
      end if
      KPD(ndet)=kd; IPD(ndet)=ip; NPD(ip+1:ip+kd)=NP(1:kd)
      Iadd_det=ndet

      End Function Iadd_det



!======================================================================
      MODULE def_list
!======================================================================
!
!     Containes the overlap factors, as the list of the number of
!     involved overlap determinants and their positions in the 
!     common det_list.
!
!     KPF(i) - number of det.s in the overlap factor 'i'  (kd)
!     IPF(i) - pointer on the list of corr. det.s in the NPF  (ip)
!     NPF(ip+1:ip+kd) - list of pointers on det.s in the det_list
!                       and their powers
!----------------------------------------------------------------------

      Use param_br, ONLY: isf,jsf

      Implicit none
      Save

      Integer(4) :: ndef = 0       ! number of determinants
      Integer(4) :: mdef = 0       ! current dimentsion of list
      Integer(4) :: idef = isf     ! supposed max. dimentsion  
      Integer(4) :: jdef = jsf     ! average number of det.s 
      Integer(4) :: kdef = 0       ! dimension of all def.s 
      
      Integer(4), Allocatable, Dimension(:) :: KPF,IPF,NPF   

      End MODULE def_list


!======================================================================
      Subroutine alloc_def(m)
!======================================================================

      USE def_list
      Use param_br, ONLY: nus

      Implicit none
      Integer(4), Intent(in) :: m
      Integer(4) :: k

      if(m.le.0) then
       if(allocated(KPF)) Deallocate(KPF,IPF,NPF)
       mdef = 0; ndef = 0; kdef = 0
      elseif(.not.allocated(KPF)) then
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef))
      elseif(m.le.mdef) then
       Return
      elseif(ndef.eq.0) then
       Deallocate (KPF,IPF,NPF)
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) KPF(1:ndef)
       write(nus) IPF(1:ndef)
       write(nus) NPF(1:kdef)
       Deallocate (KPF,IPF,NPF)
       k=kdef; mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef))
       rewind(nus)
       read(nus) KPF(1:ndef)
       read(nus) IPF(1:ndef)
       read(nus) NPF(1:k)
       Close(nus)
       write(*,*) ' Realloc_def: new dimension = ', mdef,jdef
      end if

      End Subroutine alloc_def



!----------------------------------------------------------------------
      Integer(4) Function Iadd_def (kd,NP)
!----------------------------------------------------------------------
!
!     add the overlap factor to def_list
!
!     kd    - number of det.s
!     NP(i) - pointer for the i-th det.
!
!----------------------------------------------------------------------

      USE def_list

      Implicit none
      Integer(4) , Intent(in) :: kd
      Integer(4) , Dimension(kd) :: NP
      Integer(4) :: i,j,ip

      if(kd.le.0) Return
      if(mdef.eq.0) Call Alloc_def(idef) 

! ... check: is the same def. in the list:

      Do i=1,ndef
       if(KPF(i).ne.kd) Cycle
       ip=IPF(i); Iadd_def = i
       Do j=1,kd
        if(NP(j).eq.NPF(ip+j)) Cycle; Iadd_def = 0; Exit
       End do
       if(Iadd_def.ne.0) Return
      End do

! ... Add new def.:

      ndef=ndef+1; ip=0; if(ndef.gt.1) ip=IPF(ndef-1)+KPF(ndef-1)
      if(ndef.eq.mdef.or.ip+kd.gt.kdef) then
       jdef = (ip+kd)/ndef + 1; Call Alloc_def(mdef+idef)
      end if
      KPF(ndef)=kd; IPF(ndef)=ip; NPF(ip+1:ip+kd)=NP(1:kd)
      Iadd_def=ndef

      End Function Iadd_def


!======================================================================
      MODULE NDET_list
!======================================================================
!
!     Containes the overlap determinants for all config. symmetries.
!
!----------------------------------------------------------------------

      Use param_br, ONLY: isd,jsd

      Implicit none
      Save

      Integer(4) :: ndet = 0       ! number of determinants
      Integer(4) :: mdet = 0       ! current dimentsion of list
      Integer(4) :: idet = isd     ! supposed max. dimentsion
      Integer(4) :: jdet = jsd     ! average size of one det. 
      Integer(4) :: kdet = 0       ! dimension of all det.s 
      
      Integer(4), Allocatable, Dimension(:) :: KPD,IPD,NPD   

      End MODULE NDET_list


!======================================================================
      Subroutine alloc_ndet(m)
!======================================================================

      Use param_br, ONLY: nus
      Use NDET_list

      Implicit none
      Integer(4), Intent(in) :: m
      Integer(4) :: k

      if(m.le.0) then
       if(allocated(KPD)) Deallocate (KPD,IPD,NPD)
       mdet = 0; ndet = 0; kdet =0
      elseif(.not.allocated(KPD)) then
       mdet = m; kdet = mdet*jdet
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet))
      elseif(m.le.mdet) then
       Return
      elseif(ndet.eq.0) then
       Deallocate (KPD,IPD,NPD)
       mdet = m; kdet = mdet*jdet
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) KPD(1:ndet)
       write(nus) IPD(1:ndet)
       write(nus) NPD(1:kdet)
       Deallocate (KPD,IPD,NPD)
       k=kdet; mdet = m; kdet = mdet*jdet
       Allocate(KPD(1:mdet),IPD(1:mdet),NPD(1:kdet))
       rewind(nus)
       read(nus) KPD(1:ndet)
       read(nus) IPD(1:ndet)
       read(nus) NPD(1:k)
       Close(nus)
       write(*,*) ' Alloc_ndet: new dimension = ', mdet,kdet
      end if

      End Subroutine alloc_NDET


!----------------------------------------------------------------------
      Integer(4) Function Nadd_det (kd,NP)
!----------------------------------------------------------------------
!
!     add the overlap determinant to DET_list
!
!----------------------------------------------------------------------

      Use ndet_list

      Implicit none 
      Integer(4) , Intent(in) :: kd
      Integer(4) , Intent(in), Dimension(kd) :: NP
      Integer(4) :: i,j,ip

      Nadd_det = 0
      if(kd.le.0) Return
      if(mdet.eq.0) Call Alloc_NDET(idet)

! ... check: is the same det. in the list:

      Do i=1,ndet
       if(KPD(i).ne.kd) Cycle
       ip=IPD(i); Nadd_det = i
       Do j=1,kd
        if(NP(j).eq.NPD(ip+j)) Cycle; Nadd_det = 0; Exit
       End do
       if(Nadd_det.ne.0) Return
      End do

! ... Add new det.:

      ndet=ndet+1; ip=0; if(ndet.gt.1) ip=IPD(ndet-1)+KPD(ndet-1)
      if(ndet.eq.mdet.or.ip+kd.gt.kdet) then
       jdet = (ip+kd)/ndet + 1; Call Alloc_ndet(mdet+idet)
      end if
      KPD(ndet)=kd; IPD(ndet)=ip; NPD(ip+1:ip+kd)=NP(1:kd)
      Nadd_det=ndet

      End Function Nadd_det



!======================================================================
      MODULE NDEF_list
!======================================================================
!
!     Containes the overlap factors, as the list of the number of
!     involved overlap determinants and their positions in the 
!     common det_list.
!
!     KPF(i) - number of det.s in the overlap factor 'i'  (kd)
!     IPF(i) - pointer on the list of corr. det.s in the NPF  (ip)
!     NPF(ip+1:ip+kd) - list of pointers on det.s in the det_list
!                       and their powers
!----------------------------------------------------------------------

      Use param_br, ONLY: isf,jsf

      Implicit none
      Save

      Integer(4) :: ndef = 0       ! number of determinants
      Integer(4) :: mdef = 0       ! current dimentsion of list
      Integer(4) :: idef = isf     ! supposed max. dimentsion  
      Integer(4) :: jdef = jsf     ! average number of det.s 
      Integer(4) :: kdef = 0       ! dimension of all def.s 
      
      Integer(4), Allocatable, Dimension(:) :: KPF,IPF,NPF   

      End MODULE NDEF_list


!======================================================================
      Subroutine alloc_ndef(m)
!======================================================================

      USE ndef_list
      Use param_br, ONLY: nus

      Implicit none
      Integer(4), Intent(in) :: m
      Integer(4) :: k

      if(m.le.0) then
       if(allocated(KPF)) Deallocate(KPF,IPF,NPF)
       mdef = 0; ndef = 0; kdef = 0
      elseif(.not.allocated(KPF)) then
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef))
      elseif(m.le.mdef) then
       Return
      elseif(ndef.eq.0) then
       Deallocate (KPF,IPF,NPF)
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) KPF(1:ndef)
       write(nus) IPF(1:ndef)
       write(nus) NPF(1:kdef)
       Deallocate (KPF,IPF,NPF)
       k=kdef; mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef))
       rewind(nus)
       read(nus) KPF(1:ndef)
       read(nus) IPF(1:ndef)
       read(nus) NPF(1:k)
       Close(nus)
       write(*,*) ' Realloc_ndef: new dimension = ', mdef,jdef
      end if

      End Subroutine alloc_ndef



!----------------------------------------------------------------------
      Integer(4) Function Nadd_def (kd,NP)
!----------------------------------------------------------------------
!
!     add the overlap factor to def_list
!
!     kd    - number of det.s
!     NP(i) - pointer for the i-th det.
!
!----------------------------------------------------------------------

      USE ndef_list

      Implicit none
      Integer(4) , Intent(in) :: kd
      Integer(4) , Dimension(kd) :: NP
      Integer(4) :: i,j,ip

      if(kd.le.0) Return
      if(mdef.eq.0) Call Alloc_ndef(idef) 

! ... check: is the same def. in the list:

      Do i=1,ndef
       if(KPF(i).ne.kd) Cycle
       ip=IPF(i); Nadd_def = i
       Do j=1,kd
        if(NP(j).eq.NPF(ip+j)) Cycle; Nadd_def = 0; Exit
       End do
       if(Nadd_def.ne.0) Return
      End do

! ... Add new def.:

      ndef=ndef+1; ip=0; if(ndef.gt.1) ip=IPF(ndef-1)+KPF(ndef-1)
      if(ndef.eq.mdef.or.ip+kd.gt.kdef) then
       jdef = kdef/ndef + 1; Call Alloc_ndef(mdef+idef)
      end if
      KPF(ndef)=kd; IPF(ndef)=ip; NPF(ip+1:ip+kd)=NP(1:kd)
      Nadd_def=ndef

      End Function Nadd_def



!======================================================================
     Subroutine Ndet_Idet
!======================================================================
!
!    Convert NDET and NDEF lists to the common DET and DEF lists. 
!    Connection is given in IPF array.
!    The NDET and NDEF lists are nulefied.
!----------------------------------------------------------------------

      Use param_br, ONLY: me, ibf
      Use NDET_list
      Use NDEF_list 

      IMPLICIT NONE
      Integer(4), Dimension(me) :: NP
      Integer(4) :: i,j,ip,jp,id,kd,ns
      Integer(4), External :: Iadd_det, Iadd_def, ISORT


      if(ndet.le.0) Stop ' NDET_IDET: ndet <= 0'
      Do id=1,ndet
       kd=KPD(id); ip=IPD(id); NP(1:kd)=NPD(ip+1:ip+kd)
       IPD(id) = Iadd_det(kd,NP)
      End do
      ndet = 0

      if(ndef.le.0) Stop ' NDET_IDET: ndef <= 0'
      Do id=1,ndef
       kd=KPF(id); ip=IPF(id)
       Do i=1,kd
        j=NPF(ip+i)/ibf; jp=IPD(j); ns=mod(NPF(ip+i),ibf)
        NP(i) = jp*ibf + ns 
       End do
       i = ISORT (kd,NP)                             
       IPF(id) = Iadd_def (kd,NP)
      End do
      ndef = 0

      End Subroutine Ndet_Idet


!======================================================================
      Integer(4) Function ISORT (n,NN)
!======================================================================
!
!     simple sort of integer array NN(1:n)
!
!----------------------------------------------------------------------

      IMPLICIT NONE

      Integer(4) :: n, i,j,k
      Integer(4), Dimension(n) :: NN
      
      ISORT = 0
      Do i=1,n-1
       Do j=i+1,n
        if(nn(i).le.nn(j)) Cycle
        k=nn(i); nn(i)=nn(j); nn(j)=k; ISORT=ISORT+1
       End do
      End do

      End Function ISORT 

