!======================================================================
      Subroutine pre_det_exp(nud,mkt,mkdt,pri) 
!======================================================================
!     define the det. expansions and write the information
!     in file 'nud'. ' Expansions are given for max. ML and MS,
!     and for all pairs of configurations
!-----------------------------------------------------------------------
      Use term_exp,     only: ic_case, MLT,MST
      Use spin_orbitals,only: in,md,Lsym,Msym
      Use conf_LS,      only: msh,ne,LTOTAL,STOTAL,no,nn,ln,iq,kn,LS,LS1,LSI
      Use symc_list_LS, only: nsymc,IC_need,JC_need,IC_term1,IC_term2, &
                              LT_conf, ST_conf
      Use symt_list_LS, only: IT_sort, IT_need

      Implicit none 

      Integer, intent(in) :: nud,mkt,mkdt,pri
      Integer :: i,j,k,kk,kt,ktn,ktm,kdt,kdtn,kti,it,it1,it2,ic,jc, MLTi,MSTi, m,mm
      Integer :: nua,nub
      Integer, allocatable :: IP_kt(:)
      Integer, allocatable :: IM_det(:,:),IS_det(:,:),JM_det(:,:),JS_det(:,:)
      Real(8), allocatable :: CC_det(:,:)
      Integer :: IPT_jc(nsymc)
      Integer, external :: Iglq, Iterm_LS
      Integer(8) :: ij
      Integer(8), external :: DEF_ij8

      Integer :: max_kt, max_kdt 
      Real(8) :: S, SS, SM
      Character(80) :: conf

      Call Find_free_unit(nua)
      Open(nua,form='UNFORMATTED',status='SCRATCH')
      Call Find_free_unit(nub)
      Open(nub,form='UNFORMATTED',status='SCRATCH')
      
      SM = 0
      max_kt = 0
      max_kdt = 0

      if(pri.gt.0) then
       write(pri,'(/a,i10 )') 'mkt  =',mkt 
       write(pri,'( a,i10/)') 'mkdt =',mkdt 
      end if

!----------------------------------------------------------------------
! ... loop over conf. symmeteries:

      ic_case = 0;  rewind(nud); rewind(nua)

      Do ic=1,nsymc;  if(IC_need(ic).eq.0) Cycle

      Call Symc_conf(ic,conf)

! ... define configuration ic:

       Call Get_symc_LS(ic,LTOTAL,STOTAL,no,nn,ln,iq,kn)
       k=1
       Do i=1,no; in(i)=k; k=k+iq(i); md(i)=Iglq(ln(i),iq(i)); End do

! ... number of terms:

       it1=IC_term1(ic); it2=IC_term2(ic); kti=it2-it1+1

! ... define needed values for shell and intermidiate terms (LSI)
! ... and the number of different angular symmetries (kt):

       if(Allocated(IP_kt)) Deallocate(IP_kt);  Allocate(IP_kt(kti))
       if(Allocated(LSI)) Deallocate(LSI);  Allocate(LSI(msh,5,kti))

       kt=0
       Do k =it1,it2; it=IT_sort(k); if(IT_need(it).eq.0) Cycle
        kt = kt + 1; IP_kt(kt) = it
        Call Get_symt_LS(it,ic,no,LS)
        Do i=1,no
         LS(i,1)=Iterm_LS(ln(i),iq(i),0,LS(i,1),LS(i,2),LS(i,3))  
        End do
        LSI(:,:,kt)=LS
       End do

       if(kt.eq.0) then; IC_need(ic)=0; Cycle; end if   

       Do i=1,no
       Do k=1,5
        LS1(i,k)=maxval(LSI(i,k,1:kt))
       End do
       End do

!----------------------------------------------------------------------       
! ... cycle over all config. symmetries:
   
      IPT_jc = 1   
      Do jc = 1,nsymc; if(IC_need(jc).eq.0) Cycle 

       if(IPT_jc(jc).eq.0) Cycle    
       ij=DEF_ij8(ic,jc); if(JC_need(ij).eq.0) Cycle      
 
! ... define the det. expansion:

       MLT = min(Ltotal,LT_conf(jc))     
       MST = min(Stotal,ST_conf(jc))

       rewind(nua)
       Call Det_expn (nua,kt,kdt,MLT,MST)
       rewind(nua)      

       if(kdt.eq.0) Stop 'Pre_detexp: kdt = 0'

! ... record results:

       k=1
       Do i=1,no
        Do j=k,k+iq(i)-1; Msym(j)=i; Lsym(j)=ln(i); End do
        k=k+iq(i)
       End do

       Allocate(CC_det(kt,kdt),IM_det(ne,kdt),IS_det(ne,kdt),&
                               JM_det(ne,kdt),JS_det(ne,kdt))
       Do i=1,kdt
        read(nua) CC_det(:,i),IM_det(:,i),IS_det(:,i)
       End do

       Do k = 1,kt,mkt 
         kk = k+mkt-1; if(kk.gt.kt) kk=kt; ktn=kk-k+1
        Do m = 1,kdt,mkdt 
          mm = m+mkdt-1; if(mm.gt.kdt) mm=kdt; kdtn=mm-m+1

        ic_case=ic_case+1
        write(nud) ic,ktn,kdtn,Ltotal,Stotal,MLT,MST
        write(nud) IP_kt(k:kk)
        write(nud) CC_det(k:kk,m:mm)
        write(nud) IM_det(1:ne,m:mm)
        write(nud) IS_det(1:ne,m:mm)
        write(nud) Msym(1:ne)
        write(nud) Lsym(1:ne)

        write(nub) ic 

        if(pri.gt.0) then
         S = 4.0*ktn + 8.0*ktn*kdtn + 8.0*ne*kdtn + 8.0*ne + 7*4 
         S = S / (1024*1024)
         it = count(CC_det(k:kk,m:mm).ne.0.d0)
         SS = 1.d0*it/(ktn*kdtn) 
         write(pri,'(a,2i5,i6,i10,f10.3,f10.2,a,5x,a,f10.3)')  &
           'ic_case,ic,kt,kdt',ic_case,ic,ktn,kdtn,SS,S,' Mb',trim(conf)
         if(S.gt.SM) SM=S
         if(kt.gt.max_kt) max_kt=kt
         if(kdt.gt.max_kdt) max_kdt=kdt
        end if      

        End do
       End do

       Deallocate(CC_det,IM_det,IS_det,JM_det,JS_det)

! ... mark the configurations with the same term:   

       Do i = jc,nsymc
        MLTi = min(Ltotal,LT_conf(i))
        MSTi = min(Stotal,ST_conf(i))
        if(MLTi.eq.MLT.and.MSTi.eq.MST) IPT_jc(i)=0
       End do

      End do    ! over jc
      End do    ! over ic

      if(Allocated(IP_kt)) Deallocate(IP_kt)

      Allocate(IP_kt(ic_case))
      rewind(nub)
      Do i=1,ic_case; read(nub) IP_kt(i); End do
      write(nud) IP_kt
      Deallocate(IP_kt)

      write(nud) ic_case

      if(pri.gt.0) then
       write(pri,'(/a,i10)') 'max_kt  =',max_kt 
       write(pri,'( a,i10)') 'max_kdt =',max_kdt 
       write(pri,'( a,F10.1)') 'max_mem =',SM 
      end if

      End Subroutine pre_det_exp 


!======================================================================
      Subroutine Det_expn (nua,kt,kdt,MLT,MST) 
!======================================================================
!     determined all possible determinants and their coefficients
!     for given set of terms (kt). Results are recorded to unit 'nua'
!----------------------------------------------------------------------

      USE conf_LS,       only: ne,no,ln,iq,LS,LS1,LSI
      USE spin_orbitals, only: in,md,Msym,Ssym,MS_orb,ML_orb 

      Implicit none

      Integer :: nua,kt,kdt,MLT,MST, kd,i,j,k,m,ii
      Integer :: nd(ne),idet(ne),ML(ne),MS(ne),MLp(ne),MSp(ne)  
      Real(8) :: Cdet(kt)
      Real(8) :: C
      Real(8), External :: Clebsh, DETC_sh

      kd=0; kdt=0; i=1; nd(i)=1              
    1 kd = kd + 1
      ii = in(i)
      Call DET_sh(ln(i),iq(i),nd(i),ML(i),MS(i),Idet(ii)) 

      m = iabs(ML(i)); if(ML(i).lt.0) m=m+2
      if(m.gt.LS1(i,2)) go to 2
      m = iabs(MS(i)); if(MS(i).lt.0) m=m+2
      if(m.gt.LS1(i,3)) go to 2

      if(i.eq.1) then
       MLp(1)=ML(1); MSp(1)=MS(1)
      else
       MLp(i) = MLp(i-1)+ML(i)-1
       MSp(i) = MSp(i-1)+MS(i)-1

       m = iabs(MLp(i)); if(MLp(i).lt.0) m = m + 2
       if(m.gt.LS1(i,4)) go to 2
       m = iabs(MSp(i)); if(MSp(i).lt.0) m = m + 2
       if(m.gt.LS1(i,5)) go to 2
      end if

      if(i.lt.no) then; i = i + 1; nd(i) = 1;  go to 1; end if

      if(MLp(no).ne.MLt) go to 2
      if(MSp(no).ne.MSt) go to 2

!--------------------------------------------------------------------
!                                            coefficient calculation:
      Cdet = 0.d0
      Do k=1,kt; C=1.d0; LS(:,:)=LSI(:,:,k)
       Do j=1,no
        C=C*DETC_sh(ln(j),iq(j),LS(j,1),nd(j))
        if(C.eq.0.d0) Exit
       End do
       if(C.eq.0.d0) Cycle
       Do j=2,no
        C=C*Clebsh(LS(j-1,4),MLp(j-1), &
                   LS(j  ,2),ML (j  ), &
                   LS(j  ,4),MLp(j  ))
        if(C.eq.0.d0) Exit
        C=C*Clebsh(LS(j-1,5),MSp(j-1), &
                   LS(j  ,3),MS (j  ), &
                   LS(j  ,5),MSp(j  ))
        if(C.eq.0.d0) Exit
       End do
       Cdet(k) = C
      End do

      kdt=kdt+1
      Do j = 1,ne
       Msym(j)=ML_orb(idet(j))
       Ssym(j)=MS_orb(idet(j))
      End do
      write(nua) Cdet,Msym(1:ne),Ssym(1:ne)

!--------------------------------------------------------------------
    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue

      End Subroutine Det_expn

