!======================================================================
      Subroutine pre_det_exp 
!======================================================================
!     define the det. expansions and write the information
!     in scratch file 'nud'. 'nua' is just used for temporally storing.
!     Expansion is given for max. ML and MS,
!     and for all pairs of configurations
!-----------------------------------------------------------------------

      USE bsr_breit,    only: nud,nua,ioper,mktkdt
      USE term_exp,     only: ic_case, MLT,MST
      USE spin_orbitals,only: in,md,Lsym,Msym

      USE conf_LS,      only: msh,ne,LTOTAL,STOTAL,no,nn,ln,iq,kn,LS,LS1,LSI
      USE symc_list_LS, only: nsymc,IC_need,JC_need,IC_term1,IC_term2, &
                              LT_conf, ST_conf
      USE symt_list_LS, only: IT_sort, IT_need

      Implicit none 

      Integer :: i,j,ij,k,kk,kt,ktm,kdt,kdtn,kti,it,it1,it2,ic,jc, MLTi,MSTi
      Integer, Allocatable :: IP_kt(:)
      Integer, Allocatable :: IM_det(:,:),IS_det(:,:),JM_det(:,:),JS_det(:,:)
      Real(8), Allocatable :: CC_det(:,:)
      Integer :: IPT_jc(nsymc)
      Integer, External :: Iglq, DEF_ij, Iterm_LS

!----------------------------------------------------------------------
! ... loop over conf. symmeteries:

      ic_case = 0;  rewind(nud); rewind(nua)

      Do ic=1,nsymc;  if(IC_need(ic).eq.0) Cycle

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
	        LS(i,1)=Iterm_LS(ln(i),iq(i),0,LS(i,1),LS(i,2),LS(i,3))  ! ???
        End do
	        LSI(:,:,kt)=LS
       End do

       if(kt.eq.0) then; IC_need(ic)=0; Cycle; end if   ! ???

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
       ij=DEF_ij(ic,jc); if(JC_need(ij).eq.0) Cycle      
 
! ... define the det. expansion:

       MLT = min(Ltotal,LT_conf(jc))     !  ???
	      MST = min(Stotal,ST_conf(jc))

       rewind(nua)
	       Call Det_expn (nua,kt,kdt,MLT,MST)
       rewind(nua)      

       if(kdt.eq.0) Stop 'Pre_detexp: kdt = 0'

       if(1.d0*kt*kt*kdt.gt.1d9) then
        write(*,'(a,i10)') 'kt  =',kt
        write(*,'(a,i10)') 'kdt =',kdt
!        Stop 'too extensive dimensions'
       end if

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
       

       ktm = sqrt(mktkdt/kdt+1.d0); if(ktm.gt.kt) ktm=kt
       if(ktm.eq.kt) then
        ic_case=ic_case+1
        write(nud) ic,kt,kdt,Ltotal,Stotal,MLT,MST
        write(nud) IP_kt(1:kt)
        write(nud) CC_det
        write(nud) IM_det
        write(nud) IS_det
        write(nud) Msym(1:ne)
        write(nud) Lsym(1:ne)
!        write(*,*) 'i_case, kt,kdt',ic_case,kt,kdt
       else
        Do k = 1,kt,ktm 
         kk = k+ktm-1; if(kk.gt.kt) kk=kt
         JM_det = IM_det; JS_det = IS_det

         ! compress:

         j = 0
         Do i=1,kdt
          if(SUM(abs(CC_det(k:kk,i))).eq.0.d0) Cycle 
          j = j + 1; if(i.eq.j) Cycle 
          CC_det(k:kk,j) = CC_det(k:kk,i)
          JM_det(:,j)=IM_det(:,i)
          JS_det(:,j)=IS_det(:,i)
         End do
         kdtn = j

         ic_case=ic_case+1
         write(nud) ic,kk-k+1,kdtn,Ltotal,Stotal,MLT,MST
         write(nud) IP_kt(k:kk)
         write(nud) ((CC_det(i,j),i=k,kk),j=1,kdtn)
         write(nud) JM_det(:,1:kdtn)
         write(nud) JS_det(:,1:kdtn)
         write(nud) Msym(1:ne)
         write(nud) Lsym(1:ne)
         write(*,*) 'n_case, kt,kdt',ic_case,kk-k+1,kdtn
        End do

       end if

       Deallocate(CC_det,IM_det,IS_det,JM_det,JS_det)

! ... mark the configurations with the same term:   ! ???

       Do i = jc,nsymc
        MLTi = min(Ltotal,LT_conf(i))
	       MSTi = min(Stotal,ST_conf(i))
        if(MLTi.eq.MLT.and.MSTi.eq.MST) IPT_jc(i)=0
       End do

      End do    ! over jc
      End do    ! over ic

      if(Allocated(IP_kt)) Deallocate(IP_kt)

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

