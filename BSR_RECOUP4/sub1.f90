!======================================================================
      Subroutine SUB1
!======================================================================
!     drive routine for one partial wave
!----------------------------------------------------------------------
      Use bsr_recoup

      Use channels_ion, only: nlsp_LS => nlsp, nch_LS => nch, &
                              ipar_LS => ipar, lpar_LS => lpar, ispar_LS => ispar
      
      Implicit real(8) (A-H,O-Z)

      Call R_channel(nut,klsp)
      i=len_trim(AF_mat)
      write(AF_mat(i-2:i),'(i3.3)') klsp

      open(nui,file=AF_mat,form='UNFORMATTED',action='WRITE')
      write(nui) ns,nch,npert

! ... overlaps

      Call Allocate_matrix(m)

      mui = kui
      Do ilsp = 1,nlsp_LS; kch = nch_LS(ilsp); LL = 2*lpar_LS(ilsp)+1 

       if(ipar.ne.ipar_LS(ilsp)) Cycle
       if(ITRI (jpar,LL,ispar_LS(ilsp)).eq.0) Cycle

       i=len_trim(BF_mat)
       write(BF_mat(i-2:i),'(i3.3)') ilsp
       Call Check_file(BF_mat)
       mui = mui + 1
       open(mui,file=BF_mat,form='UNFORMATTED',action='READ')
       rewind(mui)       
       read(mui) ns1, nch1, npert1

       Do 
        read(mui) ic,jc
        if(ic.le.0) Exit
        if(ic.gt.kch.or.jc.gt.kch) Stop 'not channel-channel block'
        read(mui) x(:,:)
        Call Add_block(ilsp,ic,jc,1)
         if(ic.ne.jc)  Call Add_block(ilsp,jc,ic,1) 
       End do

      End do

      Call Record_matrix(nui,1)

! ... interaction matrix

      Call Allocate_matrix(m)

      mui = kui
      Do ilsp = 1,nlsp_LS; kch = nch_LS(ilsp); LL = 2*lpar_LS(ilsp)+1 

       if(ipar.ne.ipar_LS(ilsp)) Cycle
       if(ITRI (jpar,LL,ispar_LS(ilsp)).eq.0) Cycle

       mui = mui + 1

       Do 
        read(mui) ic,jc
        if(ic.le.0) Exit
        if(ic.gt.kch.or.jc.gt.kch) Stop 'not channel-channel block'
        read(mui) x(:,:)
        Call Add_block(ilsp,ic,jc,1)
         if(ic.ne.jc)  Call Add_block(ilsp,jc,ic,1) 
       End do

      End do

      Call Record_matrix(nui,0)

!------------------------------------------------------------------------
! ... asymptotic coefficients

      ms = ns; ns=0
      Call Allocate_matrix(m)
      ns = ms
      if(allocated(acf)) deallocate(acf)
      Allocate(acf(nch,nch,0:mk));  acf = 0.d0

      mui = kui
      Do ilsp = 1,nlsp_LS; kch = nch_LS(ilsp); LL = 2*lpar_LS(ilsp)+1 

       if(ipar.ne.ipar_LS(ilsp)) Cycle
       if(ITRI (jpar,LL,ispar_LS(ilsp)).eq.0) Cycle
       mui = mui + 1

       if(allocated(bcf)) deallocate(bcf)      
       read(mui) km
       Allocate(bcf(kch,kch,0:km))                       
       read(mui) bcf

       Call print_acf (pri,kch,km,bcf,eps_acf)

       Do ic=1,kch
        Do jc=1,ic
         Call Add_block(ilsp,ic,jc,2)
         if(ic.ne.jc)  Call Add_block(ilsp,jc,ic,2) 
        End do
       End do

      End do

      Do ic=1,nch
       Do jc=1,ic
         acf(jc,ic,:) = acf(ic,jc,:)
       End do
      End do

      write(nui) mk
      write(nui) acf

      if(pri_ac.gt.0) Call print_acf (pri,nch,mk,acf,eps_acf)

      End Subroutine SUB1

!======================================================================
      Subroutine print_acf (pri,nch,mk,acf,eps_acf)
!======================================================================
      Implicit none
      Integer, intent(in) :: pri,nch,mk
      Real(8), intent(in) :: acf(nch,nch,0:mk),eps_acf
      Integer :: k,i,j,ij,i1,i2
      Character(200) :: line
      line = ' '

      write(pri,'(/a,a,i5,a,i2/)') 'Asymptotic coefficients: i,j, ACF(i,j,k)', &
        '   nch =',nch,'   mk =',mk
      Do k=0,mk
       if(SUM(acf(:,:,k)).eq.0) Cycle
       write(pri,'(a,i2)') 'k = ',k
       ij = 0
       Do i=1,nch; Do j = 1,i      
        if(abs(acf(i,j,k)).lt.eps_acf) Cycle
        i1=ij*20+1; i2=i1+19 
        write(line(i1:i2),'(2i4,E12.3)') i,j,acf(i,j,k)
        ij=ij+1
        if(ij.lt.5) Cycle
        write(pri,'(a)') line; ij=0
       End do; End do
       if(ij.eq.0) Cycle
       i1=1; i2=ij*20
       write(pri,'(a)') line(i1:i2)
      End do

      End Subroutine print_acf

