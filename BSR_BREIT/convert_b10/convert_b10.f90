!======================================================================
!     Program    convert_b10
!======================================================================
!     convert int_bnk from b10 format (ZCONF90) to b20 format (ZCONFLS)
!----------------------------------------------------------------------

      USE zconst
      USE symc_list_LS;  USE symt_list_LS;  USE conf_LS
      USE det_list; USE def_list

      Implicit none

      Integer :: nu1=1; Character(80) :: AF_inp = 'int_bnk.90' 
      Integer :: nu2=2; Character(80) :: AF_out = 'int_bnk.LS' 

      Character(26), Allocatable :: ASYM(:)
      Character(40), Allocatable :: BSYM(:)

! ... JT_conf(it) gives the configuration for given term 'it'

      Integer, Allocatable :: JT_conf(:)

      Integer, External :: Iadd_symc_LS, Iadd_symt_LS 
      Integer :: i,j, nc,nt, it,ic,ii,im, ndim1,ndim2

!-----------------------------------------------------------------------

      Call Read_aarg('inp',AF_inp)
      Call Check_file(AF_inp)
      open(nu1,file=AF_inp,form='UNFORMATTED')

      Call Read_aarg('out',AF_out)
      open(nu2,file=AF_out,form='UNFORMATTED')

! ... read old information:

      read(nu1) nt, nc
      if(Allocated(JT_conf)) Deallocate (JT_conf,ASYM,BSYM)
      Allocate(JT_conf(nt),ASYM(nc),BSYM(nt))
      Call R_a (nu1,mrecl,nc,ASYM)    
      Call R_a (nu1,mrecl,nt,BSYM)    
      Call R_i4(nu1,mrecl,nt,JT_conf) 

! ... define conf. symmetries:

      Do it = 1,nt; ic=JT_conf(it) 

       Call Decode_conf(ASYM(ic),BSYM(it))

       iconf = Iadd_symc_LS(LS(no,4),LS(no,5),no,iq,ln)
       iterm = Iadd_symt_LS(iconf,no,LS)

       if(ic.ne.iconf) Stop 'inconsitent ic & iconf'
       if(it.ne.iterm) Stop 'inconsitent it & iterm'

      End do

      rewind(nu2)
      Call Write_symc_LS(nu2)
      Call Write_symt_LS(nu2)

      ic = 4*nsymc + 2*lsymc
      it = 2*nsymt + 5*lsymt

      write(*,'(/a/)') 'INT_BNK contains:'
      write(*,'(a,i9,a,f6.2,a)') 'number of conf.symmetries = ',NC, &
      '   =>  ', 4*ic /(1024d0*1024d0), '   Mb'
      write(*,'(a,i9,a,f6.2,a)') 'number of ang. symmetries = ',NT, &
      '   =>  ', 4*it /(1024d0*1024d0), '   Mb'

! ... rewrite oper:

      Allocate(IT_oper(1:noper,nsymt*(nsymt+1)/2))

      read(nu1) ii; ndim1=noper; ndim2=nsymt*(nsymt+1)/2
      Call R_i1(nu1,mrecl,ndim1,ndim2,ii,IT_oper)
      
      Call Write_oper_LS(nu2)

      im = ndim1*ndim2 
      write(*,'(a,T40,a,f6.2,a)') 'IT_oper',' =>  ', &
        im /(1024d0*1024d0), '   Mb'

! ... rewrite determinants:

      read(nu1) ndet,ldet; jdet = ldet/ndet + 1
      Call Alloc_det(ndet)
      Call R_i4(nu1,mrecl,ndet,KPD)  ! read(nub) KPD(1:ndet)
      Call R_i4(nu1,mrecl,ndet,IPD)  ! read(nub) IPD(1:ndet)
      Call R_i4(nu1,mrecl,ldet,NPD)  ! read(nub) NPD(1:ldet)

      read(nu1) ndef,ldef; jdef = ldef/ndef + 1
      Call Alloc_def(ndef)
      Call R_i4(nu1,mrecl,ndef,KPF)  ! read(nub) KPF(1:ndef)
      Call R_i4(nu1,mrecl,ndef,IPF)  ! read(nub) IPF(1:ndef)
      Call R_i4(nu1,mrecl,ldef,NPF)  ! read(nub) NPF(1:ldef)

      write(nu2) ndet,ldet
      Call W_i4(nu2,mrecl,ndet,KPD)  ! write(nur) KPD(1:ndet)
      Call W_i4(nu2,mrecl,ndet,IPD)  ! write(nur) IPD(1:ndet)
      Call W_i4(nu2,mrecl,ldet,NPD)  ! write(nur) NPD(1:ldet)

      write(nu2) ndef,ldef
      Call W_i4(nu2,mrecl,ndef,KPF)  ! write(nur) KPF(1:ndef)
      Call W_i4(nu2,mrecl,ndef,IPF)  ! write(nur) IPF(1:ndef)
      Call W_i4(nu2,mrecl,ldef,NPF)  ! write(nur) NPF(1:ldef)

      im = 2*ndet + ldet
      write(*,'(a,i9,a,f6.2,a)') 'number of det. overlaps   = ',ndet, &
      '   =>  ', 4*im /(1024d0*1024d0), '   Mb'
      im = 2*ndef + ldef
      write(*,'(a,i9,a,f6.2,a)') 'number of det. factors    = ',ndef, &
      '   =>  ', 4*it /(1024d0*1024d0), '   Mb'

! ... rewrite coefficients:

      Call RW(nu1,nu2,nc)

      im = 20 * nc 
      write(*,'(a,i9,a,f6.2,a)') 'number of coeffocients    = ',nc, &
      '   =>  ', im /(1024d0*1024d0), '   Mb'


      End ! program convert_b10


!======================================================================
      Subroutine Decode_conf(ASYM,BSYM)
!======================================================================
!     decode the config. from 'sym' format to ZOI format
!----------------------------------------------------------------------

      USE conf_LS

      Character(26), Intent(in) ::  ASYM  ! config. symmetry
      Character(40), Intent(in) ::  BSYM  ! angular symmetry
      Integer :: i,j,k
      Integer, External :: LA
      
       no=LEN_TRIM(ASYM(1:24))/3

       k=1; j=1
       Do i=1,no
        ln(i) = LA(ASYM(k:k))
        read(ASYM(k+1:k+2),*) iq(i) 
        k=k+3
        LS(i,3) = LA(BSYM(j:j))
        LS(i,2)=2*LA(BSYM(j+1:j+1))+1
        read(BSYM(j+2:j+2),'(i1)') LS(i,1)
        LS(i,5) = LA(BSYM(j+3:j+3))
        LS(i,4)=2*LA(BSYM(j+4:j+4))+1
        j=j+5
       End do

       End Subroutine Decode_conf



!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu1,nu2
      Integer, Intent(out) :: nc
      Integer :: i,j

      Integer,Parameter :: mc = 100000
      Integer,Allocatable,Dimension(:) :: K1,K2,K3
      Real(8),Allocatable,Dimension(:) :: C

      Allocate(C(mc),K1(mc),K2(mc),K3(mc))

      nc = 0
      i = 1
    1 read(nu1,end=2) c(i),k1(i),k2(i),k3(i)
      i = i + 1; if(i.le.mc) go to 1
    2 j = i - 1
      nc = nc + j

      Do i = 1,j
       write(nu2) c(i),k1(i),k2(i),k3(i)
      End do

      i = 1;  if(j.eq.mc) go to 1

      Deallocate(C,K1,K2,K3)

      End Subroutine RW







