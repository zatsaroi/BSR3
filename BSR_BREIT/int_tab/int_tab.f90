!=====================================================================
!     UTILITY   I N T L S T
!
!               C O P Y R I G H T -- 2006
!
!     Written by:   Oleg Zatsarinny,   email: oleg_zoi@yahoo.com
!                         
!======================================================================
!
!     records to the file INT_LST the integral coefficients 
!     for channels interaction based on:
!
!     channels.lst  - channel information
!     chan_even.c   - configurations for even states
!     chan_odd.c    - configurations for odd states 
!     chan_even.bnk - angular data bank for even states
!     chan_odd.bnk  - angular data bank for odd states 
!     mult_bnk      - angular data bank for dipole interaction
!     int_lst       - angular coefficients for given case
!     
!----------------------------------------------------------------------

!=====================================================================
      Module param_intlst      
!=====================================================================
!
!     Contains specific parameters for given problem 
!
!---------------------------------------------------------------------

      Implicit none
      Save

! ... files:

      Integer(4) :: inp = 10; Character(80) :: AF_inp  = 'channels.lst'
      Integer(4) :: nco = 11; Character(80) :: AF_odd  = 'chan_odd.c'
      Integer(4) :: nce = 12; Character(80) :: AF_even = 'chan_even.c'
      Integer(4) :: nbo = 13; Character(80) :: AF_bnko = 'chan_odd.bnk'
      Integer(4) :: nbe = 14; Character(80) :: AF_bnke = 'chan_even.bnk'
      Integer(4) :: nbm = 15; Character(80) :: AF_mult = 'mult_bnk'
      Integer(4) :: out = 16; Character(80) :: AF_out  = 'channels.int'
      
! ... tollerance for coefficients:

      Real(8) :: eps_c = 1.d-7

      Integer(4) :: nch,nch1,nch2,ncomp
      Integer(4) :: nwt, nwc

      End module param_intlst

!======================================================================
!     MAIN 
!======================================================================

      Use param_intlst
      Use configs
      Use zoef_list

!----------------------------------------------------------------------
! ... open files:

      Open(inp,file=AF_inp,status='OLD')
      Open(nco,file=AF_odd,status='OLD')
      Open(nce,file=AF_even,status='OLD')
      Open(nbo,file=AF_bnko,form='UNFORMATTED',status='OLD')
      Open(nbe,file=AF_bnke,form='UNFORMATTED',status='OLD')
      Open(nbm,file=AF_mult,form='UNFORMATTED',status='OLD')
      Open(out,file=AF_out)

      Call Read_rarg('c_eps',c_eps)
!----------------------------------------------------------------------
! ... read list of one-electron orbitals:

      Call Read_ipar(inp,'nwt',nwt)
      Call Read_ipar(inp,'nwc',nwc)
      nwf = nwt + nwc

      m=nwf; Call allocate_orblist(m)

      Call Read_ipar(inp,'nwt',nwt)
      read(inp,*)
      read(inp,'(20a4)') (ELF(i),i=1,nwt)
      Call Read_ipar(inp,'nwc',nwc)
      read(inp,*)
      read(inp,'(20a4)') (ELF(i),i=nwt+1,nwf)
      Do i=1,nwf
       Call EL4_NLK(ELF(i),NEF(i),LEF(i),KEF(i))
      End do

      JORT=-1; Call Pre_iort(0)

      Call Read_ipar(inp,'nch1',nch1)
      Call Read_ipar(inp,'nch2',nch2)
      nch = nch1+nch2

!----------------------------------------------------------------------
! ... processing odd configurations:

      Call R_conf(nco,nbo)
      Call R_integrals(nbo)

! ... processing even configurations:

      Call R_conf(nce,nbe)
      Call R_integrals(nbe)

! ... dipole interaction:

      Call RR_conf(nco,nce,nbm,ncfg1,ncfg2)
      Call D_integrals(nbm)

! ... output results:

      Call Output_integrals  

      End ! program INT_TAB
 


!======================================================================
      Subroutine R_integrals(nub)
!======================================================================
!
!     extracts the data from INT.BNK for specific configurations
!
!----------------------------------------------------------------------

      Use param_intlst
      USE configs 

      IMPLICIT REAL(8) (A-H,O-Z)

      Integer(4), parameter :: idet = 2**15, idef = 2**4 
      Integer(4), Dimension(:), Allocatable :: KPD,IPD,NPD
      Integer(4), Dimension(:), Allocatable :: KPF,IPF,NPF

      Integer(4), Dimension(ne) :: IPN1,IPN2,NP1,NP2,NP

!----------------------------------------------------------------------
!                                             read determinant factors:
      read(nub) ndet,kdet
      if(allocated(KPD)) Deallocate(KPD,IPD,NPD)
      Allocate(KPD(ndet),IPD(ndet),NPD(kdet))
      read(nub) KPD; read(nub) IPD; read(nub) NPD

      read(nub) ndf,kdf
      if(allocated(KPF)) Deallocate(KPF,IPF,NPF)
      Allocate(KPF(ndf),IPF(ndf),NPF(kdf))
      read(nub) KPF; read(nub) IPF; read(nub) NPF

!----------------------------------------------------------------------
!                                                  processing the data:
    1 read(nub,end=2) C,ijt,int,idf

      Call Decode_int(icase,kpol,I1,I2,I3,I4,int)

      if(icase.ne.5.and.icase.ne.6) go to 1
      if(icase.eq.6) C = -2.d0*C

! ... determine the range of states for given coeff. 

      it = ijt/ibc;   jt = mod(ijt,ibc)
      is1 = IT_state1(it); js1 = IT_state1(jt)
      is2 = IT_state2(it); js2 = IT_state2(jt)
      if(is1.eq.0.or.js1.eq.0) go to 1

      kd = 0
      if(idf.gt.0) then
       kd=KPF(idf)
       ip=IPF(idf)
       NP(1:kd)=NPF(ip+1:ip+kd)
      end if
      
! ... loop over all relevant states: 

      Do ik=is1,is2
       is=IP_stat(ik); noi=NOCCSH(is); IPN1(1:noi)=NOCORB(1:noi,is)

      Do jk=js1,js2
       js=IP_stat(jk); noj=NOCCSH(js); IPN2(1:noj)=NOCORB(1:noj,js)

! ... consider only low-half of interaction matrix

       if(it.eq.jt.and.js.gt.is) Cycle                   

! ... define integral for specific orbitals

       j1=IPN1(i1); j2=IPN1(i2); j3=IPN2(i3); j4=IPN2(i4)

! ... define the channel index: 

       ic=IT_stat(is); jc=IT_stat(js)

! ... target energies

      if(ic.eq.jc.and.ic.le.nch.and.jc.le.nch) then
       if(j1.le.nwt.and.j2.le.nwt.and.j3.le.nwt.and.j4.le.nwt) go to 1 
      end if

! ... include the expansion coefficients 

       CC = C * WC(is) * WC(js)
       if(is.ne.js.and.ic.eq.jc) CC = CC + CC
       if(abs(CC).lt.Eps_C) Cycle

! ... find overlaps for specific orbitals with extracted continuum:  

       if(idf.gt.0) then
        kz=0; jj = 1
        Do ii = 1,kd
         id=NP(ii)/idef; kdn=KPD(id); ip=IPD(id)
         Do i=1,kdn
          k=NPD(i+ip); NP1(i)=IPN1(k/idet); NP2(i)=IPN2(mod(k,idet))
         End do
         jj=IDET_SIMP(kz,kdn,NP1,NP2)
         if(jj.eq.0) Exit
         if(jj.eq.2) Stop 'IDET_SIMP = 2'
        End do 
        if(jj.eq.0) Cycle
        CC = CC * (-1)**kz
       end if

       if(icase.eq.5) Call Add_rk_data(CC,ic,jc,kpol,j1,j2,j3,j4)
       if(icase.eq.6) Call Add_L_data(CC,ic,jc,j1,j3)

      End do    !  over js
      End do	!  over is

      go to 1
    2 Continue

      End Subroutine R_integrals


!======================================================================
      Subroutine D_integrals(nub)
!======================================================================
!
!     This routine extract the dipole matrix elements
!     from MULT_BNK (unit nub) for specific configurations
!
!----------------------------------------------------------------------

      Use param_intlst
      Use configs

      Implicit real(8) (A-H,O-Z)

      Integer(4), parameter :: idet = 2**15, idef = 2**4 
      Integer(4), Dimension(:), Allocatable :: KPD,IPD,NPD
      Integer(4), Dimension(:), Allocatable :: KPF,IPF,NPF

      Integer(4), Dimension(ne) :: NP, NP1,NP2, IPN1,IPN2

!----------------------------------------------------------------------
!                                             read determinant factors:
      read(nub) ndet,kdet
      if(allocated(KPD)) Deallocate(KPD,IPD,NPD)
      Allocate(KPD(ndet),IPD(ndet),NPD(kdet))
      read(nub) KPD; read(nub) IPD; read(nub) NPD

      read(nub) ndf,kdf
      if(allocated(KPF)) Deallocate(KPF,IPF,NPF)
      Allocate(KPF(ndf),IPF(ndf),NPF(kdf))
      read(nub) KPF; read(nub) IPF; read(nub) NPF

!----------------------------------------------------------------------
!                                                  processing the data:
    1 read(nub,end=2) C, ijt,int,idf

      it=ijt/ibc;  jt=mod(ijt,ibc)

      is1=IT_state1(it);  is2=IT_state2(it)
      js1=IT_state1(jt);  js2=IT_state2(jt)

      if(is1.eq.0.or.js1.eq.0) go to 1   ! no relevant states

      kd = 0
      if(idf.gt.0) then
       kd=KPF(idf)
       ip=IPF(idf)
       NP(1:kd)=NPF(ip+1:ip+kd)
      end if

      Call Decode_mult(icase,i1,i2,int)

!----------------------------------------------------------------------
!                                                    cycle over states:
      Do ik=is1,is2
       is=IP_stat(ik); no1=NOCCSH(is); IPN1(1:no1)=NOCORB(1:no1,is)

      Do jk=js1,js2
       js=IP_stat(jk); no2=NOCCSH(js); IPN2(1:no2)=NOCORB(1:no2,js)

       if(it.eq.jt.and.ic.gt.jc) Cycle 

! ... include the expansion coefficients 

       CC = C * WC(is) * WC(js)

       if(abs(CC).lt.Eps_C) Cycle

       isign=1
       ic=IT_stat(is); jc=IT_stat(js)
       j1 = IPN1(i1); j2 = IPN2(i2)
       if(jc.gt.ic) then
        kz = ( ILterm(is) - ILterm(js) + ISterm(is) - ISterm(js) )/2 
        isign=(-1)**kz
        ic = IT_stat(js); jc=IT_stat(is)
        j1 = IPN2(i2); j2 = IPN1(i1)
       end if

! ... check the overlap factor:  

       if(idf.gt.0) then
        kz=0; jj = 1
        Do ii = 1,kd
         id=NP(ii)/idef; kdn=KPD(id); ip=IPD(id)
         Do i=1,kdn
          k=NPD(i+ip); NP1(i)=IPN1(k/idet); NP2(i)=IPN2(mod(k,idet))
         End do
         jj=IDET_SIMP(kz,kdn,NP1,NP2)
         if(jj.eq.0) Exit
         if(jj.eq.2) Stop 'IDET_SIMP = 2'
        End do 
        if(jj.eq.0) Cycle
        CC = CC * (-1)**kz
       end if

       Call Add_D_data(CC,ic,jc,j1,j2)

      End do   ! over jc
      End do   ! over ic
      
      go to 1
    2 Continue

    
      End Subroutine D_integrals



!======================================================================
      Subroutine Output_integrals
!======================================================================

      Use configs; Use rk_data; Use L_data; Use D_data
      Use param_intlst

      IMPLICIT REAL(8) (A-H,O-Z)

! ... L-data:

      m = nlk
      Do i = 1,nlk
       if(abs(clk(i)).lt.eps_c) m=m-1
      End do

      write(out,'(i6,a)') m,'   -   number of L-integrals '       
      write(*,'(i6,a)') m,'   -   number of L-integrals '       
      Do i = 1,nlk
       if(abs(clk(i)).lt.eps_c) Cycle
!       ic = kl1(i)/ibl; jc = mod(kl1(i),ibl)        
       Call Get_ij(kl1(i),ic,jc)
       i1 = kl2(i)/ibl; i2 = mod(kl2(i),ibl)
       write(out,'(f16.8,2i6,3x,a2,a4,a1,a4,a1)') &
        clk(i),ic,jc,'L(',ELF(i1),';',ELF(i2),')'
      End do
       
! ... Rk-data:

      m = nrk
      Do i = 1,nrk
       if(abs(crk(i)).lt.eps_c) m=m-1
      End do

      write(out,'(i6,a)') m,'   -   number of Rk-integrals '       
      write(*,'(i6,a)') m,'   -   number of Rk-integrals '       
      Do i = 1,nrk
       if(abs(crk(i)).lt.eps_c) Cycle
!       ic = kr1(i)/ibr; jc = mod(kr1(i),ibr)        
       Call Get_ij(kr1(i),ic,jc)
       k  = kr2(i)
       i1 = kr3(i)/ibr; i2 = mod(kr3(i),ibr)
       i3 = kr4(i)/ibr; i4 = mod(kr4(i),ibr)
       write(out,'(f16.8,2i6,3x,a1,i2,a1,4(a4,a1))') &
        crk(i),ic,jc,'R',k,'(',ELF(i1),',',ELF(i2),';', &
                               ELF(i3),',',ELF(i4),')'
      End do

! ... D-data:

      m = ndk
      Do i = 1,ndk
       if(abs(cdk(i)).lt.eps_c) m=m-1
      End do

      write(out,'(i6,a)') m,'   -   number of D-integrals '       
      write(*,'(i6,a)') m,'   -   number of D-integrals '       
      Do i = 1,ndk
       if(abs(cdk(i)).lt.eps_c) Cycle
!       ic = kd1(i)/ibd; jc = mod(kd1(i),ibd)        
       Call Get_ij(kd1(i),ic,jc)
       i1 = kd2(i)/ibd; i2 = mod(kd2(i),ibd)
       write(out,'(f16.8,2i6,3x,a2,a4,a1,a4,a1)') &
        cdk(i),ic,jc,'d(',ELF(i1),';',ELF(i2),')'
      End do

      End Subroutine Output_integrals


!======================================================================
      Subroutine Get_ij(k,i,j)
!======================================================================

      i = sqrt(real(2*k))
      j = k - i*(i-1)/2
      if(j.gt.i) then
       i=i+1
       j = k - i*(i-1)/2
      end if

      End Subroutine Get_ij






