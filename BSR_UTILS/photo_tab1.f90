!=======================================================================
!     Processing the photoionization data
!=======================================================================

      Use target; Use channels

      Implicit real(8) (A-H,O-Z)
 
      Allocatable IP_phot(:), NU_phot(:), E_thresh(:) 

      Integer :: ke=1000, ne=0, me
      Real(8), Allocatable :: e(:),ev(:)
      Integer, Allocatable :: IPE(:)
      Real(8), Allocatable :: CSR(:,:),CSV(:,:)
      Real(8), Allocatable :: ASR(:,:,:),ASV(:,:,:)
      Character(40) :: AF_res
      Real(8) :: E_shift = 0.0
      
! ... files:

      Character(40) :: targ  = 'target';       Integer :: nut = 1
      Character(40) :: rname = 'phot_res';     Integer :: nur = 2
      Character(40) :: pname = 'bsr_phot.nnn'; Integer :: nup = 20

      me=ke; Allocate(e(me),ev(me)) 

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nut,file=targ)
      Call R_target(nut)
      Call R_channels(nut)
      E1 = etarg(1); etarg = (etarg-e1)*2  ! -> Ry
      ion = nz-nelc

!----------------------------------------------------------------------
! ... choose the partial waves:

      write(*,*) 'Number of partial waves:'
      read(*,*) nphot
      Allocate(IP_phot(nphot), NU_phot(nphot), E_thresh(nphot))
      write(*,*) 'Their indexes according to file target:'
      read(*,*) IP_phot
      ii = INDEX(pname,'.',back=.true.)

      write(*,*) 'Energy shift?'
      read(*,*) E_shift
      
      Do i=1,nphot
       klsp = IP_phot(i)
       write(pname(ii+1:ii+3),'(i3.3)') klsp
       Call Check_file(pname)
       NU_phot(i) = nup + i
       Open(NU_phot(i),file=pname)       
       it = iptar(klsp,1)
       E_thresh(i) = etarg(it)
       Call Read_energies(NU_phot(i))
      End do

      write(*,*) 'ne = ',ne
      if(ne.eq.0) Stop ' ' 

      Call Rsort(ne,e) 
      Call Rsort(ne,ev) 
      Allocate(IPE(ne))
      Call sortR(ne,e,IPE) 
      ev = ev - E_shift 

!----------------------------------------------------------------------
! ... choose the work:

    1 Continue
      write(*,*) 'What do you want:'
      write(*,*) '1  - total cross sections'
      write(*,*) '2  - channel cross sections'
      write(*,*) '3  - ionic-state cross sections'
      write(*,*) '4  - beta-parameters'
      write(*,*) '5  - dipole matrix elements'
      write(*,*) '6  - delete points'
      write(*,*) '7  - check archives'
      write(*,*) '0  - exit'

      read(*,*) icase

      if(icase.gt.0.and.icase.le.4) then
       write(*,*) 'file-name for results:'
       read(*,'(a)') AF_res
      end if

      Select case(icase)
       Case(1);    Call Total_cs
       Case(2);    Call Channel_cs
       Case(3);    Call Target_cs
       Case(4);    Call Beta_par
       Case(5);    Call Dipole_me
       Case(6);    Call Delete
       Case(7);    Call Correct
       Case Default; go to 2
      End Select

      Go to 1
    2 Continue


CONTAINS


!=======================================================================
      Subroutine Read_energies(nu)
!=======================================================================

      Implicit real(8) (A-H,O-Z)
      
      Real(8), Allocatable :: ee(:)
      Integer :: i

      rewind(nu)
    1  read(nu,'(2d15.8,4i5)',end=2,err=2) e1,e2,nopen,nwt,ikm
       read(nu,'(5d15.8)',err=2) SLP,(SL,i=1,nopen)
       read(nu,'(5d15.8)',err=2) SVP,(SV,i=1,nopen)
       read(nu,'(5d15.8)',err=2) us ,(ui,i=1,nopen)
       Do i=1,nopen
        read(nu,'(4d15.8)',err=2) dlr,dli,dvr,dvi
       End do
       if(nwt.gt.nopen) then
        read(nu,'(d15.8)',err=2) CC
        read(nu,'(5d15.8)',err=2) (WTch,i=1,nwt)
        Do i=1,nopen
         read(nu,'(5d15.8)',err=2) (AK,j=1,nwt)
        End do
       end if
       if(ikm.gt.0) read(nu,'(6d13.6)',err=2) ((C,i=1,j),j=1,nopen)       

       ie=0; Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

       if(ie.eq.0) then; ne=ne+1; ie=ne; e(ie)=e1; ev(ie)=e2; end if

       if(ne.eq.me) then
        me = me + ke
        Allocate(ee(ne))
        ee(1:ne) = e(1:ne); Deallocate(e); Allocate(e(me)); e(1:ne)=ee(1:ne)
        ee(1:ne) = ev(1:ne); Deallocate(ev); Allocate(ev(me)); ev(1:ne)=ee(1:ne)
        Deallocate(ee)
       end if

       go to 1
     2 Continue


      End Subroutine Read_energies


!=======================================================================
      Subroutine Total_cs
!=======================================================================

      Implicit real(8) (A-H,O-Z)
      Real(8), allocatable ::  x(:), y1(:), y2(:), z1(:), z2(:)
 
      if(allocated(CSR)) Deallocate(CSR,CSV)
      Allocate(CSR(ne,nphot));  CSR = 0.d0
      Allocate(CSV(ne,nphot));  CSV = 0.d0

      if(allocated(x)) Deallocate(x,y1,y2,z1,z2)
      Allocate(x(ne),y1(ne),y2(ne),z1(ne),z2(ne))
      
      Do iphot=1,nphot
       
       nu = NU_phot(iphot); rewind(nu)

    1  read(nu,'(2d15.8,4i5)',end=2) e1,e2,nopen,nwt,ikm
       read(nu,'(5d15.8)') SLP,(SL,i=1,nopen)
       read(nu,'(5d15.8)') SVP,(SV,i=1,nopen)
       read(nu,'(5d15.8)') us ,(ui,i=1,nopen)
       Do i=1,nopen
        read(nu,'(4d15.8)') dlr,dli,dvr,dvi
       End do
       if(nwt.gt.nopen) then
        read(nu,'(d15.8)') CC
        read(nu,'(5d15.8)') (WTch,i=1,nwt)
        Do i=1,nopen
         read(nu,'(5d15.8)') (AK,j=1,nwt)
        End do
       end if
       if(ikm.gt.0) read(nu,'(6d13.6)') ((C,i=1,j),j=1,nopen)       

       ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

       if(ie.eq.0) Stop 'Total_cs: problems with energies'


       CSR(ie,iphot) = SLP
       CSV(ie,iphot) = SVP

       go to 1
     2 Continue


       nx = 0
       Do ie=1,ne
        if(abs(CSR(ie,iphot)).eq.0.d0) Cycle
        nx = nx + 1;  x(nx)=e(ie);  y1(nx) = CSR(ie,iphot)
                                    z1(nx) = CSV(ie,iphot)
       End do

write(*,*) 'nx =',nx,ne,iphot

       if(nx.eq.ne) Cycle


nu = 50+iphot
Do i=1,nx
write(nu,'(2E15.4)') x(i),y1(i)
End do

       k = 2
       Call LAGRN(K,nx,ne,x,e,y1,y2)

nu = 60+iphot
Do i=1,ne
write(nu,'(2E15.4)') e(i),y2(i)
End do

       Call LAGRN(K,nx,ne,x,e,z1,z2)

       CSR(:,iphot) = y2
       CSV(:,iphot) = z2

      End do  ! over iphot
!----------------------------------------------------------------------
      Open(nur,file=AF_res)

      Do i=1,ne; ie=IPE(i)

       CL=0.d0; CV=0.d0; ii = 1 
       Do iphot=1,nphot
        if(e(ie).le.E_thresh(iphot)) Cycle
        if(CSR(ie,iphot).eq.0.d0) ii =0
        CL = CL + CSR(ie,iphot)           
        CV = CV + CSV(ie,iphot)           
       End do

!       if(ii.eq.0) Cycle

if(maxval(CSR(ie,1:nphot)).gt.200) Cycle          !  ???
       
       write(nur,'(f16.8,f15.8,8D15.7)') ev(ie),e(ie), &
                CL,CV,(CSR(ie,j),CSV(ie,j),j=1,nphot)

      End do

      End Subroutine Total_cs


!=======================================================================
      Subroutine Channel_cs
!=======================================================================

      Implicit real(8) (A-H,O-Z)
 
      Real(8), Allocatable, Dimension(:) :: SL,SV

      write(*,*) 'partial wave and range of channels (0,0 means all):'
      read(*,*) iphot,ich1,ich2

      kch = nch(IP_phot(iphot))
      if(ich1.le.0) ich1=1
      if(ich2.le.0) ich2=kch
      ich = ich2-ich1+1
      if(ich.le.0) Return
      if(Allocated(SL)) Deallocate(SL,SV)
      Allocate(SL(kch),SV(kch))

      if(Allocated(CSR)) Deallocate(CSR,CSV)
      Allocate(CSR(ne,ich),CSV(ne,ich))
      CSR = 0.d0; CSV = 0.d0

      nu = NU_phot(iphot); rewind(nu)

    1  read(nu,'(2d15.8,4i5)',end=2) e1,e2,nopen,nwt,ikm
       SL = 0.d0; SV = 0.d0
       read(nu,'(5d15.8)') SLP,(SL(i),i=1,nopen)
       read(nu,'(5d15.8)') SVP,(SV(i),i=1,nopen)
       read(nu,'(5d15.8)') us ,(ui,i=1,nopen)
       Do i=1,nopen
        read(nu,'(4d15.8)') dlr,dli,dvr,dvi
       End do
       if(nwt.gt.nopen) then
        read(nu,'(d15.8)') CC
        read(nu,'(5d15.8)') (WTch,i=1,nwt)
        Do i=1,nopen
         read(nu,'(5d15.8)') (AK,j=1,nwt)
        End do
       end if
       if(ikm.gt.0) read(nu,'(6d13.6)') ((C,i=1,j),j=1,nopen)       

       ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

       if(ie.eq.0) Stop 'Channel_cs: problems with energies'

       CSR(ie,1:ich) = SL(ich1:ich2)
       CSV(ie,1:ich) = SV(ich1:ich2)
       ev(ie) = e2

       go to 1
     2 Continue

      Open(nur,file=AF_res)

      Do i=1,ne; ie=IPE(i)
       write(nur,'(f10.6,f15.8,20D12.4)') e(ie),ev(ie), &
                  (CSR(ie,j),CSV(ie,j),j=1,ich)
      End do

      End Subroutine Channel_cs


!=======================================================================
      Subroutine Target_cs
!=======================================================================

      Implicit real(8) (A-H,O-Z)
 
      Real(8), Allocatable :: SL(:),SV(:)

      write(*,*) 'range of target states (0,0 means all):'
      read(*,*) it1,it2

      if(it1.le.0) it1=1
      if(it2.le.0) it2=ntarg
      itt = it2-it1+1
      if(itt.le.0) Return

      if(Allocated(ASR)) Deallocate(ASR,ASV)
      Allocate(ASR(ne,nphot,itt),ASV(ne,nphot,itt))
      ASR = 0.d0; ASV = 0.d0

      Do iphot=1,nphot
       
       kch = nch(IP_phot(iphot))
       if(Allocated(SL)) Deallocate(SL,SV)
       Allocate(SL(kch),SV(kch))

       nu = NU_phot(iphot); rewind(nu)

    1  read(nu,'(2d15.8,4i5)',end=2) e1,e2,nopen,nwt,ikm
       read(nu,'(5d15.8)') SLP,(SL(i),i=1,nopen)
       read(nu,'(5d15.8)') SVP,(SV(i),i=1,nopen)
       read(nu,'(5d15.8)') us ,(ui,i=1,nopen)
       Do i=1,nopen
        read(nu,'(4d15.8)') dlr,dli,dvr,dvi
       End do
       if(nwt.gt.nopen) then
        read(nu,'(d15.8)') CC
        read(nu,'(5d15.8)') (WTch,i=1,nwt)
        Do i=1,nopen
         read(nu,'(5d15.8)') (AK,j=1,nwt)
        End do
       end if
       if(ikm.gt.0) read(nu,'(6d13.6)') ((C,i=1,j),j=1,nopen)       

       ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

       if(ie.eq.0) Stop 'Target_cs: problem with energies'

       Do i=1,nopen
        it = iptar(IP_phot(iphot),i)
        if(it.lt.it1) Cycle
        if(it.gt.it2) Cycle
        jt = it-it1+1
        ASR(ie,iphot,jt) = ASR(ie,iphot,jt)+SL(i)
        ASV(ie,iphot,jt) = ASV(ie,iphot,jt)+SV(i)
       End do

       go to 1
     2 Continue

      End do  ! over iphot

      Open(nur,file=AF_res)

      Do i=1,ne; ie=IPE(i)

       ii = 1 
       Do iphot=1,nphot
        if(e(ie).le.E_thresh(iphot)) Cycle
        S = SUM(ASR(ie,iphot,:))
        if(S.eq.0.d0) ii =0
       End do      
       if(ii.eq.0) Cycle

       SL=0.d0; SV=0.d0 
       Do it = 1,itt               
        SL(it) = SUM(ASR(ie,:,it))
        SV(it) = SUM(ASV(ie,:,it))
       End do

if(SL(1).gt.80) Cycle          !  ???
if(SL(2).gt.10) Cycle          !  ???

       write(nur,'(f10.6,f15.8,40D12.5)') e(ie),ev(ie), &
                 (SL(it),SV(it),it=1,itt)   

      End do

      End Subroutine Target_cs


!=======================================================================
      Subroutine Beta_par
!=======================================================================

      Implicit real(8) (A-H,O-Z)
 
      Real(8), Allocatable :: DLR(:),DLI(:),DVR(:),DVI(:), CP(:)
      Real(8), Allocatable :: ALR(:,:,:),ALI(:,:,:),AVR(:,:,:),AVI(:,:,:)
      Real(8), Allocatable :: betaL(:),betaV(:)

! ... ask additional parameters:

      write(*,*) 'Initiasl L:'
      read(*,*) L0

      write(*,*) 'range of target states (0,0 means all):'
      read(*,*) it1,it2

      if(it1.le.0.or.it1.gt.ntarg) it1=1
      if(it2.le.0.or.it2.gt.ntarg) it2=ntarg
      if(it2.lt.it1) it2=it1
      
      mch = 0
      Do i=1,nphot
       kch = nch(IP_phot(i)); if(kch.gt.mch) mch=kch
      End do

      if(Allocated(ALR)) Deallocate(ALR,ALI,AVR,AVI)
      Allocate(ALR(ne,nphot,mch),ALI(ne,nphot,mch), &
               AVR(ne,nphot,mch),AVI(ne,nphot,mch))
      ALR = 0.d0; ALI = 0.d0; AVR = 0.d0; AVI = 0.d0

      Do iphot=1,nphot
       
       kch = nch(IP_phot(iphot))
       if(Allocated(DLR)) Deallocate(DLR,DLI,DVR,DVI)
       Allocate(DLR(kch),DLI(kch),DVR(kch),DVI(kch))
       DLR=0.d0; DLI=0.d0; DVR=0.d0; DVI=0.d0

       nu = NU_phot(iphot); rewind(nu)

    1  read(nu,'(2d15.8,3i5)',end=2) e1,e2,nopen,nwt,ikm
       read(nu,'(5d15.8)') SLP,(SL,i=1,nopen)
       read(nu,'(5d15.8)') SVP,(SV,i=1,nopen)
       read(nu,'(5d15.8)') us ,(ui,i=1,nopen)

       if(nopen.gt.kch) Stop 'nopen > kch'

       Do i=1,nopen
        read(nu,'(4d15.8)') dlr(i),dli(i),dvr(i),dvi(i)
       End do
       if(nwt.gt.nopen) then
        read(nu,'(d15.8)') CC
        read(nu,'(5d15.8)') (WTch,i=1,nwt)
        Do i=1,nopen
         read(nu,'(5d15.8)') (AK,j=1,nwt)
        End do
       end if
       if(ikm.gt.0) read(nu,'(6d13.6)') ((C,i=1,j),j=1,nopen)       

       ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

       if(ie.eq.0) Stop 'Beta_par: problems with energies'

       Do i=1,nopen
        ALR(ie,iphot,i) = DLR(i)
        ALI(ie,iphot,i) = DLI(i)         ! ???     -   
        AVR(ie,iphot,i) = DVR(i)
        AVI(ie,iphot,i) = DVI(i)
       End do

       go to 1
     2 Continue

      End do  ! over iphot


!----------------------------------------------------------------------

      lmax = 0
      Do ilsp=1,nphot
       klsp=IP_phot(ilsp); kch=nch(klsp) 
       Do ich=1,kch
        l=lch(klsp,ich)
        if(l.gt.lmax) lmax=l
       End do
      End do
      if(Allocated(CP)) Deallocate(CP); Allocate(CP(0:lmax))
      if(allocated(betaL)) Deallocate(betaL,betaV)
      Allocate(betaL(ntarg),betaV(ntarg))

      Open(nur,file=AF_res)

! ... enrgy loop:

      Do i=1,ne; ie=IPE(i)

       ii = 1 
       Do iphot=1,nphot
        if(e(ie).le.E_thresh(iphot)) Cycle
        S = SUM(ALR(ie,iphot,:))
        if(S.eq.0.d0) ii =0
       End do      
       if(ii.eq.0) Cycle

! ... Coulomb phases: cp(l+1) = s(l+1)-s(0) = cp(l) + tn^-1(q/(l+1)):

       CP = 0.d0
       if(ion.gt.0) then
        q = ion/sqrt(e(ie))
        Do l=1,lmax
         s=l; CP(l)=CP(l-1)+DATAN2(q,s)
        End do    
       end if

! ... target loop:
       
       betaL = 0.d0; betaV= 0.d0

       Do it = it1,it2

        LF = ltarg(it)

        DENL = 0.d0; DENV = 0.d0

        Do ilsp1=1,nphot; klsp1=IP_phot(ilsp1); kch1=nch(klsp1) 
        Do ich1=1,kch1; if(iptar(klsp1,ich1).ne.it) Cycle
           ILT1 = lpar(klsp1); l1 = lch(klsp1,ich1)
           a1 = ALR(ie,ilsp1,ich1); b1 = ALI(ie,ilsp1,ich1)
           c1 = AVR(ie,ilsp1,ich1); d1 = AVI(ie,ilsp1,ich1)

           denL = denL + (a1*a1+b1*b1) 
           denV = denV + (c1*c1+d1*d1) 

        Do ilsp2=1,nphot; klsp2=IP_phot(ilsp2); kch2=nch(klsp2) 
        Do ich2=1,kch2; if(iptar(klsp2,ich2).ne.it) Cycle
           ILT2 = lpar(klsp2); l2 = lch(klsp2,ich2)        
           a2 = ALR(ie,ilsp2,ich2); b2 = -ALI(ie,ilsp2,ich2)
           c2 = AVR(ie,ilsp2,ich2); d2 = -AVI(ie,ilsp2,ich2)

          S = (l1+l1+1)*(l2+l2+1)*(ILT1+ILT1+1)*(ILT2+ILT2+1)
          S = sqrt(S)
          S = S * Z_6jj(l1,l2,2,ILT2,ILT1,LF) &  
                * Z_6jj( 1, 1,2,ILT2,ILT1,L0) & 
                * Z_3j0(l1,l2,2)              &
                * Z_3j0( 1, 1,2)

          if(S.eq.0.d0) Cycle

          aa = a1*a2 - b1*b2; bb = a1*b2 + b1*a2
          cc = c1*c2 - d1*d2; dd = c1*d2 + d1*c2

          kz = (l1-l2)/2
          kz = (-1)**kz
!kz = 1
          aa = aa*kz; bb = bb*kz; cc = cc*kz; dd = dd*kz

          if(ion.ne.0) then
           Ss = CP(l2)-CP(l1)
           sr=cos(Ss); si = sin(Ss)
           a = aa*sr - bb*si; b = aa*si + bb*sr
           c = cc*sr - dd*si; d = cc*si + dd*sr
          else
           a=aa; b=bb; c=cc; d=dd
          end if


          betaL(it) = betaL(it) + a * S
          betaV(it) = betaV(it) + c * S


if(abs(e(ie)-0.1d0).lt.0.0001d0) then
write(51,'(a,6i5,3f10.5)') 'LO,LF,ILT1,ILT2,l1,l2,S', &
                     L0,LF,ILT1,ILT2,l1,l2,S*15,a,betaL(it)*15
end if


        End do; End do
        End do; End do

          if(DENL.eq.0.d0) Cycle

          S = 15 * (-1)**(L0+LF)

          betaL(it) = betaL(it) / DENL * S  
          betaV(it) = betaV(it) / DENV * S

       End do  ! over target


       write(nur,'(f10.6,f15.8,8D15.7)') e(ie),ev(ie), &
                (BETAL(it),BETAV(it),it=it1,it2)

      End do   ! over energy

      End Subroutine Beta_par



!=======================================================================
      Subroutine Dipole_me
!=======================================================================

      Implicit real(8) (A-H,O-Z)
 
      Real(8), Allocatable :: DLR(:),DLI(:),DVR(:),DVI(:)
      Real(8), Allocatable :: dd(:,:)

      write(*,*) 'partial wave and channel: '
      read (*,*) iphot,ich

      kch = nch(IP_phot(iphot))

      if(Allocated(DLR)) Deallocate(DLR,DLI,DVR,DVI)

      Allocate(DLR(kch),DLI(kch),DVR(kch),DVI(kch))
      DLR=0.d0; DLI=0.d0; DVR=0.d0; DVI=0.d0

      if(allocated(dd)) deallocate(dd); Allocate(dd(8,ne)); dd = 0.d0

      nu = NU_phot(iphot); rewind(nu)

   1  read(nu,'(2d15.8,3i5)',end=2) e1,e2,nopen,nwt,ikm
      if(nopen.gt.kch) Stop 'nopen > kch'
      read(nu,'(5d15.8)') SLP,(SL,i=1,nopen)
      read(nu,'(5d15.8)') SVP,(SV,i=1,nopen)
      read(nu,'(5d15.8)') us ,(ui,i=1,nopen)

      Do i=1,nopen
       read(nu,'(4d15.8)') dlr(i),dli(i),dvr(i),dvi(i)
      End do

      if(nwt.gt.nopen) then
       read(nu,'(d15.8)') CC
       read(nu,'(5d15.8)') (WTch,i=1,nwt)
       Do i=1,nopen
        read(nu,'(5d15.8)') (AK,j=1,nwt)
       End do
      end if
      if(ikm.gt.0) read(nu,'(6d13.6)') ((C,i=1,j),j=1,nopen)       

      if(ich.gt.nopen) go to 1

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

      if(ie.eq.0) Stop 'Dipole_me: problems with energies'

       dd(1,ie) = DLR(ich)
       dd(3,ie) = DLI(ich)         ! ???     -   
       dd(5,ie) = DVR(ich)
       dd(7,ie) = DVI(ich)

      go to 1
    2 Continue

      h1 = e(2)-e(1); hh1=h1*h1
      Do i=2,ne-1
       h = e(i+1)-e(i); hh = h*h
       dd(2,i) = (hh1*dd(1,i+1)+(hh-hh1)*dd(1,i)-hh*dd(1,i-1))/(h*h1*(h+h1))
       dd(4,i) = (hh1*dd(3,i+1)+(hh-hh1)*dd(3,i)-hh*dd(3,i-1))/(h*h1*(h+h1))
       dd(6,i) = (hh1*dd(5,i+1)+(hh-hh1)*dd(5,i)-hh*dd(5,i-1))/(h*h1*(h+h1))
       dd(8,i) = (hh1*dd(7,i+1)+(hh-hh1)*dd(7,i)-hh*dd(7,i-1))/(h*h1*(h+h1))
       h1 = h; hh1 = hh
      End do

      write(AF_res,'(a,i3.3,a,i3.3)') 'dme_',iphot,'_',ich

      Open(nur,file=AF_res)

      Do ie=2,ne-1
       write(nur,'(10E16.8)') e(ie),ev(ie), dd(1:4,ie)
      End do

      End Subroutine Dipole_me


!=======================================================================
      Subroutine Delete
!=======================================================================

      Implicit real(8) (A-H,O-Z)
 
      Real(8), Allocatable :: S(:),SS(:,:),ed(:)

! ... find delete energies:

      nud = 3
      open(nud,file='delete_list')
      nd = 0
   10 read(nud,*,end=20) x   
      nd = nd+1
      go to 10
   20 write(*,*) 'ndelete = ',nd
      if(allocated(ed)) Deallocate(ed); Allocate(ed(nd))
      rewind(nud)
      Do i=1,nd; read(nud,*) ed(i); End do
      Close(nud)
      

      Do iphot=1,nphot

       kch = nch(IP_phot(iphot))
       if(Allocated(S)) Deallocate(S,SS)
       Allocate(S(kch),SS(kch,kch))
       ndim = kch

       nu = NU_phot(iphot); rewind(nu)

       ii = INDEX(pname,'.',back=.true.)
       klsp = IP_phot(iphot)
       write(pname(ii+1:),'(i3.3,a)') klsp,'_del'
       Open(nud,file=pname)       

    1  read(nu,'(2d15.8,4i5)',end=2) e1,e2,nopen,nwt,ikm
       idel = 0; Do i=1,nd; if(ed(i).ne.e1) Cycle; idel=1; Exit; End do 
       if(idel.eq.0) write(nud,'(2d15.8,4i5)') e1,e2,nopen,nwt,ikm      
       read(nu,'(5d15.8)') SLP,(S(i),i=1,nopen)
       if(idel.eq.0) write(nud,'(5d15.8)') SLP,(S(i),i=1,nopen)
       read(nu,'(5d15.8)') SVP,(S(i),i=1,nopen)
       if(idel.eq.0) write(nud,'(5d15.8)') SVP,(S(i),i=1,nopen)
       read(nu,'(5d15.8)') us,(S(i),i=1,nopen)
       if(idel.eq.0) write(nud,'(5d15.8)') us,(S(i),i=1,nopen)

       Do i=1,nopen
        read(nu,'(4d15.8)') dlr,dli,dvr,dvi
        if(idel.eq.0) write(nud,'(4d15.8)') dlr,dli,dvr,dvi
       End do

       if(nwt.gt.nopen) then
        read(nu,'(d15.8)') CC
        if(idel.eq.0) write(nud,'(d15.8)') CC
        if(nwt.gt.ndim) then
         Deallocate(S); Allocate(S(nwt)); ndim=nwt
        end if
        read(nu,'(5d15.8)') (S(i),i=1,nwt)
        if(idel.eq.0) write(nud,'(5d15.8)') (S(i),i=1,nopen)
        Do i=1,nopen
         read(nu,'(5d15.8)') (S(j),j=1,nwt)
        if(idel.eq.0) write(nud,'(5d15.8)') (S(j),j=1,nopen)
        End do
       end if
       if(ikm.gt.0) then
        read(nu,'(6d13.6)') ((SS(i,j),i=1,j),j=1,nopen)       
        if(idel.eq.0) write(nud,'(6d13.6)') ((SS(i,j),i=1,j),j=1,nopen)
       end if
       go to 1
     2 Continue

      End do  ! over iphot
      Deallocate(S,SS)

      End Subroutine Delete



!=======================================================================
      Subroutine Correct
!=======================================================================
      Implicit real(8) (A-H,O-Z)
 
      Character(80) :: AS
      Real(8), allocatable :: S(:),SS(:,:),ed(:)

      Do iphot=1,nphot

       kch = nch(IP_phot(iphot))
       if(Allocated(S)) Deallocate(S,SS)
       Allocate(S(kch),SS(kch,kch))
       ndim = kch

       nu = NU_phot(iphot); rewind(nu)

       ii = INDEX(pname,'.',back=.true.)
       klsp = IP_phot(iphot)
       write(pname(ii+1:),'(i3.3,a)') klsp,'_corrected'
       nud=3; Open(nud,file=pname)       

       
    1  read(nu,'(a)',end=2) AS
       if(len_trim(AS).ne.45) go to 1
       if(AS(36:36).ne.' ') go to 1
       read(AS,'(2d15.8,4i5)',err=1) e1,e2,nopen,nwt,ikm
       if(nopen.le.0.or.nopen.gt.kch) go to 1
       write(nud,'(2d15.8,4i5)') e1,e2,nopen,nwt,ikm      
       read(nu,'(5d15.8)',err=1) SLP,(S(i),i=1,nopen)
       write(nud,'(5d15.8)') SLP,(S(i),i=1,nopen)
       read(nu,'(5d15.8)',err=1) SVP,(S(i),i=1,nopen)
       write(nud,'(5d15.8)') SVP,(S(i),i=1,nopen)
       read(nu,'(5d15.8)',err=1) us,(S(i),i=1,nopen)
       write(nud,'(5d15.8)') us,(S(i),i=1,nopen)

       Do i=1,nopen
        read(nu,'(4d15.8)',err=1) dlr,dli,dvr,dvi
        write(nud,'(4d15.8)') dlr,dli,dvr,dvi
       End do

       if(nwt.gt.nopen) then
        read(nu,'(d15.8)',err=1) CC
        write(nud,'(d15.8)') CC
        if(nwt.gt.ndim) then
         Deallocate(S); Allocate(S(nwt)); ndim=nwt
        end if
        read(nu,'(5d15.8)',err=1) (S(i),i=1,nwt)
        write(nud,'(5d15.8)') (S(i),i=1,nopen)
        Do i=1,nopen
         read(nu,'(5d15.8)',err=1) (S(j),j=1,nwt)
         write(nud,'(5d15.8)') (S(j),j=1,nopen)
        End do
       end if
       if(ikm.gt.0) then
        read(nu,'(6d13.6)',err=1) ((SS(i,j),i=1,j),j=1,nopen)       
        write(nud,'(6d13.6)') ((SS(i,j),i=1,j),j=1,nopen)
       end if
       go to 1
     2 Continue

      End do  ! over iphot
      Deallocate(S,SS)

      End Subroutine Correct


      End !  utility phot_tab
