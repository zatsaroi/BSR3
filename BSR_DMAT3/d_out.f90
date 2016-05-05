!=======================================================================
      Subroutine D_OUT
!=======================================================================
!     define and output dipole matrix for R-matrix solutions:
!-----------------------------------------------------------------------
      Use bsr_dmat 
      Use dmatrix
      Use target, only: etarg

      Implicit none
      Character(64) :: Label1,AS
      Real(8), allocatable ::  CL(:),CV(:)
      Real(8), allocatable ::  eval(:)
      Real(8) :: S,SL,SV,E1,E2
      Integer :: i,j,ns1,nhm1,nhm,khm,kch,kcp,isol1,isol2
      Integer :: JLT1,JST1,JLT2,JST2
      Integer, external :: Idef_st 

      if(jmode.eq.1) Stop 'D_out: jmode = 1 ' 
      
      i = INDEX(AF_d,'.'); AF = AF_d(1:i)//ALS2
      Open(nud,file=AF,form='UNFORMATTED')

!----------------------------------------------------------------------
! ... define the initial bound states:

      if(ctype1.eq.'b') then

       AF='bound.'//ALS1; Open(inb1,file=AF,status='OLD'); rewind(inb1)

       read(inb1,*) ns1,nch1,ncp1,nhm1,nstate1

       if(ns1 .ne.kns ) Stop ' ns1 <> ns '
       if(nch1.ne.kch1) Stop ' nch1 --> ?'
       if(ncp1.ne.kcp1) Stop ' ncp1 --> ?'
       if(nhm1.ne.kdm1) Stop ' nhm1 --> ? '

      elseif(ctype1.eq.'c') then
	   
	      inb1=in1; nstate1=1

      else

       i = LEN_TRIM(name1)-1; AF = name1(1:i)//ctype1
       Open(inb1,file=AF,status='OLD')
       nstate1 = Idef_st(inb1)                                  

      end if

!-----------------------------------------------------------------------
! ... read final states: inner region R-matrix solutions

      i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALS2
      Call Check_file(AF)
      Open(nuo,file=AF,form='UNFORMATTED',STATUS='OLD')
      read(nuo) nhm,khm,kch,kcp
      if(nhm.ne.kdm2) Stop 'D_out, rsol: nhm <> kdm2 '
      if(kch.ne.kch2) Stop 'D_out, rsol: kch <> kch2 '
      if(kcp.ne.kcp2) Stop 'D_out, rsol: kcp <> kcp2 '
      Allocate(eval(khm)) 
      read(nuo) eval
      nstate2 = khm

!----------------------------------------------------------------------
! ... calculation and output the dipole matrix:

      Allocate(CL(nstate2),CV(nstate2))

      Deallocate(C1,C2); Allocate(C1(kdm1),C2(kdm2))

      if(jmode.eq.0) then
       JLT1 = (ILT1-1)/2; JST1 = IST1;
       JLT2 = (ILT2-1)/2; JST2 = IST2;
      else
       JLT1 = JOT1-1; JLT2 = JOT2-1; JST1=0; JST2=0
      end if
      E2 = Etarg(1)
      write(nud) JLT2,JST2,IPT2,E2,nstate2                          

! ... read initial state:

      if(istate1.le.0.or.istate1.gt.nstate1) istate1=1
      rewind(inb1);  read(inb1,'(a)') AS
      Do isol1=1,istate1
       Call Read_sol(ctype1,inb1,kdm1,C1,Label1,E1,jot1)
      End do

! ... loop over final set:

      Do isol2=1,nstate2
       read(nuo) C2
       SL=0.d0; SV=0.d0
       Do j=1,kdm2
        S = SUM(DV(1:kdm1,j)*C1(1:kdm1));  SV = SV + C2(j)*S
        S = SUM(DL(1:kdm1,j)*C1(1:kdm1));  SL = SL + C2(j)*S
       End do
       CL(isol2)=SL; CV(isol2)=SV
      End do

      write(nud) JLT1,JST1,IPT1,E1                   
      write(nud) (cl(i),cv(i),i=nstate2,1,-1)        

      Deallocate(CL,CV)

      End Subroutine D_out