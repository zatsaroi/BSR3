!=======================================================================
      Subroutine dvec_out
!=======================================================================
!     define and output dipole vector for initial state
!-----------------------------------------------------------------------
      Use bsr_dmat
      Use dmatrix

      Implicit none
      Character(64) :: Label1
      Integer :: i,j,ns1,nhm1
      Real(8) :: E1,SL,SV 
      Real(8), allocatable ::  CL(:),CV(:)
      Integer, external :: Idef_st
 
      if(ktype.ne.'E') Stop 'dv_out: non-electric-transition case ? '
      if(kpol.ne.1)    Stop 'dv_out: kpol <> 1 --> non-dipole case ? '
      if(jmode.eq.1)   Stop 'dv_out: jmode = 1' 
      
      i = INDEX(AF_v,'.'); AF = AF_v(1:i)//ALS2
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
!----------------------------------------------------------------------
! ... calculation and output the dipole matrix:

      Allocate(CL(kdm2),CV(kdm2))
      Deallocate(C1); Allocate(C1(kdm1))

! ... define initial state:

      if(istate1.le.0.or.istate1.gt.nstate1) istate1=1
      rewind(inb1);  read(inb1,'(a)') Label1
      Do i=1,istate1
       Call Read_sol(ctype1,inb1,kdm1,C1,Label1,E1,jot1)
      End do

      if(jot1.eq.0) jot1 = ILT1

      SL=0.d0; SV=0.d0
      Do j=1,kdm2
       CV(j) = SUM(DV(1:kdm1,j)*C1(1:kdm1))
       CL(j) = SUM(DL(1:kdm1,j)*C1(1:kdm1))
      End do

      write(nud) kpol,ktype
      write(nud) E1,jot1,IPT1,Label1
      write(nud) kdm2,kch2,kcp2                  
      write(nud) CL
      write(nud) CV

      Deallocate(CL,CV)

      End Subroutine dvec_out
