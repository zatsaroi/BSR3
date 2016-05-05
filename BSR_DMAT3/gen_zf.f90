!======================================================================
      Subroutine Gen_zf
!======================================================================
! ... define the number of solutions and open relevant file,
! ... then call LS or LSJ routines:  LS_case and J_case
!----------------------------------------------------------------------
      USE bsr_dmat
      USE dmatrix

      Implicit none
      Integer :: i,ns1,ns2,nhm1,nhm2
      Integer, external ::  Idef_st     

!----------------------------------------------------------------------
! ... define the number of solutions and open relevant file:

      if(ctype1.eq.'b') then

       AF='bound.'//ALS1
       Open(inb1,file=AF,status='OLD',action='READ')
       rewind(inb1)

       read(inb1,*) ns1,nch1,ncp1,nhm1,nstate1

       if(ns1 .ne.kns ) Stop 'gen_zf:  ns1 <> ns '
       if(nch1.ne.kch1) Stop 'gen_zf:  nch1 --> ?'
       if(ncp1.ne.kcp1) then 
        write(pri,*) 'gen_zf: ncp1, kcp1 =', ncp1,kcp1
        Stop 'gen_gf:  ncp1 --> ?'
       end if
       if(nhm1.ne.kdm1) Stop 'gen_zf:  nhm1 --> ? '

      elseif(ctype1.eq.'c') then
	   
       inb1=in1; nstate1=1

      else

       i = LEN_TRIM(name1)-1; AF = name1(1:i)//ctype1
       Open(inb1,file=AF,status='OLD',action='READ')
       nstate1 = Idef_st(inb1)                                  

      end if

      if(ctype2.eq.'b') then

       AF='bound.'//ALS2
       Open(inb2,file=AF,status='OLD',action='READ')
       rewind(inb2)

       read(inb2,*) ns2,nch2,ncp2,nhm2,nstate2

       if(ns2 .ne.kns ) Stop ' ns2 <> ns '
       if(nch2.ne.kch2) Stop ' nch2 --> ?'
       if(ncp2.ne.kcp2) Stop ' ncp2 --> ?'
       if(nhm2.ne.kdm2) Stop ' nhm2 --> ? '

      elseif(ctype2.eq.'c') then
	   
       inb2=in2; nstate2=1

      else

       i = LEN_TRIM(name2)-1; AF = name2(1:i)//ctype2
       Open(inb2,file=AF,status='OLD',action='READ')
       nstate2 = Idef_st(inb2)                                  

      end if

      write(pri,*)
      write(pri,'(a,i5,a)') 'nstate1 =',nstate1, &
                           ' - number of states in 1-st set'
      write(pri,'(a,i5,a)') 'nstate2 =',nstate2, &
                           ' - number of states in 2-nd set'

      if(mstate1.le.0.or.mstate1.gt.nstate1) mstate1=nstate1
      if(mstate2.le.0.or.mstate2.gt.nstate2) mstate2=nstate2

!----------------------------------------------------------------------
! ... generation of f-values:

      Open(nur,file=AF_res,position='APPEND')

      if(jmode.eq.0)   Call LS_case
      if(jmode.ne.0)   Call J_case

      End Subroutine Gen_zf
