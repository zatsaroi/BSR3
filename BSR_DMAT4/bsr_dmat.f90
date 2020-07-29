!======================================================================
!     PROGRAM       B S R _ D M A T                        version. 3
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny (pilgrim)
!     email:        oleg_zoi@yahoo.com
!
!======================================================================
!     Provides different multipole-transition matrixes or 
!     oscillator strengths and transition probabilities   
!======================================================================
!     Multipole index is taken from mult_bnk !!!
!======================================================================
!
!     INPUT ARGUMENTS:
!
!     first four mandatory position parameters:
!
!     1.  name1.c or cfg.nnn - c-file for initial state 
!     2.  name2.c or cfg.nnn - c-file for final state 
!     3.  structure mode for initial state  (c,l,j,b)
!     4.  structure mode for final state    (c,l,j,b,p,q,d)
! 
!     additonal optional key paramters (as par=value):
!
!     istate1 [0] - index of state from 1st set of solutions 
!                   to be considered, 0 means all solutions
!     istate2 [0] - index of state from 2nd set of solutions 
!                   to be considered, 0 means all solutions
!     mstate1 [0] - restriction on the number of 1st-set solutions   
!                   to be considered, 0 means all solutions
!     mstate2 [0] - restriction on the number of 2st-set solutions   
!                   to be considered, 0 means all solutions
!     gf      [f] - output f- of gf-values
!     ialfa   [0] - additional output radial matrix elements and
!                   polarizabilities (in zf_res) 
!----------------------------------------------------------------------
!     Some name conventions are taken from MCHF package of C.F.F.,
!     1995 CPC version.
!
!     ctype = c  -  expantion coefficients are read from c-file
!     ctype = l  -  expantion coefficients are read from l-file 
!     ctype = j  -  expantion coefficients are read from j-file 
!
!     ctype = b  -  expantion coefficients are read from bound.nnn file
!                   after BSR bound-type calculations
!
!     ctype = p  -  2nd set of solutions are supposed to be R-matrix
!                   solutions in rsol.nnn file after BSR calculations 
!     ctype = q  -  cfg.nnn is expected for 2nd set but no expansion
!                   coefficients are required 
!     ctype = d  -  cfg.nnn and cfg.mmm are expected for both 1st and
!                   2nd set but no expansion coefficients are required  
!----------------------------------------------------------------------
!
!     INPUT FILES:
!
!     mult_bnk  - data bank of ang. coefficients for multipole operator
!
!     knot.dat  - B-spline information
!
!     name1.c   - c-file for initial state 
!     name1.bsw - initial state wavefunctions in B-spline representation
!
!     name2.c   - c-file for final state
!     name2.bsw - initial state wavefunctions in B-spline representation
!
!     cfg.nnn   - One (or both) c-files can be replaced by cfg.nnn files,
!                 where nnn is a partial wave index (for ctype=b).
!                 The close-coupling expansion is supposed for this case,
!                 that requires also following files (instead name.bsw):
!
!     target        - information about close-coupling expansion
!                     (target states, channels and perturbers)
!     target.bsw    - target wavefunctions in B-spline representation
!     pert_nnn.bsw  - if any (according to target information)
!
!     bound.nnn  - expansions coefficients if ctype=b  
!     rsol.nnn   - R-matrix solutions if ctype2=p
!-----------------------------------------------------------------------
!
!     OUTPUT FILES:
!
!     bsr_dmat.log - running information
!
!     zf_res    - oscilator strengths ( ctype2 = c,l,j,b )
!
!     d.nnn     - dipole matrix for R-matrix solutions   (ctype2 = p)
!
!     dv.nnn    - dipole vector for polarizability calc. (ctype2 = q)
!
!     d.nnn_mmm - clear dipole matrix without convolution with any 
!                 initial or final solution expansions   (ctype2 = d) 
!-----------------------------------------------------------------------
!     examples for calling:
!
!     bsr_dmat 1.c 2.c c c 
!     bsr_dmat 1.c 2.c j j  gf=g  ialfa=1
!     bsr_dmat 1.c cfg.001 c b  mstate2=5 
!     bsr_dmat 1.c cfg.002 c p  
!     bsr_dmat cfg.001 cfg.002 b b istate1=2 istate2=5  
!     bsr_dmat 1.c cfg.002 c p  
!     bsr_dmat 1.c cfg.002 c q  
!     bsr_dmat cfg.001.c cfg.002 b d  
!-----------------------------------------------------------------------
!     RESTRICTION:  in b->j or j->b calculations, j-file may contain
!                   data only for one value of J.
!----------------------------------------------------------------------
      Use bsr_dmat; Use dmatrix 
      Use conf_LS;  Use orb_LS; Use term_LS
      Use channels; Use target
      Use spline_param; Use spline_atomic; Use spline_orbitals 

      Implicit none
      Real(8) :: t1,t2,AWT
      Integer :: i,j,n,l,k
      Real(8), external :: DJ_fact,DJM_fact
      Integer, external :: Jfind_bsorb, IBORT

      Call CPU_time(t1)

      Call inf_bsr_dmat

      Open(pri,file=AF_log)

      write(*,'(a/a)') 'B S R _ D M A T',&
                       '***************'
!----------------------------------------------------------------------
! ... read arguments from command line and open basic files:

      Call Read_arg

      Call Check_file(name1);  Open(in1,file=name1) 
      if(name1.eq.name2) then
       in2=in1
      else
       Call Check_file(name2);  Open(in2,file=name2) 
      end if

      Call Check_file(AF_bnk)
      Open(nub,file=AF_bnk,form='UNFORMATTED') 
      read(nub) ktype,kpol

      write(*,*)
      write(*,'(a,a1,i1)') 'multipole mode:  ',ktype,kpol

! ... check mult_bnk: 

      Call R_closed(in1)
      kclosd=nclosd    
      Call R_closed(in2)
      if(kclosd.ne.nclosd) Stop 'different cores ?'

      Call R_conf(in1,in2,nub)

! ... define terms:

      Call R_term(in1)
      Call Shift_term1
      Call Def_term(in1,ILT1,IST1,IPT1)
      jot1 = 0; if(IST1.eq.0) jot1=ILT1 + 1
      Call R_term(in2)
      Call Shift_term2
      Call Def_term(in2,ILT2,IST2,IPT2)
      jot2 = 0; if(IST2.eq.0) jot2=ILT2 + 1
      Allocate( DJ(nterm1,nterm2), DJM(nterm1,nterm2) )

!----------------------------------------------------------------------
! ... read channel information from target for BSR calculations:

      if(ilsp1.gt.0.or.ilsp2.gt.0) then
       Call Check_file(AF_tar); Open(nut,file=AF_tar)
       Call R_target (nut)
       Call R_channels(nut)
       Close(nut)
      end if

      nch1=0; ncp1=ncfg1; inb1=in1; npert1=ncp1
      if(ilsp1.gt.0) then
       nch1   = nch  (ilsp1)
       ncp1   = ncp  (ilsp1)
       npert1 = npert(ilsp1)
       ipert1 = ipert(ilsp1)
       jot1   = jpar (ilsp1)
      end if

      nch2=0; ncp2=ncfg2; inb2=in2; npert2=ncp2
      if(ilsp2.gt.0) then
       nch2   = nch  (ilsp2)
       ncp2   = ncp  (ilsp2)
       npert2 = npert(ilsp2)
       ipert2 = ipert(ilsp2)
       jot2   = jpar (ilsp2)
      end if

      if(jot1.eq.0.and.jot2.gt.0 .or. jot1.gt.0.and.jot2.eq.0)  &
       Stop 'inconsistent J-dependence for initial and final states'

!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline: 
    
      CALL define_grid(z);  CALL define_spline

!----------------------------------------------------------------------
! ... read B-spline expantions for bound orbitals:
   
      Call Allocate_bsorb(nwf)
      nbf = nwf
      Do i = 1,nbf
       ebs(i)=ELF(i); nbs(i)=NEF(i); lbs(i)=LEF(i); kbs(i)=KEF(i)
       iech(i) = 0; mbs(i) = 0
      End do

! ... radial functions for initial state:

      if(ilsp1.eq.0) then

       i=INDEX(name1,'.',BACK=.TRUE.); AF=name1(1:i)//'bsw'
       Call Check_file(AF)       
       Open(nuw,file=AF,form='UNFORMATTED')
       Call Read_bsw(nuw,kset1);  Close(nuw)

      else

       Call Check_file(AF_bsw)       
       Open(nuw,file=AF_bsw,form='UNFORMATTED')
       Call Read_bsw(nuw,kset1);  Close(nuw)

       if(nwp(ilsp1).gt.0) then
        i = LEN_TRIM(AFP(ilsp1)); AF = AFP(ilsp1); AF = AF(1:i)//'.bsw'
        Call Check_file(AF)       
        Open(nuw,file=AF,form='UNFORMATTED')
        Call Read_bsw(nuw,kset1);  Close(nuw)
       end if

       Do i = 1,nch(ilsp1)
        Call EL4_nlk(ELC(ilsp1,i),n,l,k);  k = k + kset1
        j = Jfind_bsorb(n,l,k); iech(j) = i
       End do

      end if

! ... radial functions for final state:

      if(ilsp2.eq.0) then

       i=INDEX(name2,'.',BACK=.TRUE.); AF=name2(1:i)//'bsw'
       Call Check_file(AF)       
       Open(nuw,file=AF,form='UNFORMATTED')
       Call Read_bsw(nuw,kset2);  Close(nuw)

      else

       Call Check_file(AF_bsw)       
       Open(nuw,file=AF_bsw,form='UNFORMATTED')
       Call Read_bsw(nuw,kset2);  Close(nuw)

       if(nwp(ilsp2).gt.0) then
        i = LEN_TRIM(AFP(ilsp2)); AF = AFP(ilsp2); AF = AF(1:i)//'.bsw'
        Call Check_file(AF)       
        Open(nuw,file=AF,form='UNFORMATTED')
        Call Read_bsw(nuw,kset2);  Close(nuw)
       end if

       Do i = 1,nch(ilsp2)
        Call EL4_nlk(ELC(ilsp2,i),n,l,k);  k = k + kset2
      	j = Jfind_bsorb(n,l,k); iech(j) = i
       End do

      end if

! ... check the correspondence between c- and w-files: 

      j = 0
      Do i = 1,nwf
       ief(i) = iech(i)
       if(iech(i).ne.0) Cycle
       if(mbs(i).eq.0) then
        write(pri,'(a,a)') ' Absent expansion for w.f. ',ELF(i)
        j = j + 1
       end if
      End do
      if(j.gt.0) Stop 'no correspondence between c- and bsw- files'
     

! ... the < p | p > values 

      Call Alloc_orb_overlaps(nbf,lbs,iech,0)

!----------------------------------------------------------------------
! ... expansion coefficients for initial and final terms:

      Allocate(C1(ncfg1))  
      Allocate(C2(ncfg2))

      C1 = 1.d0; if(ilsp1.gt.0) C1=WC(1:ncfg1)
      C2 = 1.d0; if(ilsp2.gt.0) C2=WC(ncfg1+1:ncfg1+ncfg2)

!----------------------------------------------------------------------
! ... Rydberg constants:

      AWT=0.d0; Call Conv_au (Z,AWT,au_cm,au_eV,pri)

!----------------------------------------------------------------------
! ... generate the LSJ-factors in case of fine-structure calculations:

      if(jot1+jot2.gt.0) jmode=1      
      if(jmode.eq.1.and.(ilsp1.gt.0.or.ilsp2.gt.0)) jmode=2          

! ... jmode = 2  =>  DJ factor will be incorporated in ang. coefficients !

      if(jmode.eq.2) then

       Do i=1,nterm1; Do j=1,nterm2
        DJ(i,j) = &
        DJ_fact(ILterm1(i),ILterm2(j),ISterm1(i),ISterm2(j),JOT1,JOT2,kpol)
       End do; End do

       if(ktype.eq.'M') then
        Do i=1,nterm1; Do j=1,nterm2
         DJM(i,j)= &
         DJM_fact(ILterm1(i),ILterm2(j),ISterm1(i),ISterm2(j),JOT1,JOT2,kpol)
        End do; End do
       else
	       DJM = DJ
       end if

      end if

      if(ktype.eq.'M'.and.jmode.eq.0) &
       Stop 'Magnetic transitions are not defined for LS states'

      if(ktype.eq.'M') ialfa = 0

!----------------------------------------------------------------------
! ... print the information:

      write(pri,*)
      write(pri,'(a,a)' ) 'initial state - ',name1
      write(pri,'(a,i6)') 'ncfg1  =',ncfg1
      write(pri,'(a,i6)') 'nterm1 =',nterm1
      if(jot1.gt.0) write(pri,'(a,i6)') 'jot1   =',jot1
      write(pri,*)
      write(pri,'(a,a,i3,a,i3)') 'orbitals in first set: ', &
	        '   nwf1 = ',nwf1,'   index shift = ',kset1
      write(pri,*)
      if(debug.gt.0) write(pri,'(15a5)') ELF(1:nwf1)

      write(pri,*)
      write(pri,'(a,a)' ) 'final state   - ',name2
      write(pri,'(a,i6)') 'ncfg2  =',ncfg2
      write(pri,'(a,i6)') 'nterm2 =',nterm2
      if(jot2.gt.0) write(pri,'(a,i6)') 'jot2   =',jot2
      write(pri,'(a,a,i3,a,i3)') 'orbitals in second set:', &
	        '   nwf2 = ',nwf2,'   index shift = ',kset2
      write(pri,*)
      if(debug.gt.0) write(pri,'(15a5)') ELF(nwf1+1:nwf)

!----------------------------------------------------------------------
! ... generation D-matrix:

      Call allocate_matrix(ns,ks,nch1,nch2,npert1,npert2,ncfg1,ncfg2)

      Call D_matr

!----------------------------------------------------------------------
! ... outputs of results:

      if(dd.eq.'y') then
        write(AF,'(a,i3.3,a,i3.3)') 'dd.',ilsp1,'_',ilsp2      
        open(nudd,file=AF)
      end if

      Select case(jout)
       Case(0);                 Call Gen_zf
       Case(1);                 Call D_out
       Case(2);                 Call Dvec_out
       Case(3);                 Call DD_out
       Case default;            Call DD_out
      End Select

      Call CPU_time(t2)
      write(pri,'(/a,f10.2,a)') 'BSR_DMAT:',(t2-t1)/60,' min '
      write(*  ,'( a,f10.2,a)') 'BSR_DMAT:',(t2-t1)/60,' min '

      End ! program BSR_DMAT3
