!======================================================================
      Subroutine I_data(jtype,jpol)
!======================================================================
!     processing of I-integrals from the module 'cmdata'
!----------------------------------------------------------------------
!     we have following different structures:
!
! 1 1.0  Rk( . . . .)  ic, jc               -  bound-bound  
! 2 1.1  Rk( . . . .) < i | . > ic          -  bound-channel
! 3 1.2  Rk( . . . .) < i | . > < j | . >   -  channel-channel
! 4 1.3  Rk( . . . .) < i | j >             -  target structure
!
! 5 2.0  Rk( i . . .) < j | . >             -  channel-channel due to
! 6 3.0  Rk( . i . .) < j | . >                overlaps
! 7 4.0  Rk( . . i .) < j | . >
! 8 5.0  Rk( . . . i) < j | . >
!
! 9 2.1  Rk( i . . .)  ic                   -  bound-channel
!10 3.1  Rk( . i . .)  ic
!11 4.1  Rk( . . i .)  ic
!12 5.1  Rk( . . . i)  ic
!
!13 6.0  Rk( i . j .)                       -  direct channel-channel
!14 7.0  Rk( . i . j)
!15 8.0  Rk( i . . j)                       -  exchange channel-channel
!16 9.0  Rk( . i j .)
!
!     where .  denotes bound orbital, i,j - channels.
!
!     Rk * < i | j >  elements with i<>j are ignored becaUse
!     we assume that target states diagonalize Hamiltonian
!     These elements are included after in B(.,.) * Etarget(i)
!     where B(.,.) is B-spline overlap matrix
!     These elements are also Used for control calculation of 
!     interaction matrix between target states (Target_h).
!
!     sym_d -  symmetry for convolution
!     sym_r -  symmetry for result
!
!     jcase = 1 -->  INT( . p1; . p2)
!     jcase = 2 -->  INT( p1 .;p2 . )
!     jcase = 3 -->  INT( . p1;p2 . )
!     jcase = 4 -->  INT( p1 .; . p2)
!----------------------------------------------------------------------
      Use cmdata; Use bsr_mat; Use target
      Use spline_param; Use spline_galerkin
      Use spline_orbitals, p => pbs

      Implicit none
      Character(1) :: sym_d,sym_r,sym_v, sym_i,sym_j
      Integer :: jtype, jpol,  int1,int2,int3,int4
      Integer :: i,j, i1,i2, j1,j2, ii,jj, iii,jjj, ich,jch, ic,jc
      Integer :: jcase, io,jo
      Real(8) :: C,S,CA
      Real(8) :: v(ns), w(ns)
      Real(8) :: xx(ns,ns), dd(ns,ns)
      Real(8), external :: Sum_AmB, QUADR

      if(ncdata.eq.0) Return

      dd=0.d0; xx=0.d0; v=0.d0; w=0.d0

!----------------------------------------------------------------------
! ... prepare B_spline representation and define symmetries: 

      Select case(icase)
       Case( 3); Call MTK_cell(jpol); sym_i='n'; sym_j='n'    !  Tk
       Case( 4); Call MMK_cell(jpol); sym_i='s'; sym_j='s'    !  Mk
       Case( 5); Call MRK_cell(jpol); sym_i='s'; sym_j='s'    !  Rk
       Case( 8); Call MNK_cell(jpol); sym_i='s'; sym_j='s'    !  Nk
       Case( 9); Call MVK_cell(jpol); sym_i='n'; sym_j='s'    !  Vk
       Case(10); Call MNK_cell(jpol); sym_i='s'; sym_j='s'    !  Nk
       Case Default; Stop 'unknown case in I_data'
      End Select

      jcase = IJCASE(jtype)
      Select case(jcase)
       Case(1);      sym_d=sym_j; sym_r=sym_i
       Case(2);      sym_d=sym_i; sym_r=sym_j 
       Case Default; sym_d='x';   sym_r='x'
      End Select

!----------------------------------------------------------------------
! ... select the given structure:

      Select Case(jtype)
!----------------------------------------------------------------------
!                                                       Rk( . . ; . . )      
      Case(1)                                  

        iii=0; jjj=0
        Do j=1,ncdata;  i=IPT(j)

         if(K1(i).ne.iii) then
          i1 = K1(i)/ibi; i2 = mod(K1(i),ibi); int1=i1; int3=i2
          Call Density(ns,ks,dd,p(1,i1),p(1,i2),sym_d)
          Call Convol(ns,ks,xx,dd,jcase,sym_i,sym_j)
         end if

         if(K1(i).ne.iii.or.K2(i).ne.jjj) then
          j1 = K2(i)/ibi; j2 = mod(K2(i),ibi); int2=j1; int4=j2
          Call Density(ns,ks,dd,p(1,j1),p(1,j2),sym_r)
          S = SUM_AmB(ns,ks,xx,dd,sym_r)
          if(nud.gt.0) &
          Call pri_int(nud,icase,jpol,int1,int2,int3,int4,S)
         end if

         iii=K1(i);  jjj=K2(i);   C=S*cdata(i)

         ic=-k3(i);  jc=-k4(i);   io=-ic; jo=-jc

         if(ic.gt.0.and.jc.gt.0) then                ! R(..) (ic,jc)

          Call Update_HB(ic,jc,C)

         elseif(ic.lt.0.and.jc.gt.0) then            ! R(..) <i|.> ic

          i1 = io/ibo; i2 = mod(io,ibo)
          Call GET_V(i1,i2,v)
          ich=iech(i1)
          Call UPDATE_HV(ich,jc,ns,v,C)
  
         elseif(ic.lt.0.and.jc.lt.0) then            ! R(..) <i|.> <.|j>

          i1=io/ibo; i2=mod(io,ibo)
          Call GET_V(i1,i2,v)
          j1=jo/ibo; j2=mod(jo,ibo)
          Call GET_V(j1,j2,w)
          ich = iech(i1); jch = iech(j1)
          v = v * C
          Call UPDATE_HW(ich,jch,ns,v,w)         

         elseif(ic.lt.0.and.jc.eq.0) then            ! R(..) <i|j>  

          i1=io/ibo; i2=mod(io,ibo)
          ich = iech(i1); jch = iech(i2)

          if(coupling.eq.'LS'.and.icase.gt.7) then
           Call UPDATE_HL(ich,jch,ns,ks,sb,C)
          else
           Call Target_h(ich,jch,C,0.d0)
           if(iitar.gt.0.and.ich.ne.jch) &
            Call UPDATE_HL(ich,jch,ns,ks,sb,C)
          end if

         else

          Stop ' I_data: unknown structure for itype=1'

         end if

        End do

!----------------------------------------------------------------------
!                                                       R ( i . ; . . ) 
      Case(2,3,4,5)                               
      
        jjj=0; iii=0
        Do j=1,ncdata; i=IPT(j); jj=k1(i); ii=k2(i); ich=k3(i); io=k4(i)

         if(jj.ne.jjj.or.ii.ne.iii) then

          if(jj.ne.jjj) then
           j1=jj/ibi; j2=mod(jj,ibi)
           Call Density(ns,ks,dd,p(1,j1),p(1,j2),sym_d)
           Call Convol(ns,ks,xx,dd,jcase,sym_i,sym_j)
           jjj = jj
          end if

          sym_v = 'l';  if(icase.gt.2) sym_v = 'r'
          CALL BAV(ns,ks,xx,p(1,ii),v,sym_r,sym_v)
          iii = ii 

         end if

        if(io.gt.0) then 
          j1 = io/ibo; j2 = mod(io,ibo); jch=iech(j1)
          Call GET_V(j1,j2,w)
          w = cdata(i) * w
          Call UPDATE_HW(ich,jch,ns,v,w)
        elseif(io.lt.0) then
          ic  = -io
          Call UPDATE_HV(ich,ic,ns,v,cdata(i))
        else
          Stop ' I_data: uknown structure for itype 2-5 '
        end if

        End do     

!----------------------------------------------------------------------
!                                                      RK ( i . ; j . )
      Case(6,7,8,9)                         

       xx=0.d0;  CA=0.d0

       Do j=1,ncdata;  i=IPT(j); i1=k2(i); j1=k3(i)

        Call Density(ns,ks,dd,p(1,i1),p(1,j1),sym_d)
        xx = xx + cdata(i)*dd
        CA = CA + cdata(i)*QUADR(i1,j1,jpol)
       
        if(j.lt.ncdata) then
         if(k1(i).eq.k1(IPT(j+1))) Cycle
        end if 

        Call Convol(ns,ks,dd,xx,jcase,sym_i,sym_j)

        ich=k1(i)/ibi; jch=mod(k1(i),ibi)

        Call UPDATE_HX(ich,jch,ns,ks,dd,sym_r)

        if(icase.eq.5.and.jcase.le.2) Call UPDATE_CF(jpol,ich,jch,CA)

        xx=0.d0; CA=0.d0

       End do

      Case Default

       Stop ' I_data: unknown itype '

      End Select

      End Subroutine  I_data

