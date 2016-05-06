!======================================================================
      Subroutine  L_data(jtype,l) 
!======================================================================
!     processing of L-integrals in the module cmdata
!----------------------------------------------------------------------
!     we have following different structures:
!
!     L( . .)                        -  bound-bound
!     L( i .)                        -  bound-channel
!     L( . .) < i | . >              -  bound-channel
!     L( i j)                        -  channel-channel
!     L( i .) < j | . >              -  channel-channel
!     L( . .) < i | . > < j | . >    -  channel-channel
!
!     where .  denotes bound orbital, i,j - channels.
!
!     L( . .) < i | j > elements with i/=j are ignored becaUse
!     we assume that target states diagonalize Hamiltonian
!     L( . .) < i | i > are included after in B(.,.) * Etarget(i)
!     where B(.,.) is B-spline overlap matrix
!     These elements are Used for control calculation of interaction
!     matrix between target states (Target_h).
!     L( i j) with with i/=j are also ignored becaUse we assume
!     that target states are orthpgonal, but we have include
!     interaction with core which is included in L(i,j)
!
!     in int_bnk,  L(a,b) are represented as I(a,a;b,b)
!----------------------------------------------------------------------
      Use spline_param;  Use spline_orbitals; Use spline_galerkin
      Use channel
      Use cmdata
      Use L_core
      
      Implicit none
      Integer, Intent(in) :: jtype,l
      Integer :: i,j, i1,i2, j1,j2, ich,jch, ic,jc, io,jo
      Real(8) :: c, v(ns), w(ns)

      hl_core(:,:) = hl_full(:,:,l) 

      Select Case(jtype)

!----------------------------------------------------------------------
!                                                            L ( . . ):
      Case(1)                                

       Do j=1,ncdata;  i=IPT(j)

        if(abs(cdata(i)).lt.Eps_C) Cycle

        i1=k1(i)/ibi; i2=mod(k1(i),ibi);  C=L_int(i1,i2)*cdata(i)
        ic=-k3(i); jc=-k4(i)
       
         if(ic.gt.0.and.jc.gt.0) then                 !  L(..) (ic,jc)

          Call Update_HB(ic,jc,C)

         elseif(ic.lt.0.and.jc.gt.0) then             !  L(..) <i|.> ic

          io = -ic;  i1=io/ibo; i2=mod(io,ibo)
          Call GET_V(i1,i2,v)  
          ich=iech(i1); ic = jc
          Call UPDATE_HV(ich,ic,ns,v,C)
 
         elseif(ic.lt.0.and.jc.lt.0) then            !  L(..) <i|.> <.|j>

          io = -ic;  i1=io/ibo; i2=mod(io,ibo)
          Call GET_V(i1,i2,v)
          jo = -jc;  j1=jo/ibo; j2=mod(jo,ibo)
          Call GET_V(j1,j2,w)
          ich = iech(i1); jch = iech(j1)
          v = v * C
          Call UPDATE_HW(ich,jch,ns,v,w)         

         elseif(ic.lt.0.and.jc.eq.0) then            !  L(..) <i|j>   
                                                
          io = -ic;  i1=io/ibo; i2=mod(io,ibo)
          ich = iech(i1); jch = iech(i2)

         if(iitar.gt.0.and.ich.ne.jch) &
            Call UPDATE_HL(ich,jch,ns,ks,sb,C)

          Call Target_h(ich,jch,C,0.d0)

         else

          Stop ' L_data: unknown structure '

         end if

        End do

!----------------------------------------------------------------------
!                                                            L ( . i ):
      Case(2,3,4,5)                                   

       Do j=1,ncdata; i=IPT(j); C=cdata(i)  

         if(abs(C).lt.Eps_C) Cycle

         v=L_vec(:,k2(i)); ich=k3(i); io=k4(i) 

         if(io.gt.0) then                               !  L(.i) <.|j>
          j1=io/ibo; j2=mod(io,ibo); jch=iech(j1)
          Call GET_V(j1,j2,w)
          w=C*w
          Call UPDATE_HW(ich,jch,ns,v,w)
         elseif(io.lt.0) then                           !  L(.i)  ic 
          ic=-io
          Call UPDATE_HV(ich,ic,ns,v,C)
         else
          Stop ' L_data: uknown structure for itype 2-5 '
         end if

       End do     

!----------------------------------------------------------------------
!                                                            L ( i j ):
      Case(6,7,8,9)                         

       Do j=1,ncdata;  i=IPT(j); C=cdata(i)

        if(abs(C).lt.Eps_C) Cycle

        ich=k1(i)/ibi; jch=mod(k1(i),ibi)
        if(ich.eq.jch.or.iitar.gt.1) then
          xx = C*hl_core
          Call UPDATE_HX(ich,jch,ns,ks,xx,'x')
        else
          Call Target_h(ich,jch,0.d0,C)
        end if

       End do

      Case Default

       Stop ' L_data: uknown itype '

      End Select

      End Subroutine L_data

