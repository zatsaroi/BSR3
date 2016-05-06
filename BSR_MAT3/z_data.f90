!======================================================================
      Subroutine Z_data (jtype,l) 
!======================================================================
!     processing the Z-integrals from the module 'cmdata'
!----------------------------------------------------------------------
!     we have following different structures:
!
!     Z( . .)                        -  bound-bound
!     Z( i .)                        -  bound-channel
!     Z( . .) < i | . >              -  bound-channel
!     Z( i i)                        -  channel-channel
!     Z( i .) < j | . >              -  channel-channel
!     Z( . .) < i | . > < j | . >    -  channel-channel
!
!     where .  denotes bound orbital, i,j - channels.
!----------------------------------------------------------------------
      Use Z_core; Use bsr_mat
      Use cmdata; Use target;  Use channel
      Use spline_param; Use spline_orbitals; Use spline_galerkin
     
      Implicit none
      Integer, intent(in) :: jtype,l
      Integer :: i,j, i1,i2, j1,j2, ich,jch, ic,jc, io,jo
      Real(8) :: c, v(ns), w(ns)

      zl_core(:,:) = zl_full(:,:,l)                

      Select Case(jtype)
!----------------------------------------------------------------------
!                                                            Z ( . . ):
      Case(1)                                

       Do j=1,ncdata; i=IPT(j)
  
        if(abs(cdata(i)).lt.Eps_C) Cycle

        i1=k1(i)/ibi; i2=mod(k1(i),ibi);  C = Z_int(i1,i2)*cdata(i)

        ic=-k3(i); jc=-k4(i)
       
         if(ic.gt.0.and.jc.gt.0) then                ! Z(..) (ic,jc)

          Call Update_HB(ic,jc,C)

         elseif(ic.lt.0.and.jc.gt.0) then            ! Z(..) <i|.> ic

          io = -ic;  i1=io/ibo; i2=mod(io,ibo)
          Call GET_V(i1,i2,v)  
          ich=iech(i1); ic = jc
          Call UPDATE_HV(ich,ic,ns,v,C)

         elseif(ic.lt.0.and.jc.lt.0) then            ! Z(..) <i|.> <.|j>

          io = -ic;  i1=io/ibo; i2=mod(io,ibo)
          Call GET_V(i1,i2,v)
          jo = -jc;  j1=jo/ibo; j2=mod(jo,ibo)
          Call GET_V(j1,j2,w)
          ich = iech(i1); jch = iech(j1)
          v = v * C
          Call UPDATE_HW(ich,jch,ns,v,w)         

         elseif(ic.lt.0.and.jc.eq.0) then            ! Z(..) <i|j>                            !   < i | j >

          io = -ic;  i1=io/ibo; i2=mod(io,ibo)
          ich = iech(i1); jch = iech(i2)

          if(coupling.eq.'LS') then
		         Call UPDATE_HL(ich,jch,ns,ks,sb,C)
		        elseif(iitar.gt.0.and.ich.ne.jch) then
		         Call UPDATE_HL(ich,jch,ns,ks,sb,C)
          else
           Call Target_h(ich,jch,C,0.d0)
		        end if

         else

          Stop ' Z_data: unknown structure '

         end if

       End do

!----------------------------------------------------------------------
!                                                            Z ( . i ):
      Case(2,3,4,5)                                   

       Do j=1,ncdata; i=IPT(j); C=cdata(i)  

         if(abs(C).lt.Eps_C) Cycle

         v=Z_vec(:,k2(i)); ich=k3(i); io=k4(i) 

         if(io.gt.0) then                           !  Z(.i) <.|j>

          j1=io/ibo; j2=mod(io,ibo); jch=iech(j1)

          Call GET_V(j1,j2,w)
          w=C*w
          Call UPDATE_HW(ich,jch,ns,v,w)

         elseif(io.lt.0) then                       !  Z(.i)  ic 

          ic=-io
          Call UPDATE_HV(ich,ic,ns,v,C)

         else
          Stop ' Z_data: uknown structure for itype 2-5 '
         end if

       End do     
!----------------------------------------------------------------------
!                                                            Z ( i j ):
      Case(6,7,8,9)                         

       Do j=1,ncdata; i=IPT(j); C=cdata(i)
        if(abs(C).lt.Eps_C) Cycle
        ich=k1(i)/ibi; jch=mod(k1(i),ibi)
        xx=C*zl_core; Call UPDATE_HX(ich,jch,ns,ks,xx,'x')
       End do

      Case Default

       Stop ' Z_data: uknown itype '

      End Select

      End Subroutine  Z_data

