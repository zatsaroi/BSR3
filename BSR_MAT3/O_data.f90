!======================================================================
      Subroutine O_data
!======================================================================
!     processing of total overlaps
!----------------------------------------------------------------------
!     possible structures:
!
!     c,  (ic,jc)                        -  bound-bound
!     c < i | . >, ic                    -  bound-channel
!     c < i | . > < j | . >              -  channel-channel
!
!     where .  denotes bound orbital, i,j - channels.
!
!     < i | j > elements are ignored becaUse we assume that target
!     states are orthogonal. These elements are Used to check this 
!     orthogonality (see array 'to' in bsr_matrix).
!
!----------------------------------------------------------------------
      Use bsr_mat
      Use cmdata  
      Use spline_param; Use spline_orbitals; Use spline_galerkin

      Implicit none
      Integer :: i,j, io,jo, ich,jch, ic,jc, i1,i2, j1,j2
      Real(8) :: C, v(ns),w(ns)

!----------------------------------------------------------------------
      if(itype.ne.1) Stop ' O-data: itype <> 1 '                                  

      Do j=1,ncdata;  i=IPT(j)

       ic=-k3(i); jc=-k4(i); C=cdata(i)

       if(abs(C).lt.Eps_c) Cycle

       if(ic.gt.0.and.jc.gt.0) then                !   (ic,jc)

        Call Update_HB(ic,jc,C)

       elseif(ic.lt.0.and.jc.gt.0) then            !   <i|.> ic

        io = -ic;  i1=io/ibo; i2=mod(io,ibo)

        Call GET_V(i1,i2,v)  
        ich=iech(i1); ic = jc
        Call UPDATE_HV(ich,ic,ns,v,C)

       elseif(ic.lt.0.and.jc.lt.0) then            !   <i|.> <.|j>

        io = -ic;  i1=io/ibo; i2=mod(io,ibo)
        Call GET_V(i1,i2,v)
        jo = -jc;  j1=jo/ibo; j2=mod(jo,ibo)
        Call GET_V(j1,j2,w)
        ich = iech(i1); jch = iech(j1)
        v = v * C

        Call UPDATE_HW(ich,jch,ns,v,w)         

       elseif(ic.lt.0.and.jc.eq.0) then            !    <i|j>         

        io = -ic;  i1=io/ibo; i2=mod(io,ibo)

        ich = iech(i1); jch = iech(i2)

        Call Target_h(ich,jch,0.d0,C)
        if(iitar.gt.1.and.ich.ne.jch) &
         Call UPDATE_HL(ich,jch,ns,ks,sb,C)
       else

        Stop ' O_data: unknown structure '

       end if

      End do

      End Subroutine O_DATA
