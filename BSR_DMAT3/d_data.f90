!======================================================================
      Subroutine  D_data 
!======================================================================
!     processing of d-integrals in the module 'cmdata'
!----------------------------------------------------------------------
!
!    we have following 10 different structures for radial integrals:
!
! 1  d( . .)  ic, jc               -  bound-bound  

! 2  d( i .)  jc                   -  bound-channel
! 3  d( . j)  ic                    
! 4  d( . .) < i | . > jc           
! 5  d( . .) < . | j > ic           

! 6  d( i j)                       -  channel-channel
! 7  d( i .) < . | j >              
! 8  d( . j) < i | . >               
! 9  d( . .) < i | . > < . | j >    
!10  d( . .) < i | j >              
!
! where . denotes bound orbital, i,j - channels, ic,jc - configurations
!
!----------------------------------------------------------------------
      Use spline_param,    only: ns,ks
      Use spline_orbitals, only: iech,qbs
      Use bsr_dmat,        only: ibo,bbbs,pri
      Use cmdata

      Implicit none
      Integer :: i,j, i1,i2, ich,jch, ic,jc, ib,jb, io,jo
      Real(8) :: v(ns),w(ns), v1(ns),v2(ns), w1(ns),w2(ns), &
	                xl(ns,ks+ks-1),xv(ns,ks+ks-1)

      Select Case(itype)

      Case(1)                            !  d( . .)  ic, jc                             

       Do j=1,ncdata;  i=IPT(j)
        Call Update_DB(k1(i),k2(i),CLDATA(i),CVDATA(i))
       End do

      Case(2)                            !  d( i .)  jc                             

       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); ich=iech(i1); jch=0; ic=0; jc=k3(i)
        Call Get_dvr(i1,i2,cldata(i),cvdata(i),v,w)
        Call UPDATE_DV(ich,jch,ic,jc,ns,v,w)
       End do

      Case(3)                             !  d( . j)  ic                              

       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); ich=0; jch=iech(i2); ic=k3(i); jc=0
        Call Get_dvl(i1,i2,cldata(i),cvdata(i),v,w)
        Call UPDATE_DV(ich,jch,ic,jc,ns,v,w)
       End do

      Case(4)                             !  d( . .) < i | . > jc                            

       Do j=1,ncdata;  i=IPT(j)
        io=k1(i); ich=iech(io/ibo); ib=mod(io,ibo); jch=0 
        if(ich.eq.0) Stop 'D_data: itype=4, but ich=0'
        ic=0; jc=k2(i)
        v = cldata(i) * qbs(:,ib); w = cvdata(i) * qbs(:,ib)
        Call UPDATE_DV(ich,jch,ic,jc,ns,v,w)
       End do

      Case(5)                             !  d( . .) < . | j > ic                            

       Do j=1,ncdata;  i=IPT(j)
        io=k1(i); ich=0; jch=iech(mod(io,ibo)); jb=io/ibo
        if(jch.eq.0) Stop 'D_data: itype=5, but jch=0'
        ic=k2(i); jc=0
        v = cldata(i) * qbs(:,jb); w = cvdata(i) * qbs(:,jb)
        Call UPDATE_DV(ich,jch,ic,jc,ns,v,w)
       End do

      Case(6)                             !  d( i j)                            
  
       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); ich=iech(i1); jch=iech(i2)
        if(ich.eq.0) Stop 'D_data: itype=6, but ich=0'
        if(jch.eq.0) Stop 'D_data: itype=6, but jch=0'
        Call Get_dm(i1,i2,cldata(i),cvdata(i),xl,xv)
        Call UPDATE_DX(ich,jch,ns,ks,xl,xv)
       End do

      Case(7)                              !  d( i .) < . | j >     

       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); ich=iech(i1) 
        if(ich.le.0) Stop 'D_data: itype=7, but ich=0'
        io=k3(i); jch=iech(mod(io,ibo)); jb=io/ibo
        if(jch.le.0) Stop 'D_data: itype=7, but jch=0'
        w1 = qbs(:,jb) 
	       w2 = qbs(:,jb) 
        Call Get_dvr(i1,i2,cldata(i),cvdata(i),v1,v2)
        Call UPDATE_DW(ich,jch,ns,v1,w1,v2,w2)         
       End do

      Case(8)                              !  d( . j) < i | . >     

       Do j=1,ncdata;  i=IPT(j)
        i1=k1(i); i2=k2(i); jch=iech(i2) 
        if(jch.le.0) Stop 'D_data: itype=4, but jch=0'
        io=k3(i); ich=iech(io/ibo); ib=mod(io,ibo) 
        if(ich.le.0) Stop 'D_data: itype=4, but ich=0'
        v1 = qbs(:,ib) 
	       v2 = qbs(:,ib) 
        Call Get_dvl(i1,i2,cldata(i),cvdata(i),w1,w2)
        Call UPDATE_DW(ich,jch,ns,v1,w1,v2,w2)         
       End do

      Case(9)                              !  d( . .) < i | . > < . | j >                       

       Do j=1,ncdata;  i=IPT(j)

        io=k1(i); ich=iech(io/ibo); ib=mod(io,ibo)
        jo=k2(i); jch=iech(mod(jo,ibo)); jb=jo/ibo

        if(ich.eq.0) Stop 'D_data: itype=9, but ich=0'
        if(jch.eq.0) Stop 'D_data: itype=9, but jch=0'
        v1 = cldata(i)*qbs(:,ib);  w1 = qbs(:,jb)
        v2 = cvdata(i)*qbs(:,ib);  w2 = qbs(:,jb)
        Call UPDATE_DW(ich,jch,ns,v1,w1,v2,w2)         

       End do

      Case(10)                           !   d( . .) < i | j >                            

       Do j=1,ncdata;  i=IPT(j)
        io=k1(i); ich=iech(io/ibo); jch=iech(mod(io,ibo))
        if(ich.eq.0) Stop 'D_data: itype=10, but ich=0'
        if(jch.eq.0) Stop 'D_data: itype=10, but jch=0'
        xl = cldata(i) * bbbs; xv = cvdata(i) * bbbs
        Call UPDATE_DX(ich,jch,ns,ks,xl,xv)
       End do

      Case Default

       Stop ' D_data: uknown itype '

      End Select
      
      End Subroutine D_data



