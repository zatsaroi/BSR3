!======================================================================
      Function Idef_itype(i1,i2,i3,i4,io,jo,ic,jc,k1,k2,k3,k4,ibi)
!======================================================================
!
!     we have following 16 different structures for radial integrals:
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
!                                        ibi         ibi
! 1.0 Rk( . . . .)  ic, jc       -  k1=(i1,i3)  k2=(i2,i4)  k3=-ic  k4=-jc 
! 1.1 Rk( . . . .) < i | . > ic  -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4=-ic 
! 1.2 Rk( . . . .) <i|.> <j|.>   -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4= jo
! 1.3 Rk( . . . .) < i | j >     -  k1=(i1,i3)  k2=(i2,i4)  k3= io  k4= 0 
!
! 2.0 Rk( i . . .) < j | . >     -  k1=(i2,i4)  k2=i3  k3=ich  k4=io
! 3.0 Rk( . i . .) < j | . >     -  k1=(i1,i3)  k2=i4  k3=ich  k4=io
! 4.0 Rk( . . i .) < j | . >     -  k1=(i2,i4)  k2=i1  k3=ich  k4=io
! 5.0 Rk( . . . i) < j | . >     -  k1=(i1,i3)  k2=i2  k3=ich  k4=io
!
! 2.1 Rk( i . . .)  ic           -  k1=(i2,i4)  k2=i3  k3=ich  k4=-ic
! 3.1 Rk( . i . .)  ic           -  k1=(i1,i3)  k2=i4  k3=ich  k4=-ic           
! 4.1 Rk( . . i .)  ic           -  k1=(i2,i4)  k2=i1  k3=ich  k4=-ic
! 5.1 Rk( . . . i)  ic           -  k1=(i1,i3)  k2=i2  k3=ich  k4=-ic
!
! 6.0 Rk( i . j .)               -  k1=(ich1,ich2)  k2=i2  k3=i4  k4=0
! 7.0 Rk( . i . j)               -  k1=(ich1,ich2)  k2=i1  k3=i3  k4=0
! 8.0 Rk( i . . j)               -  k1=(ich1,ich2)  k2=i2  k3=i3  k4=0
! 9.0 Rk( . i j .)               -  k1=(ich1,ich2)  k2=i1  k3=i4  k4=0
!
!----------------------------------------------------------------------

      USE bsr_mat, only: pri, debug
      USE orb_LS, ONLY: ief

      Integer, INTENT(in)  :: i1,i2,i3,i4,io,jo,ic,jc,ibi
      Integer, INTENT(out) :: k1,k2,k3,k4

      Idef_itype=0

      if(jc.gt.0) then

       k1=i1*ibi+i3; k2=i2*ibi+i4; k3=-ic; k4=-jc; Idef_itype = 1

      elseif(ic.gt.0) then

       if(io.eq.0) then

        if(ief(i1).gt.0) then
         k1=i2*ibi+i4; k2=i3; k3=ief(i1); k4=-ic; Idef_itype = 2
        elseif(ief(i2).gt.0) then
         k1=i1*ibi+i3; k2=i4; k3=ief(i2); k4=-ic; Idef_itype = 3
        elseif(ief(i3).gt.0) then
         k1=i2*ibi+i4; k2=i1; k3=ief(i3); k4=-ic; Idef_itype = 4
        elseif(ief(i4).gt.0) then
         k1=i1*ibi+i3; k2=i2; k3=ief(i4); k4=-ic; Idef_itype = 5
        end if

       elseif(io.gt.0) then

        k1=i1*ibi+i3; k2=i2*ibi+i4; k3=io; k4=-ic; Idef_itype = 1

       end if

      elseif(io.eq.0) then
       if(ief(i1).gt.0.and.ief(i3).gt.0) then
        k1=ief(i1)*ibi+ief(i3); k2=i2; k3=i4; k4=0; Idef_itype = 6
       elseif(ief(i2).gt.0.and.ief(i4).gt.0) then
        k1=ief(i2)*ibi+ief(i4); k2=i1; k3=i3; k4=0; Idef_itype = 7
       elseif(ief(i1).gt.0.and.ief(i4).gt.0) then
        k1=ief(i1)*ibi+ief(i4); k2=i3; k3=i2; k4=0; Idef_itype = 8
       elseif(ief(i2).gt.0.and.ief(i3).gt.0) then
        k1=ief(i3)*ibi+ief(i2); k2=i1; k3=i4; k4=0; Idef_itype = 9
       end if                                       

      elseif(io.gt.0.and.jo.eq.0) then

        if(ief(i1).gt.0) then
         k1=i2*ibi+i4; k2=i3; k3=ief(i1); k4=io; Idef_itype = 2
        elseif(ief(i2).gt.0) then
         k1=i1*ibi+i3; k2=i4; k3=ief(i2); k4=io; Idef_itype = 3
        elseif(ief(i3).gt.0) then
         k1=i2*ibi+i4; k2=i1; k3=ief(i3); k4=io; Idef_itype = 4
        elseif(ief(i4).gt.0) then
         k1=i1*ibi+i3; k2=i2; k3=ief(i4); k4=io; Idef_itype = 5
       else
         k1=i1*ibi+i3; k2=i2*ibi+i4; k3=io; k4=0; Idef_itype = 1
       end if

      elseif(io.gt.0.and.jo.gt.0) then

       k1=i1*ibi+i3; k2=i2*ibi+i4; k3=io; k4=jo; Idef_itype = 1

      end if

      if(Idef_itype.eq.0) Stop 'Idef_itype=0'

      End Function Idef_itype
