!======================================================================
      Subroutine Add_mvc
!======================================================================
! ... add mass-velocity corrections:
!----------------------------------------------------------------------
      Use bsr_hd
      Use blacs
      Use channel,      only: iptar, lch
      Use target,       only: Etarg
      Use spline_param, only: ns,ks

      Implicit none
      Real(8) :: S, vc(ns,ks)
      Real(8), external :: BVMV
      Integer :: i,j, i1,i2, it,l, idim, ich, ii

      Call CPU_time(t0)

      Do ich = 1,kch

       if(io_processor) then
        it=iptar(ich);  l=lch(ich);  Call mvcv(l,vc)
        i1=ipsol(ich-1)+1; i2=ipsol(ich); ii = ipsol(ich-1) 
        add = zero
        Do i=i1,i2; if(i-i1.gt.jmvc) Cycle; ! if(bval(i).gt.Etarg(it)) Cycle 
        Do j=i1,i; if(j-i1.gt.jmvc) Cycle;  ! if(bval(j).gt.Etarg(it)) Cycle 
          S =  - 0.5*BVMV(ns,ks,vc,'s',bb(:,i),bb(:,j))
          add(i-ii,j-ii) = S
          add(j-ii,i-ii) = S
        End do; End do
       end if

       idim=ipsol(ich)-ipsol(ich-1); i1=ipsol(ich-1)+1

       Call pdgeadd ('notrans', idim, idim,     &
                      one,  add,  1,  1, descadd, &
                      one,    a, i1, i1, desca  )
      End do

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'Add_mvc:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'Add_mvc:,', (t1-t0)/60, ' min.'
      end if

      End Subroutine Add_mvc
