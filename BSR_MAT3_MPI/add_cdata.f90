!======================================================================
      Subroutine Add_cdata(ip,jp,C,j1,j2,j3,j4)
!======================================================================
!     add new data (C,j1,j2,j3,j4) to block (ip:jp) in module cmdata
!----------------------------------------------------------------------
      Use cmdata
      Implicit none
      Integer :: ip,jp, j1,j2,j3,j4, k,l,m, i,ii
      Real(8) :: C

! ... case when block is empty (jp < ip):

      if(jp.lt.ip) then
       cdata(ip)=C; k1(ip)=j1; k2(ip)=j2; k3(ip)=j3; k4(ip)=j4
       jp=ip; Return
      end if       

! ... search position (k) for new integral

      k=ip; l=jp;
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (j1.lt.K1(m)) then;       l = m - 1
      elseif(j1.gt.K1(m)) then;       k = m + 1
      else
       if    (j2.lt.K2(m)) then;      l = m - 1
       elseif(j2.gt.K2(m)) then;      k = m + 1
       else
        if    (j3.lt.K3(m)) then;     l = m - 1
        elseif(j3.gt.K3(m)) then;     k = m + 1
        else
         if    (j4.lt.K4(m)) then;    l = m - 1
         elseif(j4.gt.K4(m)) then;    k = m + 1
         else
          cdata(m)=cdata(m)+C; Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue 

! ... shift up the rest of data:

      Do i=jp,k,-1
       ii = i + 1
       cdata(ii)=cdata(i)
       K1(ii)=K1(i); K2(ii)=K2(i); K3(ii)=K3(i); K4(ii)=K4(i)
      End do

! ... add new integral:

      Cdata(k)=C; K1(k)=j1; K2(k)=j2; K3(k)=j3; K4(k)=j4; jp=jp+1

      End Subroutine Add_cdata
