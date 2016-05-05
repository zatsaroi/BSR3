!======================================================================
      Subroutine Add_cdata(ip,jp,CL,CV,j1,j2,j3)
!======================================================================
!     add new coefficients with indexes (j1,j2,j3) to the list (ip:jp)
!----------------------------------------------------------------------
      Use cmdata
      Implicit none
      Integer:: ip,jp, j1,j2,j3, k,l,m, i,ii
      Real(8), intent(in) :: CL,CV

! ... first element in the list:

      if(jp.lt.ip) then
       cldata(ip) = CL; cvdata(ip) = CV;
       k1(ip)=j1; k2(ip)=j2; k3(ip)=j3; jp=ip; Return
      end if       

! ... search position for new integral

      k=ip; l=jp;
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if(j1.lt.K1(m)) then
       l = m - 1
      elseif(j1.gt.K1(m)) then
       k = m + 1
      else
       if(j2.lt.K2(m)) then
        l = m - 1
       elseif(j2.gt.K2(m)) then
        k = m + 1
       else
        if(j3.lt.K3(m)) then
         l = m - 1
        elseif(j3.gt.K3(m)) then
         k = m + 1
        else
         cldata(m) = cldata(m) + CL
         cvdata(m) = cvdata(m) + CV
         Return
        end if
       end if
      end if
      go to 1

    2 Continue 

      Do i=jp,k,-1
       ii = i + 1
       cldata(ii)=cldata(i); cvdata(ii)=cvdata(i)
       K1(ii)=K1(i); K2(ii)=K2(i); K3(ii)=K3(i)
      End do
      CLdata(k)=CL; CVdata(k)=CV;
      K1(k)=j1; K2(k)=j2; K3(k)=j3
      jp = jp + 1

      End Subroutine Add_cdata
