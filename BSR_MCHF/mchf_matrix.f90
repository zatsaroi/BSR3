!=====================================================================
      Subroutine mchf_matrix(i,hfm,rhs,hx) 
!=====================================================================
!     Set up the mchf-matrix hfm for orbital 'i' 
!     hx  is the matrix of extra terms for hf_nr;
!     rhs is the right-hand-side for orbitals that appear 
!     only once in the integral.
!----------------------------------------------------------------------
!     Integral          hfm               + hx                rhs
!                      
!     I(i,i)           I(., .)
!     I(i,j)                                            (1/2) I(.,.)j  
!                      
!     R(i,i;i,i)      2R(.,i;.,i)       4R(.,.;i,i)
!                      
!     R(i,i;i,j)       R(.,i;.,j)   (3/2)R(.,.;i,j)
!                 (1/2)R(.,.;i,j)                
!
!     R(i,i;a,b)       R(.,.;a,b)
!     R(i,a;i,b)       R(.,a;.,b)
!                      
!     R(i,a;b,c)                                       (1/2)R(.,a;.,c)b
!
!
!     all is devided on factor 2  ???
!----------------------------------------------------------------------
      Use bsr_mchf
      Use id4_data
       
      Implicit none
      Integer, intent(in)  :: i
      Real(8), intent(out) :: hfm(ns,ns),rhs(ns),hx(ns,ns)
      Real(8) :: d(ns,ks),dd(ns,ks), x(ns,ns),xx(ns,ns),  hlx(ns,ns), &
                 v(ns), CI(nbf), vv(ns), v1(ns), aa(ns,ks,nbf), iaa(nbf)

      Real(8) :: S,C,SS,t1,t2
      Integer :: int,io,j,ii,jj,k,i1,i2,i3,i4, m, met,ip
      Integer, external :: Ipointer
      Real(8), external :: WW, quadr, ZCB
  
      Call CPU_time(t1)

      hfm = 0.d0; rhs =0.d0; hx = 0.d0

! ... core-orbital case:

      if(i.le.ncore) then
       Call core_matrix (i,hfm,hx)
       Call CPU_time(t2); time_matrix = time_matrix + (t2-t1)
       Return
      end if

!----------------------------------------------------------------------
! ... one-electron HL integrals:

      Call hlm(lbs(i))
      Call Full_mat_sym(ns,ks,hl,hlx,'l') 

! ... interaction with core:

      if(ncore.gt.0) then

! ... direct interaction:

       Call MRK_cell(0)
       dd = 0.d0
       Do ip = 1,ncore
        C = -2.d0*(4*lbs(ip)+2)
        Call Density (ns,ks,d,p(1,ip),p(1,ip),'s')
        dd = dd + C*d
       End do

       Call Convol(ns,ks,d,dd,1,'s','s')
       Call UPDATE_HS(ns,ks,hlx,d,'s')          

! ... exchange interaction:

       Do k = 0, kmax
        Call MRK_cell(k)
        xx = 0.d0
        met = 0
        Do ip = 1,ncore
         C = (4*lbs(ip)+2) * ZCB(lbs(i),k,lbs(ip))
         if(C.eq.0.d0) Cycle
         Call Density (ns,ks,x,p(1,ip),p(1,ip),'x')
         xx = xx + C*x                    
         met = 1
        End do
        if(met.eq.0) Cycle
        Call Convol(ns,ks,x,xx,4,'s','s')
        Call UPDATE_HS(ns,ks,hlx,x,'x')          
       End do    

     end if  ! over ncore

!----------------------------------------------------------------------
! ... L-integrals:   to right side

      S = 0.d0; CI = 0.d0
      Do int = 1,Lint
       i1 = i1_Lint(int); i2 = i2_Lint(int)
       Do m = ip_Lint(int-1)+1,ip_Lint(int) 
        C = WW(ic_Lcoef(m),jc_Lcoef(m)) * L_coef(m)   
        if     (i1.eq.i.and.i2.eq.i) then;    S = S + C
        elseif (i1.eq.i.and.i2.ne.i) then;    CI(i2) = CI(i2) + C  
        elseif (i1.ne.i.and.i2.eq.i) then;    CI(i1) = CI(i1) + C 
        else;                                 Cycle
        end if
       End do
      End do

s = -2*s

      if(abs(S).lt.qsum_min) S = qsum_min  
      qsum(i) = S

      if(sum(abs(CI)).ne.0.d0) then             
       vv = 0.d0
       Do io = 1,nbf
        if(CI(io).eq.0.d0) Cycle
        vv = vv + p(:,io) * CI(io)/S/2.d0             
       End do
       rhs = rhs + matmul(hlx,vv)                      
      end if

!-----------------------------------------------------------------------------
! ... Collect all relevant Rk-integrals into module id4_data.
! ... Flags: -1 -> direct; -2 -> exchange; -3 -> hx;  >0 -> rhs.
! ... All angular coefficients are modified due to corresponding expansion 
! ... coefficients. We can do that for each k with id3_data ???

      Do k=0,kmax
       Call MRK_cell(k)

      nid = 0
      Do int = nk_int(k-1)+1,nk_int(k)
       i1=i1_int(int); i2=i2_int(int); i3=i3_int(int); i4=i4_int(int)
       Do m = ip_int(int-1)+1,ip_int(int) 
        C = WW(ic_coef(m),jc_coef(m))*Rk_coef(m)/S          
        if(abs(C).lt.eps_c) Cycle
        ii=0
        if(i.eq.i1) ii=ii+1000
        if(i.eq.i2) ii=ii+100
        if(i.eq.i3) ii=ii+10
        if(i.eq.i4) ii=ii+1
        if(ii.eq.0) Cycle

       Select case(ii)

        Case(1111);   Call Add_id4_data(k,-1,i1,i3,2*C)                      
                      Call Add_id4_data(k,-3,i1,i3,4*C)                     

        Case(1010);   Call Add_id4_data(k,-1,i2,i4,C)      
        Case(0101);   Call Add_id4_data(k,-1,i1,i3,C)

        Case(1001);   Call Add_id4_data(k,-2,i2,i3,C)      
        Case(1100);   Call Add_id4_data(k,-2,i3,i4,C)     
        Case(0110);   Call Add_id4_data(k,-2,i1,i4,C)
        Case(0011);   Call Add_id4_data(k,-2,i1,i2,C)

        Case(0111);   Call Add_id4_data(k,-1,i1,i3,C)
                      Call Add_id4_data(k,-2,i1,i3,C/2)
                      Call Add_id4_data(k,-3,i1,i3,3*C/2)
        Case(1011);   Call Add_id4_data(k,-1,i2,i4,C)
                      Call Add_id4_data(k,-2,i2,i4,C/2)
                      Call Add_id4_data(k,-3,i2,i4,3*C/2)
        Case(1101);   Call Add_id4_data(k,-1,i1,i3,C)
                      Call Add_id4_data(k,-2,i1,i3,C/2)
                      Call Add_id4_data(k,-3,i1,i3,3*C/2)
        Case(1110);   Call Add_id4_data(k,-1,i2,i4,C)
                      Call Add_id4_data(k,-2,i2,i4,C/2)
                      Call Add_id4_data(k,-3,i2,i4,3*C/2)

        Case(1000);   Call Add_id4_data(k,i3,i2,i4,C/2)         ! rhs          
        Case(0100);   Call Add_id4_data(k,i4,i1,i3,C/2)                
        Case(0010);   Call Add_id4_data(k,i1,i2,i4,C/2)                
        Case(0001);   Call Add_id4_data(k,i2,i1,i3,C/2)                

       End Select
       End do
      End do  ! over int

      if(nid.eq.0) Cycle

! ... add direct integrals

      dd = 0.d0; met=0
      Do j = 1,nid   
       if(id2(j).ne.-1) Cycle
       if(cid(j).eq.0.d0) Cycle 
       Call density(ns,ks,d,p(1,id3(j)),p(1,id4(j)),'s')

       dd = dd + cid(j)*d
       met = 1
      End do

      if(met.eq.1) then
       Call Convol(ns,ks,d,dd,1,'s','s')
       Call UPDATE_HS(ns,ks,hfm,d,'s')          
      end if

! ... add exchange integrals

      xx = 0.d0; met=0
      Do j = 1,nid   
       if(id2(j).ne.-2) Cycle
       if(cid(j).eq.0.d0) Cycle 
       Call density(ns,ks,x,p(1,id3(j)),p(1,id4(j)),'x')
       xx = xx + cid(j)*x                    

       met=1
      End do

      if(met.eq.1) then
       Call Convol(ns,ks,x,xx,4,'s','s')
       Call UPDATE_HS(ns,ks,hfm,x,'x')          
      end if

! ... add rhs 

      aa = 0.d0; iaa = 0
      Do j = 1,nid   
       if(id2(j).le. 0) Cycle
       Call density(ns,ks,d,p(1,id3(j)),p(1,id4(j)),'s')
       aa(:,:,id2(j)) = aa(:,:,id2(j)) + d*cid(j) 
       iaa(id2(j))=iaa(id2(j))+1
      End do

      Do j = 1,nbf 
       if(iaa(j).eq.0) Cycle
       d = aa(:,:,j)
       Call Convol(ns,ks,dd,d,1,'s','s')
       Call Full_mat_sym(ns,ks,dd,xx,'s')    ! ???
       v1 = MATMUL(xx,p(1:ns,j))
       rhs(1:ns) = rhs(1:ns) + v1(1:ns)
      End do

! ... add hx (for NR method):

      if(newton.ne.0) then
       xx = 0.d0; met=0
       Do j = 1,nid   
        if(id2(j).ne.-3) Cycle
        if(cid(j).eq.0.d0) Cycle
        Call density(ns,ks,x,p(1,id3(j)),p(1,id4(j)),'x')
        xx = xx + cid(j)*x
        met=1
       End do
       if(met.eq.1) then
        Call Convol(ns,ks,xx,x,4,'s','s')
        Call UPDATE_HS(ns,ks,hx,xx,'x')          
       end if
      end if

      End do ! over k

!----------------------------------------------------------------------
! ... convert rhs to hfm:

      srhs = SUM(abs(rhs))
      if(srhs.lt.1.d-8) then; rhs=0.d0; srhs=0.d0; end if

      if(irhs.gt.0.and.srhs.gt.0.d0) then
       v = p(:,i)
       SS = -SUM(rhs*v) 
       hfm = hfm + SS*BB
       vv = Matmul(BB,v)  
       Do ii = 1,ns; Do jj=1,ns 
        hfm(ii,jj) = hfm(ii,jj) + rhs(ii)*vv(jj) + rhs(jj)*vv(ii)
       End do; End do
       rhs = 0.d0
       srhs = 0.d0
      end if

! ... final matrix:

      hfm = -0.5d0*hlx + hfm ; rhs = -rhs                 !  -0.5  ???
      v = p(:,i)
      e(i)  = dot_product(p(:,i), matmul(hfm,p(:,i))-rhs)

      Call CPU_time(t2); time_matrix = time_matrix + (t2-t1)

      End Subroutine mchf_matrix



!=====================================================================
      Subroutine core_matrix(i,hfm,hx) 
!=====================================================================
!     Set up the mchf_matrix, hfm, for the core orbital i,
!     hx is required for the nr procedure.
!-----------------------------------------------------------------------
      Use bsr_mchf
      
      Implicit none
      Integer, intent(in) :: i
      Real(8), intent(out) :: hfm(ns,ns), hx(ns,ns)
      Real(8) :: d(ns,ks), dd(ns,ks), x(ns,ns), xx(ns,ns)
      Real(8) :: c
      Real(8), external :: a,b
      Integer :: j,k,met

! ... one-electron integral contribution

      Call hlm(lbs(i))
      Call Full_mat_sym(ns,ks,hl,hx,'l') 

      hfm = -0.5d0*hx;  hx = 0.d0

! ... two-electron contribution for different integrals

      Do k=0,kmax

!       if(meth.eq.'d') then
!         Call MRK_diff(k)
!       else
         Call MRK_cell(k)
!       end if

! ... Add the direct contribution   

       dd = 0.d0; met=0
       Do j = 1,nbf   
        c = a(i,j,k)
        if(i.ne.j) c = c *  qsum(j)
        if(i.eq.j) c = c *  (qsum(i)-1)
        if (c == 0.d0) cycle
        Call density(ns,ks,d,p(1,j),p(1,j),'s')

        dd = dd + c*d
        met = 1
       End do

       if(met.eq.1) then 
        Call Convol(ns,ks,d,dd,1,'s','s')
        Call UPDATE_HS(ns,ks,hfm,d,'s')          
       end if

! ... Add the exchange contribution 

       xx = 0.d0; met=0
       Do j = 1,nbf   
        c = b(i,j,k)
        if(i.ne.j) c = c *  qsum(j)
        if(i.eq.j) c = c *  (qsum(i)-1) 
        if (c == 0.d0) cycle
        Call density(ns,ks,x,p(1,j),p(1,j),'x')
        xx = xx + x*c

        met = 1
       End do

       if(met.eq.1) then 
        Call Convol(ns,ks,x,xx,4,'s','s')
        Call UPDATE_HS(ns,ks,hfm,x,'x')          
       end if

! ... add H_aa contributions
        
      if(newton.ne.0) then
       c = 2.d0 * a(i,i,k); if(c.eq.0.d0) Cycle
       c = c *  (qsum(i)-1)
       Call density(ns,ks,x,p(1,i),p(1,i),'x')
       xx = c*x
       Call Convol(ns,ks,x,xx,4,'s','s')
       Call UPDATE_HS(ns,ks,hx,x,'x')          
      end if

      End do ! over k

!      hfm = hfm / qsum(i); hx = hx / qsum(i)   

      End Subroutine core_matrix


!======================================================================
      Real(8) Function WW(ic,jc)
!======================================================================
! ... generalized value of WC(ic)*WC(jc) - product of expansion
! ... coeficients for configurations ic and jc
!======================================================================
      Use bsr_mchf
      
      Implicit none
      Integer, intent(in) :: ic,jc 
      Integer :: i,ip
      WW = 0.d0
      Do i=1,nlevels
       ip = ip_level(i)
       WW = WW + coefs(ip+ic)*coefs(ip+jc)*weight(i)  
      End do 
      if(ic.ne.jc) WW = WW + WW  

      End Function WW


