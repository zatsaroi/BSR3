!=====================================================================
      Subroutine hf_matrix(i,hfm) 
!=====================================================================
!     Set up the hf_matrix "hfm" for orbital i
!     Call(s):  a,b - angular coeficient routines
!-----------------------------------------------------------------------
      Use bsr_hf 
      Use hf_orbitals
      
      Implicit none
      Integer, intent(in) :: i
      Real(8), intent(out) :: hfm(ns,ns)
      Real(8) :: d(ns,ks),dd(ns,ks),x(ns,ns),xx(ns,ns)
      Real(8) :: c,t1,t2
      Real(8), external :: a,b,   BVMV, rk, SUM_AmB
      Integer :: j,k,met

      Call CPU_time(t1)

! ... one-electron integral contribution

      Call hlm(lbs(i))
      Call Full_mat_sym(ns,ks,hl,x,'l') 
      hfm = -qsum(i)*x/2

! ... two-electron contribution for direct and exchage integrals
! ... this contribution is devided on 4 parts: pppp, pqpq, qpqp, qqqq
! ... and added separately

      Do k=0,kmax
                                    
      Call mrk_cell(k)

! ... direct contribution    

       dd = 0.d0; met=0
       Do j = 1,nbf   
        c = a(i,j,k); if(c.eq.0.d0) cycle; if(i.eq.j) c=c+c      
        Call density(ns,ks,d,p(1,j),p(1,j),'s')
        dd = dd + c*d; met=met+1
       End do
       if(met.gt.0) Call Convol(ns,ks,d,dd,1,'s','s')
       if(met.gt.0) Call Update_hs(ns,ks,hfm,d,'s')          

! ... exchange contribution 

       xx = 0.d0; met=0
       Do j = 1,nbf   
        c = b(i,j,k); if (c.eq.0.d0) cycle
        Call density(ns,ks,x,p(1,j),p(1,j),'x')
        xx = xx + x*c; met=met+1
       End do
       if(met.gt.0) Call Convol(ns,ks,x,xx,4,'s','s')
       if(met.gt.0) Call Update_hs(ns,ks,hfm,x,'x')          

      End do ! over k

      hfm = hfm / qsum(i)

      Call CPU_time(t2)
      time_hf_matrix = time_hf_matrix + t2-t1

      End Subroutine hf_matrix


