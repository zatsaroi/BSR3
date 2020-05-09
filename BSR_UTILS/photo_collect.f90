!======================================================================
!     UTILITY   P H O T O _ C O L L E C T
!=====================================================================
 
      IMPLICIT REAL(8) (A-H,O-Z)

      Integer, parameter :: ma = 80
      Integer :: ics = 16; Character(ma) :: AF_out  = 'bsr_phot.nnn'
      Integer :: iph = 17; Character(ma) :: AF_ph   = 'photo.nnn'
      Character(ma) :: AF,BF,AS   

      Character(60) :: AFORM

      klsp = 1
      Call Read_iarg('klsp',klsp)
      nprocs = 1024
      Call Read_iarg('np',nprocs)

      ii=INDEX(AF_out,'.',BACK=.TRUE.)
      write(AF,'(a,i3.3)') AF_out(1:ii),klsp
      Do i=0,nprocs-1
        write(BF,'(a,i4.4)') AF_out(1:ii),i
        write(AS,'(a,a,a,a)') 'cat ',trim(BF),' >> ',trim(AF)
        Call System(AS)
      End do
     
       ii=INDEX(AF_ph,'.',BACK=.TRUE.)
       write(AF,'(a,i3.3)') AF_ph(1:ii),klsp
       Do i=0,nprocs-1
         write(BF,'(a,i4.4)') AF_ph(1:ii),i
         write(AS,'(a,a,a,a)') 'cat ',trim(BF),' >> ',trim(AF)
         Call System(AS)
       End do

       open(iph,file=AF)
       Call Sort_photo(iph,AFORM)
 
      END  ! program photo_collect 
 

!======================================================================
      Subroutine sort_photo (nu,AF)
!======================================================================
!     Sorting results in  photo.out
!     (from repeated results save only new ones ! )
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Integer, Intent(in) :: nu
      Character(60), Intent(in) :: AF
      Character(200) :: AS

      Integer,Allocatable, Dimension(:) :: ip
      Real(8),Allocatable, Dimension(:) :: E,EV,C,SL,SV,P,D

! ... find number of energies

      rewind(nu)
      me = 0
      Do 
       read(nu,'(a)',end=2) AS
       read(AS,*,err=1) E1,EV1,C1,SL1,SV1,P1,D1 
       me=me+1
    1  Continue
      End do
    2 if(me.eq.0) Return 
      Allocate(E(me),EV(me),C(me),SL(me),SV(me),P(me),D(me),IP(me))

      rewind(nu)
      ne = 0
      Do 
       read(nu,'(a)',end=20) AS
       read(AS,*,err=10) E1,EV1,C1,SL1,SV1,P1,D1 
       ie=0
       Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit;  End do
       if(ie.eq.0) then; ne=ne+1; ie=ne;  end if
       E(ie)=E1; EV(ie)=EV1; C(ie)=C1; SL(ie)=SL1;SV(ie)=SV1;
       P(ie)=P1; D(ie)=D1
   10 Continue
      End do
   20 if(ne.eq.0) Return

      Call SORTR(ne,E,IP)

! ... adjust the phase:

      Do j=2,ne; i=IP(j); i1=IP(j-1)

    3 PP = P(i) + 1.d0
      a1=abs(P(i)-P(i1))
      a2=abs(PP  -P(i1)) 
      if(a2.lt.a1) then
       P(i) = PP
       go to 3
      end if

    4 PP = P(i) - 1.d0
      a1=abs(P(i)-P(i1))
      a2=abs(PP  -P(i1))
      if(a2.lt.a1) then
       P(i) = PP
       go to 4
      end if

      End do

! ... output results:

      rewind(nu)
      Do j=1,ne; i=ip(j)
       write(nu,'(f14.8,f14.7,f12.2,2f15.5,f12.6,f8.3)') & 
         E(i),EV(i),C(i),SL(i),SV(i),P(i),D(i) 
      End do

      Deallocate(E,EV,C,SL,SV,P,D,IP)

      End Subroutine sort_photo
