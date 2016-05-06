!======================================================================
      Subroutine sort_photo (nu,AF)
!======================================================================
!     Sorting results in photo.out  -->     photo.new
!     (from repeated results save only new !)
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)

      Integer, intent(in) :: nu
      Character(60), intent(in) :: AF

      Integer, allocatable :: ip(:)
      Real(8), allocatable :: E(:),EV(:),C(:),SL(:),SV(:),P(:),D(:)

! ... find number of energies
  
      rewind(nu)
      me = 0
    1 read(nu,AF,end=2,err=1) E1,EV1,C1,SL1,SV1,P1,D1 
      me=me+1
      go to 1
    2 write(*,*) ' me =',me
      Allocate(E(me),EV(me),C(me),SL(me),SV(me),P(me),D(me),IP(me))

      rewind(nu)
      ne = 0
   10 read(nu,AF,end=20,err=10) E1,EV1,C1,SL1,SV1,P1,D1 

      ie=0
      Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit;  End do
      if(ie.eq.0) then; ne=ne+1; ie=ne;  end if
      E(ie)=E1; EV(ie)=EV1; C(ie)=C1; SL(ie)=SL1;SV(ie)=SV1;
      P(ie)=P1; D(ie)=D1
 
      go to 10
   20 write(*,*) ' ne =',ne

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
       write(nu,AF) E(i),EV(i),C(i),SL(i),SV(i),P(i),D(i) 
      End do

      Deallocate(E,EV,C,SL,SV,P,D,IP)

      End Subroutine sort_photo
