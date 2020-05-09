!======================================================================
!     bsr_dd - utility for calculation of dipole matrix elements
!              between R-matrx states saved in rsol.nnn files
!======================================================================
!
!     OUTPUT:   dd_nnn_mmm
!
!     SYSTEM CALLS:   MULT3,  BSR_DMAT3
!
!---------------------------------------------------------------------
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Integer :: nut=7; Character(20) :: AF_tar = 'target'
      Character(1)  :: blank = ' '
      Character(80) :: AS,BS,BI,BJ,AF

! ... read target and channels information:

      Call Check_file(AF_tar)
      Open(nut,file=AF_tar)
      Call R_target (nut)
      Call R_channels(nut)
      Close(nut)

      kpol = 1;  icount = 0

      klsp = 0; Call Read_iarg('klsp',klsp)
      if(klsp.le.0.or.klsp.gt.nlsp) klsp=nlsp

!----------------------------------------------------------------------

      Do i=1,klsp

       if(ipar(i).ne.1) Cycle

       write(BI,'(a,i3.3)') 'rsol.' ,i 
       if(Icheck_file(BI).eq.0) Cycle
       write(BI,'(a,i3.3)') 'cfg.' ,i 

      Do j=1,klsp

       write(BJ,'(a,i3.3)') 'rsol.' ,j 
       if(Icheck_file(BJ).eq.0) Cycle
       write(BJ,'(a,i3.3)') 'cfg.' ,j 

       if(ipar(i).eq.ipar(j)) Cycle

       if(ispar(i).ne.ispar(j)) Cycle

       if(ITRA(lpar(i),kpol,lpar(j)).eq.0) Cycle

       if(lpar(i)+lpar(j).lt.kpol) Cycle      

       write(AS,'(3(a,1x),a,2i3.3)') 'mult3',trim(BJ),trim(BI),'E1 mult_bnk_',i,j
       write(*,*) trim(AS)
       Call System(AS)

       write(AS,'(a,a,2i3.3,a)') 'cp ','mult_bnk_',i,j,' mult_bnk'
       Call System(AS)
   
       icount = icount + 1 

       write(AF,'(a,i3.3)') 'D',icount
       write(AS,'(a,a,a,a,a,a,a)')  &
         'bsr_dmat3 ',trim(BJ),' ',trim(BI),' b d  AF_dd=',trim(AF)
       write(*,*) trim(AS)

       Call System(AS)

       write(81,*)  icount, ispar(j),lpar(j),(1-ipar(j))/2, &
                            ispar(i),lpar(i),(1-ipar(i))/2

      End do; End do 

      open(82,file='D00',form='UNFORMATTED')

      write(82) icount
      rewind(81)
      Do i = 1,icount
       read(81,*) i0,i1,i2,i3,i4,i5,i6
       write(82)     i1,i2,i3,i4,i5,i6
      End do
      Close(82)
      Close(81,status='DELETE') 

      End ! program dd_bsr

