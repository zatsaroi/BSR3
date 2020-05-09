!----------------------------------------------------------------------
!     target, rsol.nnn  --> rsol.all
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Character(80) :: AF_r = 'rsol.nnn',  AF_all = 'rsol.all' 
      Character(80) :: AF, AF1,AF2
      Integer, allocatable :: nchn(:), nstk(:)
      Real(8), allocatable :: v(:)

      Call Read_iarg('nlsp',nlsp)
      if(nlsp.le.0) Stop 'nlsp <= 0'

      iout = 2
      open(iout,file=AF_all,form='UNFORMATTED')

      Allocate(nchn(nlsp), nstk(nlsp))

      Do ilsp = 1,nlsp

       write(AF,'(a,i3.3)') 'rsol.',ilsp
       Call Check_file(AF)
       inp=1; open(inp,file=AF,form='UNFORMATTED')

       read(inp)  nhm,khm,kch,kcp,ns 
       read(inp)  (e,i=1,khm)

       Allocate(v(kch*ns))

       Do  j = 1,khm 
        read(inp) v
        write(iout) v
       End do 

       Deallocate(v)

       nchn(ilsp) = kch
       nstk(ilsp) = khm

      End do

      ipri = 3; open(ipri,file='rsol_collect.out')

      write(ipri,'(a)') 'nlsp (inast), number of partial waves:'
      write(ipri,'(i5)') nlsp

      write(ipri,'(a)') 'bspl_ndim:'
      write(ipri,'(i5)') ns

      nchmx = maxval(nchn)
      write(ipri,'(a)') 'nchmx, maximum number of channels:'
      write(ipri,'(i5)') nchmx

      nstmx = maxval(nstk)
      write(ipri,'(a)') 'nstmx, maximum number of RM solutions:'
      write(ipri,'(i5)') nstmx

      write(ipri,'(a)') 'nchn(1:nlsp) - number of channels for each partial wave:'
      write(ipri,'(10i10)') nchn

      write(ipri,'(a)') 'nstk(1:nlsp) - number of RM solutions for each partial wave:'
      write(ipri,'(10i10)') nstk


      End ! program

      