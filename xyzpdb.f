c     This program converts xyz*.dat files into appropriate .pdb format
c     Max. number of beads = 101

      program xyzpdb

      implicit none
      integer nbeads, ibead
      double precision x(101), y(101), z(101)

 1    format('ATOM   ',i4,'  N   ALA     1     ',f7.2,' ',f7.2,' ',
     +     f7.2,'  1.00  7.00      MAIN')
 2    format('PSF')
 3    format(1x,i5,' !NATOM')
 4    format(2x,'  ',i4,' MAIN 1    ALA  N    N   0   45           0')
 5    format(2x,i4,' !NTITLE')
 6    format('END')
 7    format('ATOM  ',i5,'  C   ALA     1     ',f7.2,' ',f7.2,' ',
     +     f7.2,'  1.00  7.00      MAIN')
 8    format(1x,'  ',i5,' MAIN 1    ALA  C    C   0   45           0')
 999  format(2x,'  ',i4,' MAIN 1    ALA  O    O   0   45           0')
 998  format('ATOM   ',i4,'  O   ALA     1     ',f7.2,' ',f7.2,' ',
     +     f7.2,'  1.00  7.00      MAIN')

      open(unit=1,file='xyz61p10fdwt.dat',status='old')
      open(unit=2,file='movie61p10fdwt.psf',status='new')
      open(unit=3,file='movie61p10fdwt.pdb',status='new')

      read(1,*) nbeads
      backspace(unit=1)
      if(nbeads.gt.101)then
         write(*, *) 'Number of beads exceeds maximum limit (101)'
         write(*, *) 'bye'
         stop
      endif

********************************************************
*     PSF file
      write (2,2)
      write (2,*) " "
      write (2,5) 1
      write (2,*) 'REMARKS ==========='
      write (2,3) nbeads
      do 15 ibead = 1, nbeads
         if (ibead.le.nbeads) then
              write (2,4) ibead
         else
              write (2,8) ibead
         endif
 15   continue
********************************************************

c     read current frame from xyz file
 20   read(unit = 1, fmt = *, end = 200)
      do 25 ibead = 1, nbeads
            read(1, *) x(ibead), y(ibead), z(ibead)
c           write xyz coords to .pdb file
            write(3, 1) ibead, x(ibead), y(ibead), z(ibead) 
 25   continue
c     write 'END' to signify the end of the frame
c     and continue to next frame
      write(3,6)
      goto 20
 200  continue

      close(unit = 1)
      close(unit = 2)
      close(unit = 3)

      stop      
      end

