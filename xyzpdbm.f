c     This program generates .pdb and .psf files for VMD
c     Used for the simulation of chain motion in a magnetic obstacle array
c     (lens is the length scale in microns)
c     obstacle coords from movie.xyz are in units of obstacle dia (1 micron)
c     Max. number of beads = 1001
c     Max. number of obstacles = 10,000
c     beads = N atoms
c     obstacles = C atoms

      program xyzpdbm

      implicit none
      integer maxnbeads, maxnobs
      parameter (maxnbeads = 1001)
      parameter (maxnobs = 10000)
      integer nbeads, Nobs, Npart, i
      character*32 header
      real*8 lens 
      parameter (lens = 0.554d0)
      double precision xb(maxnbeads), yb(maxnbeads), zb(maxnbeads)
      double precision xobs(maxnobs), yobs(maxnobs), zobs(maxnobs)
      double precision xcm, ycm, zcm, xtemp, ytemp

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

      open(unit = 1, file = 'rla3d1p5_1.dat', status = 'OLD')
      open(unit = 2, file = 'movie.xyz', status = 'OLD')
      open(unit = 3, file = 'mla3d1p5m.psf', status = 'NEW')
      open(unit = 4, file = 'mla3d1p5m.pdb', status = 'NEW')

      read(1, *) nbeads
      backspace(unit = 1)
      if(nbeads.gt.maxnbeads) then
         write(*, *) 'Number of beads exceeds maximum limit'
         write(*, *) 'bye'
         stop
      endif

 10   read(2, *) header
      if(header.eq.'#') goto 10
      backspace(unit = 2)
      read(2, *) Nobs
      if(Nobs.gt.maxnobs) then
         write(*, *) 'Number of obstacles exceeds max limit'
         write(*, *) 'bye'
         stop
      endif
      read(2, *) header
      do 12 i = 1, Nobs
         read(2, *) header, xtemp, ytemp, zobs(i)
         xobs(i) = xtemp/lens
         yobs(i) = ytemp/lens
 12   continue
      Npart = nbeads + Nobs

********************************************************
*     PSF file
      write (3, 2)
      write (3, *) " "
      write (3, 5) 1
      write (3, *) 'REMARKS ==========='
      write (3, 3) Npart
      do 15 i = 1, Npart
         if (i.le.nbeads) then
              write (3, 4) i
         else
              write (3, 8) i
         endif
 15   continue
********************************************************
*     PDB file

c     read current frame from file containing bead coords
 20   read(unit = 1, fmt = *, end = 200)
c     compute chain center of mass coords
      xcm = 0.d0
      ycm = 0.d0
      zcm = 0.d0
      do 25 i = 1, nbeads
            read(1, *) xb(i), yb(i), zb(i)
            xcm = xcm + xb(i)
            ycm = ycm + yb(i)
            zcm = zcm + zb(i)
 25   continue      
      xcm = xcm/nbeads
      ycm = ycm/nbeads
      zcm = zcm/nbeads
c     write bead coords wrt chain cm to pdb file
      do 30 i = 1, nbeads
            write(4, 1) i, xb(i)-xcm, yb(i)-ycm, zb(i)-zcm 
 30   continue
c     write obstacle coords wrt chain cm to pdb file
      do 40 i = 1, Nobs
            write(4, 7) nbeads+i, xobs(i)-xcm, yobs(i)-ycm,
     +                              zobs(i)-zcm
 40   continue      
c     write 'END' to signify the end of the frame
c     and continue to next frame
      write(4, 6)
      goto 20
 200  continue

      close(unit = 1)
      close(unit = 2)
      close(unit = 3)
      close(unit = 4)

      stop      
      end

