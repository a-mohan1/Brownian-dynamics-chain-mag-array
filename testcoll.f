c     this program tests collcount.f and collcheck.f
c     by comparing with xyzla3d1p1 movie (la3d1p1s50p.mpg)
c     (total dist covered = xstart + 50 lattice spacings)

      program testcoll

      implicit none
      integer nbeads, lattype, numcoll, maxcoll, ncols
      parameter(nbeads = 38)
      parameter(ncols = 101)
      parameter(maxcoll = 1000)
      parameter(lattype = 0)
      integer ibead, icol
      real*8 ntrap(ncols)
      real*8 x(nbeads), y(nbeads), z(nbeads)
      real*8 xmin, xmax, ymin, ymax, lspac, xstart
      real*8 obsx(maxcoll), obsy(maxcoll)
      parameter(lspac = 5.42d0)
      parameter(xstart = 3.d0)
      
      numcoll = 0

      do 5 icol = 1, ncols
         ntrap(icol) = 0.d0
 5    continue
      
      open(unit = 1, file = 'fmtxyzla3d1p1.dat', status = 'OLD')
  
 10   read(unit = 1, fmt = *, end = 200)    
      backspace(unit = 1)
      do 20 ibead = 1, nbeads
         read(1, *) x(ibead), y(ibead), z(ibead)
 20   continue

      xmin = x(1)
      xmax = x(1)
      ymin = y(1)
      ymax = y(1)    
      do 50 ibead = 1, nbeads
         if(xmin.gt.x(ibead)) then
            xmin = x(ibead)
         elseif(xmax.lt.x(ibead)) then
            xmax = x(ibead)
         endif
         if(ymin.gt.y(ibead)) then
            ymin = y(ibead)
         elseif(ymax.lt.y(ibead)) then
            ymax = y(ibead)
         endif
 50   continue     

      call collcount(x, y, nbeads, xmin, xmax, ymin, ymax, 
     &               lspac, lattype, xstart, numcoll, 
     &               obsx, obsy, maxcoll, ntrap, ncols)
  
      
      goto 10
 200  continue

      write(*, *) numcoll

      close(unit = 1)
      open(unit = 2, file = 'ntrapla3d1p1s50.dat', status = 'NEW')
      write(2, 201) (ntrap(icol), icol = 1, ncols)
      close(unit = 2, status = 'KEEP')

 201  format(1e12.6)
      stop 
      end
      
  
