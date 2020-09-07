c     This program tests subroutine hsevhexarray

c     Feb 25 07

c     x, y = array of bead xy coords
c     nbeads = number of beads
c     a = bead radius
c     Robs = obstacle radius
c     lspac = lattice spacing
c     xstart = x coord of plane at which obstacle array starts

      program testhexarray

      implicit none
      integer nbeads  
      
      parameter (nbeads = 1)
      
      real*8 x(nbeads), y(nbeads), a, Robs, lspac, xstart
      real*8 PI, theta, jbead, cx, cy, x0
      integer ibead, jvals(2), k, m

      parameter (PI = 3.14159)
      parameter (theta = PI/3.d0)
      
      x(1) = 9.0
      y(1) = -1
      lspac = 1.2
      xstart = 4.0;
      a = 0.2;
      Robs = 1.0;

c     ibead = counter over beads
c     jbead = +ve lattice row corresponding to bead pos (nonintegral)
c     jvals = lattice rows adjacent to bead 
c     k = int position/lattice spacing
c     m = counter over jvals
c     cx, cy = current obstacle coords
c     x0 = starting x coord for current row
c     theta = lattice angle (subtended at hexagon center by a side)

c     iterate over beads
      do 500 ibead = 1, nbeads
         if(y(ibead).ge.0.d0) then
            jbead = y(ibead)/lspac/sin(theta) + 1.d0
         else
            jbead = 1.d0 - y(ibead)/lspac/sin(theta) 
         endif

         jvals(1) = int(jbead)
         jvals(2) = int(jbead) + 1

         do 400 m = 1, 2
            cy = (jvals(m) - 1.d0)*lspac*sin(theta)
            if(y(ibead).lt.0.d0) then
              cy = -1.d0*cy
            endif
            if(mod(jvals(m), 2).ne.0) then
              x0 = xstart 
            else 
              x0 = xstart + lspac*cos(theta)
            endif
            if(x(ibead).le.(x0+a)) then
              cx = x0;
              write(*, *) cx, cy
            else
              k = int((x(ibead)-x0)/lspac)
              cx = x0 + lspac*k
              write(*, *) cx, cy
              cx = cx + lspac
              write(*, *) cx, cy
            endif
 400     continue 
 500  continue 

      stop
      end


