c     This subroutine iterates over beads, identifies obstacles
c     with which an overlap may occur for each bead, and
c     calls subroutine hsev1obst to implement hard sphere EV
c     interactions
c     for a SQUARE array of obstacles

c     Feb 19 07

c     x, y = array of bead xy coords
c     nbeads = number of beads
c     a = bead radius
c     Robs = obstacle radius
c     lspac = lattice spacing
c     xstart = x coord of plane at which obstacle array starts

      subroutine hsevsqarray(x, y, nbeads, a, Robs, lspac, xstart)

      implicit none
      integer nbeads
      real*8 x(nbeads), y(nbeads), a, Robs, lspac, xstart
      real*8 jbead, kbead, cx1, cx2, cy1, cy2
      integer ibead, j1, j2, k1, k2

c     ibead = counter over beads
c     jbead, kbead = +ve lattice row/ lattice col corresponding to bead pos (nonintegral)
c     j1, j2, k1, k2 = lattice rows/cols adjacent to bead 
c     cx1, cx2, cy1, cy2 = obstacle coords

c     iterate over beads
      do 500 ibead = 1, nbeads
         if(x(ibead).ge.xstart) then
           kbead = (x(ibead) - xstart)/lspac + 1.d0
           k1 = int(kbead)
           k2 = k1 + 1
         else
           k1 = 1
           k2 = 1
         endif
         cx1 = (k1 - 1.d0)*lspac + xstart
         cx2 = (k2 - 1.d0)*lspac + xstart
         if(y(ibead).ge.0.d0) then
            jbead = y(ibead)/lspac + 1.d0
         else 
            jbead = 1.d0 - y(ibead)/lspac
         endif
         j1 = int(jbead)
         j2 = j1 + 1
         cy1 = (j1 - 1.d0)*lspac 
         cy2 = (j2 - 1.d0)*lspac
         if(y(ibead).lt.0.d0) then
           cy1 = -1.d0 * cy1
           cy2 = -1.d0 * cy2
         endif
         call hsev1obst(x, y, cx1, cy1, a, Robs, nbeads, ibead)
         call hsev1obst(x, y, cx1, cy2, a, Robs, nbeads, ibead)
         call hsev1obst(x, y, cx2, cy1, a, Robs, nbeads, ibead)
         call hsev1obst(x, y, cx2, cy2, a, Robs, nbeads, ibead)
 500  continue    

      return
      end


