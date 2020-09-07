c     This subroutine takes arrays for the x, y and z coords of beads 
c     at the current time step
c     and combines them into a 3nbeads-by-1 array called r
c     r = x,y,z coordinates of beads, one bead after the other
c     Arguments:
c     xc, yc, zc = bead coords 
c     nbeads = number of beads (parameter)

      subroutine xyz2r(xc, yc, zc, nbeads, r)

      implicit none
      integer nbeads, ibead
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 r(3*nbeads) 

      do 10 ibead = 1, nbeads
            r((ibead-1)*3+1) = xc(ibead)
            r((ibead-1)*3+2) = yc(ibead)
            r((ibead-1)*3+3) = zc(ibead)
10    continue
       
      return
      end

