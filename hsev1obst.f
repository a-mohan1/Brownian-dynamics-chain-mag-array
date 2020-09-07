c     This subroutine implements hard sphere excluded volume interactions
c     between the given bead and the given obstacle

c     subroutines hsevhexarray, hsevsqarray and mhsevobs iterate over beads 
c     and determine the xy coords of all obstacles with which an overlap 
c     may be possible

c     Feb 18 2007

c     Inputs
c     x, y = arrays of bead xy coords
c     nbeads = number of beads
c     ibead = current bead index
c     a = bead radius
c     cx, cy = obstacle axes xy coords
c     Robs = obstacle radius

c     all distances in units of l (max spring length)

      subroutine hsev1obst(x, y, cx, cy, a, Robs, nbeads, ibead)

      implicit none
      integer nbeads, ibead
      real*8 x(nbeads), y(nbeads), cx, cy
      real*8 a, Robs
      real*8 dist, ccdist, delx, dely

c     dist = distance by which bead must be moved in case of a bead-obstacle overlap
c     ccdist = shortest distance between bead center and obstacle axis
c     delx, dely = displacement to be made to bead center, if required

      ccdist = sqrt( (cx - x(ibead))**2 + (cy - y(ibead))**2 )
      dist = Robs + a - ccdist
      if(dist.gt.0.d0) then
         delx = dist * (cx - x(ibead))/ccdist
         dely = dist * (cy - y(ibead))/ccdist
         x(ibead) = x(ibead) - delx
         y(ibead) = y(ibead) - dely
      endif

      return
      end
