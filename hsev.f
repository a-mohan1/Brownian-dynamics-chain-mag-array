c     This subroutine implements hard sphere excluded volume interactions
c     between beads and obstacles, and between beads and channel walls

c     useful if the xy coordinates of the axes of all obstacles are available
c     in arrays cx, cy, each of size nobst

c     Feb 16 2007

c     Inputs
c     x, y, z = bead coords
c     nbeads = number of beads
c     a = bead radius
c     cx, cy = obstacle axes xy coords
c     nobst = number of obstacles
c     Robs = obstacle radius
c     hchann = channel height (lower surface: z=-hchann/2, top surface: z=hchann/2) 

c     all distances in units of l (max spring length)

      subroutine hsev(x, y, z, nbeads, a, cx, cy, nobst, Robs, hchann)

      implicit none
      integer nbeads, nobst
      integer ibead, iobst
      real*8 x(nbeads), y(nbeads), z(nbeads), cx(nobst), cy(nobst)
      real*8 a, Robs, hchann
      real*8 dist, ccdist, delx, dely

c     ibead = bead counter
c     iobst = obstacle counter
c     dist = distance by which bead must be moved in case of a bead-obstacle overlap
c     ccdist = shortest distance between bead center and obstacle axis
c     delx, dely = displacement required to be made to bead center

c     iterate over beads

      do 100 ibead = 1, nbeads
c     check for bead-wall overlap
         if(z(ibead).lt.(-hchann/2.d0+a)) then   
            z(ibead) = -hchann/2.d0 + a
         elseif(z(ibead).gt.(hchann/2.d0-a)) then
            z(ibead) = hchann/2.d0 - a
         endif
c        iterate over obstacles and check for bead-obstacle overlaps
         do 50 iobst = 1, nobst
            ccdist = sqrt( (cx(iobst) - x(ibead))**2 +
     &              (cy(iobst) - y(ibead))**2 )
            dist = Robs + a - ccdist
            if(dist.gt.0.d0) then
              delx = dist * (cx(iobst) - x(ibead))/ccdist
              dely = dist * (cy(iobst) - y(ibead))/ccdist
              x(ibead) = x(ibead) - delx
              y(ibead) = y(ibead) - dely
            endif
 50      continue
 100  continue

      return
      end
