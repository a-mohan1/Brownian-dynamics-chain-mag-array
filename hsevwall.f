c     This subroutine implements hard sphere excluded volume interactions
c     between beads and channel walls

c     Feb 18 2007

c     Inputs
c     z = bead z-coords
c     nbeads = number of beads
c     a = bead radius
c     hchann = channel height (lower surface: z=-hchann/2, top surface: z=hchann/2) 

c     all distances in units of l (max spring length)

      subroutine hsevwall(z, nbeads, a, hchann)

      implicit none
      integer nbeads
      integer ibead
      real*8 z(nbeads), a, hchann

c     ibead = bead counter

c     iterate over beads

      do 100 ibead = 1, nbeads
c     check for bead-wall overlap
         if(z(ibead).lt.(-hchann/2.d0+a)) then   
            z(ibead) = -hchann/2.d0 + a
         elseif(z(ibead).gt.(hchann/2.d0-a)) then
            z(ibead) = hchann/2.d0 - a
         endif
 100  continue

      return
      end
