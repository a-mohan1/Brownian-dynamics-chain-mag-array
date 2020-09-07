c     This subroutine assembles the full D tensor from the Djkmns
c     taking as inputs:
c     xc, yc, zc = bead coords at current time step istep
c     nbeads = number of beads (parameter)
c     istep = current time step
c     whicHI = which Green's tensor to use
c      0 = Free-draining
c      1 = Oseen-Burgers
c      2 = RPY
c      3 = Blake
c     a = bead radius
c     The result is stored in Dful
 
      subroutine Dfull(xc, yc, zc, nbeads, istep, whicHI, 
     &                 a, Dful)
      implicit none
      integer nbeads, istep, whicHI, m, n, j, k
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 a, Dful(3*nbeads, 3*nbeads)
      real*8  Djkmn
     
c     m = 'destination' bead index
c     n = 'source' bead index
c     j = direction of velocity at m
c     k = direction of force at n
c     Djkmn = one element of the diffusion tensor 
c             returned by function Djkmn

      do 40 m = 1, nbeads
       do 30 n = 1, nbeads
        do 20 j = 1, 3
         do 10 k= 1, 3
           Dful(3*(m-1)+j, 3*(n-1)+k) =
     &     Djkmn(xc, yc, zc, nbeads, istep, m, n, j, k, whicHI, a)  
 10      continue
 20     continue
 30    continue
 40   continue
 
      return
      end
