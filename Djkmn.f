c     This function returns an element Djkmn of
c     of the 3nbeads-by-3nbeads diffusion tensor
c     corresponding to the velocity induced in the j direction
c     at the location of bead m due to a stokeslet acting 
c     in the k direction
c     at the center of bead n

c     whicHI = which Green's tensor to use
c     0 = free-draining
c     1 = Oseen-Burgers
c     2 = RPY
c     3 = Blake

c     xc, yc, zc = bead positions at current time step istep
c     nbeads = number of beads (parameter)
c     m = 'destination' bead at which velocity is induced
c     n = 'source' bead approximated by a stokeslet acting at its center
c     j = direction of velocity
c     k = direction of force
c     istep = current time step

      real*8 function Djkmn(xc, yc, zc, nbeads, istep, m, n, j, k,   
     &                 whicHI, a)
      
      implicit none
      integer nbeads, istep, m, n, j, k, whicHI
      integer kdelta 
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 a, h, r(3), rmag, Rim(3), Rimmag, blake
      logical overlap      

      overlap = .FALSE.

c     if m = n, then Djkmn is given by Stokes' law  

      if (m.eq.n) then
       Djkmn = kdelta(j, k) 
c      return control    
       goto 200           
      endif

c     bead n to bead m relative vector 

      r(1) = xc(m) - xc(n)
      r(2) = yc(m) - yc(n)
      r(3) = zc(m) - zc(n)

c     center-to-center distance between beads m and n

      rmag = sqrt(r(1)**2 + r(2)**2 + r(3)**2)

c     check if the beads overlap

      if (rmag.lt.(2.d0*a)) then
          overlap = .TRUE.
c          write(*, *) 'overlapping beads detected: ', m, n, istep
          if (whicHI.ne.0.AND.whicHI.ne.2) then
             write(*, *) 'EV has not been implemented correctly'
             write(*, *) 'bye'
             stop
          endif
      endif

c     image of bead n to bead m relative vector
c     image of bead n is the reflection of bead n wrt the xy plane

      Rim(1) = xc(m) - xc(n)
      Rim(2) = yc(m) - yc(n)
      Rim(3) = zc(m) + zc(n)

c     distance between bead m and image of bead n
      Rimmag = sqrt(Rim(1)**2 + Rim(2)**2 + Rim(3)**2)

c     z-displacement of bead n stokeslet from the wall 
      h = zc(n) 

c     if m is not equal to n, calculate the specified tensor element
      if (whicHI.eq.0) then 
            Djkmn = 0
      elseif (whicHI.eq.1) then
            Djkmn = (3.d0*a)/(4.d0*rmag)*
     &      (kdelta(j, k) + r(j)*r(k)/rmag**2)
      elseif (whicHI.eq.2) then
            if(overlap) then
               Djkmn = (1-9.d0*rmag/(32.d0*a))*kdelta(j, k)
     &         + 3.d0/(32.d0*a*rmag)*r(j)*r(k)
            else
               Djkmn = 3.d0*a/(4.d0*rmag) *
     &         ((1+2.d0*a*a/(3.d0*rmag*rmag))*kdelta(j, k)
     &          + (1-2.d0*a*a/(rmag*rmag))*r(j)*r(k)/(rmag*rmag) )
            endif
      elseif (whicHI.eq.3) then
            Djkmn = blake(j, k, r, rmag, Rim, Rimmag, h, a)
      else
            write(*, *) 'unknown HI'
            write(*, *) 'bye'
            stop
      endif                

 

 200  continue

      return
      end


      

      
      
