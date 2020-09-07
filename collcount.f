c     This subroutine checks if the current bead coords
c     indicate that a collision has taken place with an obstacle.
c     If so, the collision count is incremented and 
c     the obstacle coords are saved.
c     A collision is said to have occurred if the chain is wrapped 
c     around an obstacle, i.e., a portion of the chain is present in
c     all 4 quadrants centered at the obstacle center.
c     The obstacle coords from previous collisions are saved 
c     and checked to avoid counting a collision multiple times.
c     A collision with an obstacle can occur only once per trajectory.
c     Avoid counting a single collision with multiple obstacles
c     multiple times. 

c     Mar 8 07

c     x, y = bead x, y coords at current time step
c     nbeads = number of beads
c     xmin, xmax = min and max chain x coords
c     ymin, ymax = min and max chain y coords
c     lspac = lattice spacing
c     lattype = lattice type (hexagonal = 0, square = 1)
c     eqmdist = x coord of plane at which viewing area starts
c     numcoll = # collisions in current traj (size of obsx and obsy)
c     obsx, obsy = array of obstacle xy coords with which 
c            collisions have occurred in current traj
c     maxcoll = max size of arrays obsx, obsy
c     ncols = number of cols in viewing area
c     ntrap = array of # collisions with cols in viewing area

c     NOTE: int() rounds towards 0

      subroutine collcount(x, y, nbeads, xmin, xmax, ymin, ymax, 
     &                     lspac, lattype, eqmdist, numcoll, 
     &                     obsx, obsy, maxcoll, ntrap, ncols)

      implicit none
      integer nbeads, lattype, numcoll, maxcoll, ncols
      integer colind
      integer ibead, mulim, mllim, m, nulim, nllim, n, i
      real*8 ntrap(ncols)
      real*8 x(nbeads), y(nbeads), xmin, xmax, ymin, ymax
      real*8 obsx(maxcoll), obsy(maxcoll), xc, yc
      real*8 lspac, eqmdist, x0, vdist
      logical issaved, ishooked
  
c     x0 = starting x value for current obstacle row in viewing area
c     vdist = distance between lattice rows
c     xc, yc = current obstacle coords
c     ibead = bead counter
c     mulim, mllim = rows corresponding to min, max chain y coords
c     m = counter over obstacle rows
c     defined as yc = m*vdist
c     nulim, nllim = cols corresponding to min, max chain x coords
c     n = counter over obstacle cols in current row
c     defined as xc - x0 = n*lspac
c     colind = viewing area col index at which collision is detected
c     i = array counter
c     issaved = indicates whether current obstacle is already counted
c     ishooked = indicates whether a hooked config exists

      colind = -1
      
      if(lattype.eq.0) then
         vdist = lspac * sqrt(3.d0)/2.d0
      elseif(lattype.eq.1) then
         vdist = lspac
      else
         write(*, *) 'undefined lattice; bye'
         stop
      endif

c     calculate upper and lower lims for obstacle y coords
c     3 cases: the chain may be in the half plane y>0 or y<0 or 
c     one end may be in y>0 and one in y<0

c     calculate range of obstacle y coords
c     will use y = m*vdist
      if(ymax.le.0.d0) then
         mulim = int(ymax/vdist)
         mllim = int(ymin/vdist) - 1
      elseif(ymin.ge.0.d0) then
         mllim = int(ymin/vdist)
         mulim = int(ymax/vdist) + 1
      else
         mllim = int(ymin/vdist) - 1
         mulim = int(ymax/vdist) + 1
      endif

c     iterate over obstacle rows to be checked for collisions
      do 100 m = mllim, mulim
c        calculate y coord of current obstacle
         yc = m*vdist
c        calculate x coords of obstacles corresp to current y coord
                    
c        calculate starting x coord for current row
c        n=0 is starting obstacle in viewing area for all rows
c        note: m=0 for the central plane (previously, j=1)
         if( (mod(m, 2).ne.0).AND.(lattype.eq.0) ) then
           x0 = eqmdist + lspac/2.d0
         else
           x0 = eqmdist
         endif

c        calculate upper and lower obstacle x coord lims
         if(xmax.ge.x0) then
            nulim = int((xmax - x0)/lspac) + 1
         else
            nulim = 0
         endif
         if(xmin.ge.x0) then
            nllim = int((xmin - x0)/lspac)
         else
            nllim = 0
         endif

c        iterate over obstacle x coords
         do 80 n = nllim, nulim
            xc = n*lspac + x0

c           check if xc, yc are already saved 
            issaved = .FALSE.
            do 50 i = 1, numcoll
                  if( (xc.eq.obsx(i)).AND.(yc.eq.obsy(i)) ) then
                    issaved = .TRUE.
                  endif
 50         continue
        
c           if not, check if a hooked config has been formed
            if(.NOT.issaved) then
               call collcheck(x, y, nbeads, xc, yc, ishooked)
c              if so, increment # collisions and save obstacle coords
               if(ishooked) then
                 numcoll = numcoll + 1
                 if(numcoll.gt.maxcoll) then
                   write(*, *) 'max # collisions exceeded; bye'
                   stop
                 endif
                 obsx(numcoll) = xc
                 obsy(numcoll) = yc
                 if(colind.lt.0) then
                  if((mod(m, 2).eq.0).AND.(lattype.eq.0)) then
                    colind = 2*n + 1
                  elseif((mod(m, 2).ne.0).AND.(lattype.eq.0)) then
                    colind = 2*(n+1)
                  elseif(lattype.eq.1) then
                    colind = n+1
                  endif
                 endif
               endif
            endif
 80      continue
 100  continue

c     if a collision is detected in current step then increment collision count                 
c     if a collision involves multiple simulataneous obstacles only the lowest obstacle upstream
c     is counted

c     Note: does not work if a collision with 1 obst occured at one time step
c     and one chain end collided with another obst in next step
c     both collisions are counted

      if(colind.gt.0) then
         ntrap(colind) = ntrap(colind) + 1.d0
      endif

      return
      end



            
            
           

           
