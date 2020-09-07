c     This subroutine counts collisions of the chain
c     with obstacles in a magnetic lattice, called at each time step
c
c     March 23 2007
c
c     x, y = arrays of bead xy coords
c     nbeads = number of beads
c     xmin, xmax = min and max chain x coords
c     ymin, ymax = min and max chain y coords
c     Nobs = number of obstacles
c     xobs, yobs = array of obstacle center xy coords
c     Nbinsx, Nbinsy = # bins along x and y coords
c     Nbins = total # bins
c     binsizex, binsizey = dimensions of bin along x, y 
c     binhead = array containing 1 obstacle of each bin
c     binlist = array containing remaining obstacles in bins
c     pointed to by binhead, separated by 0s
c     numcoll = # collisions in current traj, 
c     initialized to 0 for current traj in main program
c     maxcoll = max allowable #collisions per trajectory
c     (max size of arrays sobsx, sobsy)
c     sobsx, sobsy = saved array of obstacle xy coords with which
c     collisions have occurred in current traj.
c 
c     Periodic BCs are implemented by repeating the unit cell;
c     chain coords are unaffected
c     Origin is at channel lower left hand corner
c     A collision is defined to have occurred if a portion of the chain
c     is present in all 4 quadrants centered at the obst (subroutine collcheck)
c
c     A simultaneous collision with multiple obsts must not be counted mult times
c     To avoid this: check obstacles for hooked configs from min to max x
c     (dir of chain motion) in outer loop and for each fixed x, 
c     from min to max y in inner loop
c     and quit the subroutine after at most one collision per time step
c     is detected

      subroutine mcollcount(x, y, nbeads, xmin, xmax, ymin, ymax, 
     &                      Nobs, xobs, yobs, Nbinsx, Nbinsy, Nbins,
     &                      binsizex, binsizey, binhead, binlist, 
     &                      numcoll, maxcoll, sobsx, sobsy)

      implicit none
      integer nbeads, Nobs, Nbinsx, Nbinsy, Nbins
      integer binhead(Nbins), binlist(Nobs)
      integer numcoll, maxcoll
      integer mmin, mmax, nmin, nmax, j, k, jp, kp
      integer testbin, cobsind, icoll
      real*8 x(nbeads), y(nbeads), xmin, xmax, ymin, ymax
      real*8 xobs(Nobs), yobs(Nobs), binsizex, binsizey
      real*8 sobsx(maxcoll), sobsy(maxcoll)
      real*8 xcobs, ycobs
      logical ishooked, issaved
   
c     mmin, mmax, nmin, nmax = bin index row and col vals
c     of chain xy min and max coords (without periodic BCs)
c     j = counter from mmin to mmax (without periodic BCs)
c     k = counter from nmin to nmax (without periodic BCs)
c     jp, kp = j, k after implementing periodic BCs 
c     testbin = current bin index to be tested corresp to jp, kp
c     given by 1 + jp + kp*Nbinsx
c     cobsind = current obstacle index, 
c               ie, binhead(testbin), binlist(binhead(testbin)) etc.
c     icoll = counter over sobsx, sobsy elements
c     xcobs, ycobs = current obst coords corresp to cobsind
c     after translating the lattice to match j, k position
c     ishooked indicates whether a hooked config exists around xcobs, ycobs
c     (initialized in subroutine collcheck)
c     issaved indicates whether or not xcobs, ycobs is already saved

c     check if the chain max x coord has entered the array; if not, do nothing      
      if(xmax.lt.0.d0) then
        goto 2000
      endif

c     calculate mmin, mmax, nmin, nmax without periodic BCs
c     mmin must be at least 0 (since lattice starts at x=0)
c     if xmin is less than zero, must still check the obstacles 
c     along x = 0, i.e., m=0

      if(xmin.ge.0.d0) then
         mmin = int(xmin/binsizex)
      else
         mmin = 0
      endif
      mmax = int(xmax/binsizex)
      nmin = int(ymin/binsizey) 
      nmax = int(ymax/binsizey)

c     if ymax or ymin lie below 0, subtract 1 from nmin, nmax
c     (since int() rounds towards 0 and 
c      numbering of bin y-index k = -1, -2, .... below x-axis)

      if(ymax.lt.0.d0) then
        nmax = nmax - 1
        nmin = nmin - 1
      elseif(ymin.lt.0.d0) then
        nmin = nmin - 1
      endif

c     loop over all bins enclosed by mmin, mmax in x direction (j)
c     and nmin, nmax in y direction (k)
c     if j falls out of unit cell (ie, exceeds Nbinsx-1) then
c     implement periodic BCs to determine jp
c     j will be at least 0
c     if k falls out of unit cell (ie, exceeds Nbinsy-1 or falls below 0)
c     then implement periodic BCs to determine kp
c     determine bin index testbin corresp. to each jp, kp
c     check binhead(testbin), binlist(binhead(testbin)),.. etc. until 
c     all obstacles in bin have been checked for hooked configs
c     after first translating obstacle coords to correspond to j, k values
c     by periodically repeating the unit cell
c     if a hooked chain config is detected, check if the obst coords have already 
c     been saved in current traj; if not, save them and increment numcoll
c     repeat for all bins to be tested

      do 400 j = mmin, mmax
         jp = j
         if(j.gt.(Nbinsx-1)) then
           jp = mod(j, Nbinsx)
         endif
         do 200 k = nmin, nmax
            kp = k
            if(k.lt.0) then
              kp = k + (abs(int((k+1)/Nbinsy)) + 1) * Nbinsy
            elseif(k.gt.(Nbinsy-1)) then 
              kp = mod(k, Nbinsy)
            endif
            testbin = 1 + jp + kp*Nbinsx
            cobsind = binhead(testbin)
 100        continue
            if(cobsind.gt.0) then 
              xcobs = xobs(cobsind) + (j-jp)*binsizex
              ycobs = yobs(cobsind) + (k-kp)*binsizey
              call collcheck(x, y, nbeads, xcobs, ycobs, ishooked)
              if(ishooked) then
                issaved = .FALSE.
                icoll = numcoll
 110            continue
                if((icoll.ge.1).AND.(.NOT.issaved)) then
                   if((xcobs.eq.sobsx(icoll)).AND.
     &                (ycobs.eq.sobsy(icoll))) then
                     issaved =.TRUE.
                   endif
                   icoll = icoll - 1
                   goto 110
                endif
                if(.NOT.issaved) then 
                  numcoll = numcoll + 1
                  if(numcoll.gt.maxcoll) then
                    write(*, *) 'max number of collisions exceeded; bye'
                    stop
                  endif
                  sobsx(numcoll) = xcobs
                  sobsy(numcoll) = ycobs
                endif
c               hooked config. detected, quit for this time step
                goto 2000
              else
                cobsind = binlist(cobsind)
                goto 100
              endif
            endif
 200     continue
 400  continue   
  
 2000 continue
      return 
      end
