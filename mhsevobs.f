c     This subroutine iterates over beads,
c     identifies obstacles with which each bead might overlap,
c     and calls subroutine hsev1obst to implement HS EV with obstacles
c     in the given lattice of binned magnetic beads  
c
c     March 19 2007
c 
c     x, y = arrays of bead x, y coords
c     nbeads = number of beads
c     a = bead radius
c     binhead, binlist = array of obstacles in bins
c     Nbinsx = number of bins along x coord
c     Nbinsy = number of bins along y coord
c     Nbins = total # bins = size of array binhead
c     binsizex, binsizey = bin dimensions along x and y
c     xobs, yobs = arrays of obstacle center x, y coords
c     Nobs = number of obstacles
c     Robs = obstacle radius
c
c     Periodic BCs are implemented not by changing the bead coords  
c     but by identifying the periodic lattice beyond the channel
c     Note- origin is at channel lower left hand corner
c     int() rounds towards 0

      subroutine mhsevobs(x, y, nbeads, a, binhead, binlist, 
     &                    Nbinsx, Nbinsy, Nbins, binsizex, 
     &                    binsizey, xobs, yobs, Nobs, Robs)
      
      implicit none
      integer nbeads, Nobs, Nbinsx, Nbinsy, Nbins
      integer binhead(Nbins), binlist(Nobs)
      integer ibead, m, n, j, k, jpbc, kpbc, testbin, cobsind
      real*8 x(nbeads), y(nbeads), a, xobs(Nobs), yobs(Nobs), Robs
      real*8 binsizex, binsizey, xcobs, ycobs

c     ibead = bead counter
c     m, n = bin in which current bead would be placed without periodic BCs
c     j, k = counter over m-1, m, m+1 and n-1, n, n+1
c     jpbc, kpbc = j, k values after implementing periodic BCs
c     testbin = index of current bin to be tested
c     cobsind = obstacle indices in current bin being tested
c     xcobs, ycobs = current obstacle x, y coords for HSEV implementation
c
c     If current bead x coord < -Robs then the bead has not yet entered the array
c     and is not interacting with any obstacle, so do nothing
c     If current bead x coord exceeds -Robs but is below 0, the bead is not yet in
c     the channel but may interact with the obstacles on the y axis 
c     (m=0, j=0, jpbc=0)
c     If the bead has entered the array, then:
c     calculate m, which can exceed Nbinsx-1, and will be at least 0,
c     and n, which can exceed Nbinsy - 1 or be -ve
c     loop over the bin and all surrounding bins, i.e., 
c     j = {m-1, m, m+1} and k = {n-1, n, n+1} 
c     identify j and k with jpbc and kpbc within unit cell by implementing periodic
c     BCs st jpbc is in [0, Nbinsx-1] and kpbc in [0, Nbinsy-1]
c     Next: calculate bin index testbin and binhead obst index cobsind corresp
c     to jpbc, kpbc, and corresp. obstacle center coords, which may lie outside
c     the channel (depending on j, k)
c     test binhead(cobsind), binlist(binhead(cobsind)), ... for HS EV
c     until binlist element becomes 0, and then move to next j, k bin 
c     note: int rounds towards 0, so subtract 1 from n if y(ibead)<0
c     m is at least 0 so j is at least -1
c     but n may be a large -ve number in principle
c     so if k is -ve, simply adding Nbinsy will not give kpbc

      ibead = 1
 10   continue
      if(x(ibead).ge.0.d0) then
         m = int(x(ibead)/binsizex)
         n = int(y(ibead)/binsizey)
         if(y(ibead).lt.0.d0) then
           n = n - 1
         endif         
         do 400 j = m-1, m+1
            jpbc = j
            if(j.lt.0) then
               jpbc = j + Nbinsx
            elseif(j.gt.(Nbinsx-1)) then
               jpbc = mod(j, Nbinsx)
            endif
            do 200 k = n-1, n+1
               kpbc = k
               if(k.lt.0) then 
                 kpbc = k + (abs(int((k+1)/Nbinsy)) + 1) * Nbinsy
               elseif(k.gt.(Nbinsy-1)) then
                 kpbc = mod(k, Nbinsy)
               endif
               testbin = 1 + jpbc + kpbc*Nbinsx
               cobsind = binhead(testbin)
 100           continue
               if(cobsind.gt.0) then
                  xcobs = xobs(cobsind) + (j-jpbc)*binsizex
                  ycobs = yobs(cobsind) + (k-kpbc)*binsizey
                  call hsev1obst(x, y, xcobs, ycobs, a, Robs, nbeads,
     &                           ibead)
                  cobsind = binlist(cobsind)             
                  goto 100
               endif
 200        continue
 400     continue
      elseif(x(ibead).ge.(-Robs)) then
        j = 0
        n = int(y(ibead)/binsizey)
        if(y(ibead).lt.0.d0) then
          n = n - 1
        endif
        do 600 k = n-1, n+1
           kpbc = k
           if(k.lt.0) then 
             kpbc = k + (abs(int((k+1)/Nbinsy)) + 1) * Nbinsy
           elseif(k.gt.(Nbinsy-1)) then
             kpbc = mod(k, Nbinsy)
           endif
           testbin = 1 + j + kpbc*Nbinsx
           cobsind = binhead(testbin)
 500       continue
           if(cobsind.gt.0) then 
             xcobs = xobs(cobsind)
             ycobs = yobs(cobsind) + (k-kpbc)*binsizey
             call hsev1obst(x, y, xcobs, ycobs, a, Robs, nbeads,
     &                      ibead)
             cobsind = binlist(cobsind)
             goto 500
           endif
 600    continue
      endif
      ibead = ibead + 1
      if(ibead.le.nbeads) then
        goto 10
      endif

      return
      end
       
  
