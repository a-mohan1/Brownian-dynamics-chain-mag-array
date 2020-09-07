c     This subroutine calculates the x, y and z coords of beads
c     at time step tstep 
c     from the combined 3nbeads-by-one array r
c     r = x,y,z coordinates of beads in sequence
c     at time tstep

      subroutine r2xyz(r, nbeads, xt, yt, zt)

c     nbeads is a fixed parameter defined in the program
c     nbeads = number of beads

      implicit none
      integer nbeads, ibead
      real*8 xt(nbeads), yt(nbeads), zt(nbeads)
      real*8 r(3*nbeads)
      
      do 10 ibead = 1, nbeads
            xt(ibead) = r((ibead-1)*3+1)
            yt(ibead) = r((ibead-1)*3+2)
            zt(ibead) = r((ibead-1)*3+3) 
10    continue
        
      return
      end

