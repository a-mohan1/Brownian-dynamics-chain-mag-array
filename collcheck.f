c     This subroutine checks whether the chain has
c     formed a hooked config around the given obstacle
c
c     Feb 28 07
c
c     x, y = current bead xy coords
c     nbeads = number of beads
c     xc, yc = current obstacle xy coords
c     ishooked indicates whether the chain is hooked

      subroutine collcheck(x, y, nbeads, xc, yc, ishooked)

      implicit none
      integer nbeads
      integer flaglr, flagur, flagll, flagul, ibead
      real*8 x(nbeads), y(nbeads), xc, yc
      logical ishooked
      
c     ibead = bead counter
c     flaglr indicates whether a portion of the chain is in the 
c     lower right hand corner of the obstacle
c     flagur indicates whether a portion of the chain is in the 
c     upper right hand corner of the obstacle
c     flagll indicates whether a portion of the chain is in the
c     lower left hand corner of the obstacle
c     flagul indicates whether a portion of the chain is in the 
c     upper left hand corner of the obstacle
c     All 4 conditions hold simultaneously in a hooked config

c     initialization
      ishooked = .FALSE.
      flaglr = 0
      flagur = 0
      flagll = 0
      flagul = 0

c     iterate over beads until either a hooked config is detected
c     or all beads have been checked

      ibead = 1
 10   continue
      if(x(ibead).gt.xc) then
        if(y(ibead).gt.yc) then
            flagur = 1
        elseif(y(ibead).lt.yc) then
            flaglr = 1
        endif
      else
        if(y(ibead).gt.yc) then
            flagul = 1
        elseif(y(ibead).lt.yc) then
            flagll = 1
        endif 
      endif
        
      if(((flagur*flaglr*flagul*flagll).eq.0).AND.(ibead.lt.nbeads)) 
     &  then
        ibead = ibead + 1
        goto 10
      endif

      if((flaglr*flagur*flagll*flagul).eq.1) then
        ishooked = .TRUE.
      endif

      return
      end
