c     This function returns a uniform random deviate in [0, 1]
      real*8 function rannum(whicRNG, ranflag, jseed)
c      
c     whicRNG = which random number generator to use
c     0 = Ottinger's ran2 (ranils to initialize, ranuls 
c     to generate a random number)
c     1 = Knuth's ran3
c     ranflag = 0 if this is the first call to RNG
c     ranflag = 1 otherwise
c    
      implicit none
      integer whicRNG, ranflag, jseed
      integer idum, idumss
      real*8 ranuls
      real ran3 
c     
c     if whicRNG = 0, jseed must be an integer between 0 and 2E9
c     if whicRNG = 1, jseed must be a negative integer
c     for first call to ran3, idum is jseed, a negative integer
c     for subsequent calls, idum is idumss, a positive integer
c     the numerical value of idumss is never used by ran3
  
      parameter (idumss = 120)

c     check to see if this is the first call to RNG
c     if so, check whether to use ranuls or ran3
c     if ranuls, initialize and generate
c     if ran3, set idum to a negative integer (jseed)
c     and then call ran3
c     finally, increment ranflag 
c     and return control
c 
      if(ranflag.eq.0) then
        if(whicRNG.eq.0) then
           call ranils(jseed)
           rannum = ranuls()
        else
           idum = jseed
           rannum = dble(ran3(idum))
        endif
        ranflag = 1
c       emergency exit, return control to calling routine        
        return
      endif
c
c     not the first call to RNG
c       
      if (whicRNG.eq.0) then
c     already initialized above 
         rannum = ranuls()
      else
c        MUST NOT REINITIALIZE
c        must now set idum to a positive integer (idumss)
         idum = idumss
         rannum = dble(ran3(idum))
      endif  

      return 
      end
