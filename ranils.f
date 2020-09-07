c     p 309 of Stochastic Processes in Polymeric Fluids by Ottinger
c     based on ran2 of Numerical Recipes in Fortran 
c     by Press, Teukolsky, Vetterling and Flannery, p 272 
c     combination of two generators using the multiplicative congruential algorithm
c
      subroutine ranils(iseed)
c     ranils initializes function ranuls
c     choice of iseed 0 <= iseed <= 2E9, iseed is an integer
c      
      parameter (in=2147483563, ik=40014, iq=53668, ir=12211, ntab=32)
      integer iv(ntab)
c
      common /ranbls/ idum, idum2, iy, iv
c
c     initial seeds for two random number generators
c 
      idum = iseed + 123456789
      idum2 = idum
c
c     load the shuffle table after 8 warm-ups
c      
      do 10 j = ntab+8, 1, -1
         k = idum/iq
         idum = ik*(idum-k*iq) - k*ir
         if (idum.lt.0) idum = idum + in
         if (j.le.ntab) iv(j) = idum
10    continue 
c      
      iy = iv(1)
c 
      return
      end
