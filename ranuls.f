c     p 309 of Stochastic Processes in Polymeric Fluids by Ottinger
c     based on ran2 of Numerical Recipes in Fortran 
c     by Press, Teukolsky, Vetterling and Flannery, p 272 
c     combination of two generators using the multiplicative congruential algorithm
c       
      real*8 function ranuls()
c    
c     returns a random number uniformly distributed in [0, 1]
c
      parameter (in1=2147483563, ik1=40014, iq1=53668, ir1=12211,
     &           in2=2147483399, ik2=40692, iq2=52774, ir2=3791,
     &           ntab=32, an=1./in1, inm1=in1-1, ndiv=1+inm1/ntab)
      integer iv(ntab) 
      common /ranbls/ idum, idum2, iy, iv
c     
c     linear congruential generator 1
c
      k = idum/iq1
      idum = ik1*(idum-k*iq1) - k*ir1
      if (idum.lt.0) idum = idum+ in1
c
c     linear congruential generator 2
c
      k = idum2/iq2
      idum2 = ik2*(idum2-k*iq2) - k*ir2
      if (idum2.lt.0) idum2 = idum2 + in2
c
c     shuffling and subtracting
c
      j = 1 + iy/ndiv
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy.lt.1) iy = iy + inm1
c 
      ranuls = an*iy
      return
      end
