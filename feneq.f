c     This function calculates the root of the equation
c     Q**3 - R*Q**2 - (1+2*delt*Hspr)*Q + R == 0
c     lying in (0, 1)
c     Inputs:
c     R = magnitude of right hand side vector
c     Hspr = spring constant
c     delt = time step
      
      real*8 function feneq(R, Hspr, delt)
      implicit none
      integer maxit, iter
      real*8 R, Hspr, delt
      real*8 relerr, Qi, Qf, delQ, c3, c2, c1, c0
      parameter (maxit = 20, relerr = 1.d-10)

c     maxit = maximum number of Newton's iterations allowed
c     iter = Newton's iterations counter
c     relerr = maximum acceptable relative error in root
c     Qi = initial guess for root
c     Qf = root obtained after one Newton's iteration
c     delQ from Newton's method leads to Qf = Qi - delQ
c     c3, c2, c1, c0 = coeffs in cubic eq
      c3 = 1.0 
      c2 = -R
      c1 = -(1.0 + 2.0*delt*Hspr)
      c0 = R

c     Initialize starting guess: 
c     For R<1, use asymptotic soln for R<<1
c     For R>1, use asymptotic soln for R>>1
c     (if the latter becomes negative, use the former)
      
      if(R.lt.1.0) then
       Qi = R/(1.0 + 2.0*delt*Hspr)
      else
       Qi = 1.0 - delt*Hspr/R
       if(Qi.le.0.0) then
        Qi = R/(1.0 + 2.0*delt*Hspr)
       endif
      endif

c     Perform Newton's iterations
c     Qf = Qi - f(Qi)/fprime(Qi)
c     where f(Q) == 0 is the cubic eq and fprime is df/dQ
c     delQ = f(Q)/fprime(Q)

      do 10 iter = 1, maxit
         delQ = (c3*Qi**3 + c2*Qi**2 + c1*Qi + c0)/
     &           (3.0*c3*Qi**2 + 2.0*c2*Qi + c1)
         Qf = Qi - delQ
c        Must ensure that 0<Qf<1:
         if(Qf.ge.1.0) then
          delQ = (Qi-1.0)/2.0
          Qf = (Qi+1.0)/2.0
         elseif(Qf.le.0.0) then
          delQ = Qi/2.0
          Qf = Qi/2.0
         endif
c        Check if convergence has been reached:
         if(abs(delQ/Qi).lt.relerr) goto 100
c        Transfer Qf to Qi and continue iterations
         Qi = Qf          
 10   continue

c     Reaching here means we have failed to converge
      write(*, *) 'feneq failed; max iterations exceeded'
      write(*, *) 'bye'
      stop

c     Reaching here means we converged
 100  continue  

c     Return converged value
      feneq = Qf

      return
      end      
