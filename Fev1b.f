c     This subroutine calculates the vector of Gaussian EV
c     forces acting on a specified bead, jbead
c     in a WLC
c     and takes as inputs:
c     xc, yc, zc = bead coordinates at current time step
c     nbeads = number of beads (param)
c     jbead = current bead  
c     v = excluded volume parameter 
c     Nks = number of Kuhn segments per spring
c     The result is stored in vector evforc

      subroutine Fev1b(xc, yc, zc, nbeads, jbead, v, Nks,
     &                 evforc)

      implicit none
      integer nbeads, jbead
      integer lbead, icoord
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 v, Nks, evforc(3)
      real*8 PI, relx, rely, relz, reldsq, const
      parameter (PI = 3.1416)
      
c     lbead, icoord = counters for bead indices and coords
c     (relx, rely, relz) = center-to-center vector from lbead to jbead
c     reldsq = square of the distance between lbead and jbead

c     const = PI*(9/2/PI)**2.5 * Nks**4.5 * v
      const = PI*(4.5d0/PI)**2.5d0 * Nks**4.5 * v     

c     initialize vector evforc
      do 10 icoord = 1, 3
         evforc(icoord) = 0.d0
 10   continue   

c     Loop over beads
c     There will be no contribution from the lbead=jbead term

      do 20 lbead = 1, nbeads
       relx = xc(jbead) - xc(lbead) 
       rely = yc(jbead) - yc(lbead)
       relz = zc(jbead) - zc(lbead)
       reldsq = relx**2 + rely**2 + relz**2
       evforc(1) = evforc(1) + const*exp(-4.5d0*Nks*reldsq)*relx
       evforc(2) = evforc(2) + const*exp(-4.5d0*Nks*reldsq)*rely   
       evforc(3) = evforc(3) + const*exp(-4.5d0*Nks*reldsq)*relz
 20   continue
 
      return
      end



            
           
