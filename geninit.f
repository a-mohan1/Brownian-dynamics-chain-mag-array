c     This subroutine generates initial configurations
c     by generating a Gaussian coil and evolving it for 5 rel times
c     when called for the first time
c     and subsequently evolving the previous init config for 1 rel time 

c     Variables:
c     xeq, yeq, zeq = eqm bead coords to be used as init config in program
c     Rg = eqm radius of gyration
c     nbeads = number of beads 
c     sprtype = spring force law (no springs=0, Hookean=1, Fraenkel=2, 
c                                 FENE=3, WLS=4)
c     Hspr = spring constant = 3*Nks/lambda
c     whicRNG = which RNG to use to generate uniform deviates (Knuth = 1)
c     ranflag = 0 at first call to RNG, 1 otherwise
c     nseed = seed for RNG (negative integer for whicRNG = 1)
c     delt = time step
c     relt = Rouse relaxation time
c     Nks = number of Kuhn steps per spring
c     v = EV parameter, used if isgev is .TRUE.
c     isgev is .TRUE. if a Gaussian EV force is included
c     istab indicates whether or not the look-up table for FENE and WLC springs 
c     has been created and saved
c     iseqconfig indicates whether an eqm config has been previously generated 
c     by calling subroutine geninit
c     initialized to .FALSE. in program; 
c     changed to .TRUE. during first call to geninit

c     all eqs and vars are non-dimensional

      subroutine geninit(xeq, yeq, zeq, Rg, nbeads, sprtype, Hspr,
     &                  whicRNG, ranflag, nseed, delt, relt, 
     &                  Nks, v, isgev, istab, iseqconfig)

c     ntsteps = number of time steps
c     kappa = 3-by-3 transposed vel gradient tensor (non-dimensionalized by 
c     convective time scale) times Peclet number, not relevant for eqm calc.,
c     set to 0
c     tether = .FALSE. (pc1sfd.f does not work correctly for tethered chains)
c     overext is set to .TRUE. if subroutine calcFspr encounters 
c     an overextended FENE spring or WLS
c     isflow is .FALSE. since there is no imposed linear fluid flow (eqm calc.)
c     ibead, istep = counters for beads/ time steps
c     i, j = counters for coords
c     stdev = std dev of spring end-to-end x, y, z coords
c     xcm, ycm, zcm = center of mass coords

      implicit none
      integer nbeads, ntsteps, sprtype, whicRNG, ranflag, nseed
      integer ibead, istep, i, j
      real*8 xeq(nbeads), yeq(nbeads), zeq(nbeads)
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 xn(nbeads), yn(nbeads), zn(nbeads)
      real*8 xcm, ycm, zcm
      real*8 Hspr, delt, relt, Nks, v
      real*8 kappa(3, 3)
      real*8 stdev, gauss, Rg
      logical isgev, istab, iseqconfig
      logical tether, overext, isflow

      parameter (tether = .FALSE., isflow = .FALSE.)
      
c     initialize vars
      overext = .FALSE.
      Rg = 0.d0
c     define vel gradient transpose for no linear flow       
      do 20 i = 1, 3
       do 10 j = 1, 3
          kappa(i, j) = 0.d0
 10    continue
 20   continue

c     check if this is the first call to subroutine geninit; 
c     if not, initial eqm config exists already
c     if init eqm config does not exist, generate gaussian coil
c     and evolve for 5 rel times
c     if eqm config exists then evolve it for 1 rel time

      if(.NOT.iseqconfig) then 
         iseqconfig = .TRUE.
c        run for five rel times
         ntsteps = 5*relt/delt 
c        standard dev. of gaussian spring x, y, z end-to-end coords
         stdev = 1.d0/sqrt(3.d0*Nks)
c        generate gaussian coil as starting config 
c        function gauss returns a std normal deviate
         xc(1) = 0.d0
         yc(1) = 0.d0
         zc(1) = 0.d0
         do 30 ibead = 2, nbeads
            xc(ibead) = xc(ibead-1) +   
     &      gauss(whicRNG, ranflag, nseed)*stdev
            yc(ibead) = yc(ibead-1) + 
     &      gauss(whicRNG, ranflag, nseed)*stdev
            zc(ibead) = zc(ibead-1) + 
     &      gauss(whicRNG, ranflag, nseed)*stdev
 30      continue
      else
c        run for 1 rel time
         ntsteps = relt/delt
c        initialize starting config
         do 50 ibead = 1, nbeads 
            xc(ibead) = xeq(ibead)
            yc(ibead) = yeq(ibead)
            zc(ibead) = zeq(ibead)
 50      continue
      endif   
       
c     run FD simulation for given init config  
c     iterate over time steps 
c     calc bead coords xn, yn, zn at istep from xc, yc, zc at istep-1
      do 100 istep = 1, ntsteps
c           time-stepping
            call pc1sfd(xc, yc, zc, xn, yn, zn, nbeads, istep-1, 
     &      sprtype, Hspr, whicRNG, ranflag, nseed, delt, Nks, v,
     &      kappa, tether, overext, isgev, isflow, istab)

c           check for overextended springs
            if(overext) then
              write(*, *) 'overextended spring encountered; bye'
              stop
            endif

c           transfer xn, yn, zn to xc, yc, zc for next time step
            do 60 ibead = 1, nbeads
             xc(ibead) = xn(ibead)
             yc(ibead) = yn(ibead)
             zc(ibead) = zn(ibead)
 60         continue

c           finish iterations over time steps
 100  continue

c     save config at final time step 
      do 120 ibead = 1, nbeads
            xeq(ibead) = xn(ibead)
            yeq(ibead) = yn(ibead)
            zeq(ibead) = zn(ibead) 
 120  continue

c     calculate cm
      xcm = xeq(1)
      ycm = yeq(1)
      zcm = zeq(1)
      do 140 ibead = 2, nbeads
         xcm = xcm + xeq(ibead)
         ycm = ycm + yeq(ibead)
         zcm = zcm + zeq(ibead)
140   continue
      xcm = xcm/nbeads
      ycm = ycm/nbeads
      zcm = zcm/nbeads 
            
c     calc Rg   
      do 160 ibead = 1, nbeads
               Rg = Rg + (xeq(ibead)-xcm)**2 + (yeq(ibead)-ycm)**2 + 
     &                   (zeq(ibead)-zcm)**2
 160  continue
      Rg = sqrt(Rg/nbeads)
      write(*, *) 'Rg = ', Rg
      
c     shift origin to the location of the 1st bead
      do 200 ibead = 2, nbeads
         xeq(ibead) = xeq(ibead) - xeq(1)
         yeq(ibead) = yeq(ibead) - yeq(1)
         zeq(ibead) = zeq(ibead) - zeq(1)
 200  continue   
      xeq(1) = 0.d0
      yeq(1) = 0.d0
      zeq(1) = 0.d0

      return 
      end
