c     This program simulates the motion of a chain
c     in an obstacle array (hexagonal or square)
c     Initially, the chain is in a coil config, with the first bead at (0,0,0)
c     Free draining conditions are assumed
c
c     Feb 27 07

      program obstar

c     xc, yc, zc, xn, yn, zn = bead coords at successive time steps
c     nbeads = number of beads
c     ntsteps = number of time steps
c     eqmdist = distance covered by chain before program is terminated
c               (if ntsteps not specified)
c     a = bead radius
c     sprtype = type of spring (WLS = 4)
c     nens = ensemble size
c     Hspr = spring constant
c     whicRNG = RNG (Knuth = 1)
c     nseed = seed for RNG (-ve integer for whicRNG=1)
c     ranflag = indicates whether or not RNG has been initialized; initialized to 0
c     delt = time step size
c     Nks = number of Kuhn segments per spring
c     lambda = ratio of effective to true persistence length
c     v = EV param
c     Pe = Peclet number (\mu_0 E N \zeta A)/(k_B T)
c     uve = unit vector in field direction, x: (1,0,0)
c     Robs = obstacle radius
c     xstart = x coord of plane at which lattice starts
c     lspac = lattice spacing
c     lattype = lattice pattern (hexagonal = 0, square = 1)
c     hchann = channel height (bottom: z=-hchann/2, top: z=hchann/2)
c              (if the chain is confined)
c     overext indicates whether a spring has been overstretched
c     isgev = .TRUE. if gaussian EV is included
c     istab indicates whether or not a lookup table exists 
c     (reqd by predictor-corrector scheme)
c     iseqconfig indicates whether or not an initial eqm config
c     has been generated by subroutine geninit
c     isateq indicates whether or not eqmdist has been covered
c     istep = counter for time stepping
c     ibead = counter for beads
c     iens = counter over trajectories in ensemble
c     sstep = interval of number of steps after which config is sampled
c     xeq, yeq, zeq = initial eqm bead coords
c     Rousetime = Rouse rel time
c     tempRg = eqm rad of gyration returned by geninit, never used
c     xmin, xmax = min and max chain coords in field direction
c     stret = ensemble-averaged stretch in field direction
c     isamp = sampling index
c     nsamp = sample size
c     maxnsamp = maximum allowable sample size
c     inovercnt = # times pred corr method gets overextended spring
c     outovercnt = # times pred corr method returns overextended spring
c                  (never used, since program is stopped if this happens)

c     all eqs and variables are nondimensional
c     nondimensional max spring length = 1.0

      real etime, elapsed(2)

      integer nbeads, ntsteps, sprtype, whicRNG, nseed, ranflag
      integer nens, lattype
      integer istep, ibead, iens, sstep, isamp, nsamp, maxnsamp
      integer inovercnt, outovercnt

      parameter (nbeads = 38)
      parameter (whicRNG = 1)
      parameter (nseed = -17)
      parameter (sprtype = 4)
      parameter (nens = 100)
      parameter (lattype = 0)
c      parameter (ntsteps = 1.d5)
      parameter (sstep = 50)
      parameter (maxnsamp = 20000)

      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 xn(nbeads), yn(nbeads), zn(nbeads)
      real*8 xeq(nbeads), yeq(nbeads), zeq(nbeads)
      real*8 a, Hspr, delt, Nks, lambda, v, Pe, uve(3)
      real*8 Robs, xstart, lspac, hchann, eqmdist
      real*8 Rousetime, tempRg
      real*8 xmin, xmax, stret(maxnsamp)
       
      parameter (a = 0.d0)
      parameter (delt = 1.d-3)
      parameter (Nks = 5.23d0)
      parameter (lambda = 1.91d0)
      parameter (Hspr = 3.d0*Nks/lambda)
      parameter (v = 2.35d-3)
      parameter (Pe = 5.d0)
      parameter (Robs = 0.903d0)
      parameter (xstart = 2.d0)
      parameter (lspac = 5.42d0)
c      parameter (hchann = 3.61d0)
      parameter (eqmdist = 50.d0*lspac)
      parameter (Rousetime = 8.91d0)

      logical overext, isgev, istab, iseqconfig
      logical isateq

      parameter (isgev = .TRUE.)

      write(*, *) 'welcome'

c     initialization
      ranflag = 0
      overext = .FALSE.
      istab = .FALSE.
      iseqconfig = .FALSE.
      nsamp = maxnsamp
      do 10 isamp = 1, maxnsamp
         stret(isamp) = 0.d0
 10   continue
      inovercnt = 0
      outovercnt = 0

c     unit vector in electric field direction
      uve(1) = 1.d0
      uve(2) = 0.d0
      uve(3) = 0.d0

c     open file for saving xyz coords
      open(unit = 5, file = 'xyz.dat', status = 'NEW')

c     iterate over trajectories in ensemble
      do 1000 iens = 1, nens
       write(*, *) 'trajectory ', iens

c      initialization for current trajectory
       isateq = .FALSE.
       istep = 0
       isamp = 0
      
c      generate initial eqm config for current trajectory
c      first bead is at the origin
       call geninit(xeq, yeq, zeq, tempRg, nbeads, sprtype, Hspr,
     &             whicRNG, ranflag, nseed, delt, Rousetime, Nks, v, 
     &             isgev, istab, iseqconfig) 

       do 300 ibead = 1, nbeads
         xc(ibead) = xeq(ibead)
         yc(ibead) = yeq(ibead)
         zc(ibead) = zeq(ibead)
 300   continue

c      implement HS EV
c      bead-wall EV (if confined)
c       call hsevwall(zc, nbeads, a, hchann)
c      bead-obstacle EV
       if(lattype.eq.0) then
         write(*, *) 'hexagonal lattice'
         call hsevhexarray(xc, yc, nbeads, a, Robs, lspac, xstart)
       elseif(lattype.eq.1) then
         write(*, *) 'square lattice'
         call hsevsqarray(xc, yc, nbeads, a, Robs, lspac, xstart)
       else
         write(*, *) 'undefined lattice; bye'
         stop
       endif 

c      time stepping
c      compute bead positions at istep from their values at istep-1

 400   continue
       if(isamp.ge.maxnsamp) then
          write(*,*) 'maximum sample size reached; bye'
          stop
       endif
       istep = istep + 1 
       call pc1sfdEext(xc, yc, zc, xn, yn, zn, nbeads, istep-1, sprtype,
     &        Hspr, whicRNG, ranflag, nseed, delt, Nks, v, Pe, uve, 
     &        overext, isgev, istab, inovercnt, outovercnt)
 
       if(overext) then
            write(*, *) 'overstretched spring; bye'
            stop
       endif  
         
c      implement HS EV
c       call hsevwall(zn, nbeads, a, hchann)
       if(lattype.eq.0) then
            call hsevhexarray(xn, yn, nbeads, a, Robs, lspac, xstart)
       elseif(lattype.eq.1) then
            call hsevsqarray(xn, yn, nbeads, a, Robs, lspac, xstart)
       else
            write(*, *) 'undefined lattice; bye'
            stop
       endif 

c      calculate max and min x coords
       xmin = xn(1)
       xmax = xn(1)
       do 800 ibead = 2, nbeads
          if(xn(ibead).lt.xmin) then
             xmin = xn(ibead)
          elseif(xn(ibead).gt.xmax) then
             xmax = xn(ibead)
          endif
 800   continue       
       
c      sampling
       if(mod(istep-1, sstep).eq.0) then
         isamp = isamp + 1
c        save stretch
         stret(isamp) = stret(isamp) + xmax - xmin           
             
c        save coords to xyz file for any one trajectory, say iens=nens
         if(iens.eq.nens) then
           write(5, 2003) nbeads
           do 850 ibead = 1, nbeads
              write(5, 2002) xn(ibead), yn(ibead), zn(ibead)
 850       continue   
         endif
       endif    

c      transfer xn, yn, zn to xc, yc, zc for next time step
       do 900 ibead = 1, nbeads
          xc(ibead) = xn(ibead)
          yc(ibead) = yn(ibead)
          zc(ibead) = zn(ibead)
 900   continue   

c      check if equilibration dist has been covered
c      if not, continue iterations
       if(xmin.ge.eqmdist) then
          isateq = .TRUE.
       endif

       if(.NOT.isateq) then
          goto 400
       endif

c      identify smallest sample size among trajectories
       if(nsamp.gt.isamp) then
         nsamp = isamp
         ntsteps = istep
       endif
c      finish iterations over trajectories
 1000 continue

      close(unit = 5, status = 'KEEP')

c     sample size and #steps correspond to smallest 
c     among all trajectories
c     ensemble-averaging
      do 1020 isamp = 1, nsamp
         stret(isamp) = stret(isamp)/nens
 1020 continue   
      
c     write output to files
      open(unit = 10, file = 'stret.dat', status = 'NEW')
      write(10, 2001) (stret(isamp), isamp = 1, nsamp)
      close(unit = 10, status = 'KEEP')

 2001 format(1X, 1e12.6)
 2002 format(3(1X, 1e12.6))
 2003 format(i4)

      write(*, *) '# input overext springs = ', inovercnt
      write(*, *) '# output overext springs = ', outovercnt
      write(*, *) 'number of time steps = ', ntsteps
      write(*, *) 'sample size = ', nsamp
      write(*, *) 'output written. bye'

      write(*, *) 'total time taken = ', etime(elapsed)
     
      stop
      end


    
