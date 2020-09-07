c     This subroutine either reads or generates an equilibrated chain config
c     in the obstacle array
c 
c     xc, yc, zc = bead coords to be determined
c     xeq, yeq, zeq = bead coords in equilibrium coil config
c     xn, yn, zn = bead coords determined by pred-corr during equilibration
c     xstart = abs bead 1 x coord upstream of channel 
c     W = channel width
c     nbeads = number of beads
c     sprtype = type of spring
c     Hspr = spring constant
c     Nks = # Kuhn steps per spring
c     v = EV parameter
c     Rousetime = rouse relaxation time of the chain
c     whicRNG = RNG 
c     ranflag = indicates whether RNG has been initialized
c     nseed = seed for RNG
c     deltinit = time step for generating equilibrium coil config
c     delt = time step for array simulation
c     ntstepseq = # time steps to equilibration in array
c     overext = indicates whether a spring has been overstretched
c     isgev = indicates whether Gaussian EV is included
c     istab = indicates whether lookup table for pred-corr exists
c     iseqconfig = indicates whether an equilibrium coil config has
c     already been generated
c     iscoords = indicates whether bead coords file exists
c     readcoords = indicates whether to read bead coords from file
c     nchainfile = unit# of file from which bead coords are to be read
c     tempRg = rad of gyration of eqm coil config; never used
c     a = bead radius for HSEV implementation
c     Nobs = number of obstacles
c     Robs = obstacle radius
c     xobs, yobs = obst center x,y coords
c     binhead, binlist = obstacle bin head and list
c     Nbins = # bins
c     Nbinsx, Nbinsy = number of bins along x, y directions
c     binsizex, binsizey = bin dims along x, y directions
c     Pe = Peclet number = (\mu_0 E N \zeta A)/(k_BT)  
c     uve = unit vector in electric field direction
c     inovercnt = # overext springs as input to pred-corr
c     outovercnt = # overext springs as output from pred-corr
c     estret = chain stretch during equilibration
c     sstep = interval of #steps after which config is sampled
c     maxind = max array size
c     inde = array index for current traj
c     xmax, xmin = max and min chain x coords
c     emsdx = mean sq displacement of center of mass in x dir
c     during equilibration wrt position of initial coil config
c     ecmx = ens-averaged center of mass x position
c     during equilibration wrt position of initial coil config
c     excmi = initial c.m. x coord
c     (reference value during equilibration)
c     xcm = cm x coord at subsequent time points during equilibration

c     Apr 17 07

      subroutine inchainconf(nchainfile, xc, yc, zc, xeq, yeq, zeq, 
     &                       xstart, W, nbeads, sprtype, Hspr, Nks, v,
     &                       Rousetime, whicRNG, ranflag, nseed, delt,
     &                       deltinit, ntstepseq, isgev, overext,
     &                       istab, iseqconfig, iscoords, readcoords, a,
     &                       Nobs, Robs, xobs, yobs, binhead, binlist,
     &                       Nbins, Nbinsx, Nbinsy, binsizex, binsizey,
     &                       inovercnt, outovercnt, Pe, uve, 
     &                       estret, sstep, maxind, inde, emsdx, ecmx)
      
      implicit none
      integer nchainfile, nbeads, sprtype, whicRNG, ranflag, nseed
      integer ntstepseq, Nobs, Nbins, Nbinsx, Nbinsy
      integer binhead(Nbins), binlist(Nobs)
      integer inovercnt, outovercnt, sstep, maxind, inde
      integer istep, ibead
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 xeq(nbeads), yeq(nbeads), zeq(nbeads), xstart, W
      real*8 Hspr, Nks, v, Rousetime, deltinit, delt, a
      real*8 Robs, xobs(Nobs), yobs(Nobs), binsizex, binsizey
      real*8 Pe, uve(3), estret(maxind), emsdx(maxind), ecmx(maxind)
      real*8 xmax, xmin, xcm, excmi
      real*8 tempRg, xn(nbeads), yn(nbeads), zn(nbeads)
      logical isgev, overext, istab, iseqconfig, iscoords, readcoords
 
c     initialize array index to 0 for current traj.
c     index is incremented as and when values are saved to arrays
      inde = 0
    
c     if bead coords file exists and is to be read, then read file,
c     initialize bead coords and return

      if(iscoords.AND.readcoords) then
        do 10 ibead = 1, nbeads
           read(nchainfile, *) xc(ibead), yc(ibead), zc(ibead)
 10     continue
        goto 1000
      endif 
     
c     if not, then generate initial eqm coil and 
c     simulate in array until equilibration
      
      if(readcoords.AND.(.NOT.iscoords)) then
        write(*, *) 'coords file does not exist'
      endif
      write(*, *) 'generating equilibrated chain for current traj.'

c     generate coil; first bead is at the origin
c     (origin = channel lower left hand corner)
      call geninit(xeq, yeq, zeq, tempRg, nbeads, sprtype, Hspr,
     &             whicRNG, ranflag, nseed, deltinit, Rousetime, Nks, v, 
     &             isgev, istab, iseqconfig) 
 
c     transfer xeq, yeq, zeq to xc, yc, zc
c     after shifting bead x coords upstream by subtracting xstart
c     and bead y coords to channel center line by adding W/2

      do 20 ibead = 1, nbeads
          xc(ibead) = xeq(ibead) - xstart
          yc(ibead) = yeq(ibead) + W/2.d0
          zc(ibead) = zeq(ibead)
 20   continue
c     bead 1 is now at (-xstart, W/2.d0, 0)

c     implement HS EV
c     bead-obstacle EV 
      call mhsevobs(xc, yc, nbeads, a, binhead, binlist,
     &              Nbinsx, Nbinsy, Nbins, binsizex,
     &              binsizey, xobs, yobs, Nobs, Robs)

c     calc initial cm x coord of coil config upstream of array
      excmi = xc(1)
      do 30 ibead = 2, nbeads
         excmi = excmi + xc(ibead)
 30   continue
      excmi = excmi/nbeads

c     time stepping to equilibrate in array
c     compute bead positions at istep from their values at istep-1

c     iterate over time steps until equilibration
      do 100 istep = 1, ntstepseq
         
c        predictor-corrector method, accepts overstretched springs
         call pc1sfdEext(xc, yc, zc, xn, yn, zn, nbeads, istep-1, 
     &                   sprtype, Hspr, whicRNG, ranflag, nseed, delt,
     &                   Nks, v, Pe, uve, overext, isgev, istab, 
     &                   inovercnt, outovercnt)
 
c        if overextended input spring remains overstretched after pred-corr,
c        reduce time step and restart program
         if(overext) then
            write(*, *) 'pred-corr returned overstretched spr during eq'
            write(*, *) 'reduce time step and restart; bye'
            stop
         endif  
         
c        implement HS EV with obstacles
         call mhsevobs(xn, yn, nbeads, a, binhead, binlist,
     &                  Nbinsx, Nbinsy, Nbins, binsizex,
     &                  binsizey, xobs, yobs, Nobs, Robs)

c        save data
         if(mod(istep-1, sstep).eq.0) then
          inde = inde + 1
          if(inde.gt.maxind) then
            write(*, *) 'max array sizes exceeded; bye'
            stop
          endif
          xmax = xn(1)
          xmin = xn(1)
          xcm = xn(1)
          do 40 ibead = 2, nbeads
             xcm = xcm + xn(ibead)
             if(xmax.lt.xn(ibead)) then
               xmax = xn(ibead)
             elseif(xmin.gt.xn(ibead)) then
               xmin = xn(ibead)
             endif
 40       continue   
          estret(inde) = estret(inde) + xmax - xmin
          xcm = xcm/nbeads
          emsdx(inde) = emsdx(inde) + (xcm-excmi)*(xcm-excmi)
          ecmx(inde) = ecmx(inde) + (xcm-excmi)
         endif

c        transfer xn, yn, zn to xc, yc, zc for next time step
         do 60 ibead = 1, nbeads
            xc(ibead) = xn(ibead)
            yc(ibead) = yn(ibead)
            zc(ibead) = zn(ibead)
 60      continue   

c        finish iterations over time steps
 100  continue
c     chain config (xc, yc, zc) has been equilibrated in array and 
c     satisfies HS EV constraints
      write(*, *) 'equilibrated chain generated'

 1000 continue
      
      return
      end
