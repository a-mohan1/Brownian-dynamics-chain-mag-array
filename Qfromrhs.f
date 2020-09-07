c     This subroutine calculates the spring vector from the right hand side 
c     vector at the end of the 1st or 2nd corrector steps 
c     and is called by subroutine pc1s
c     Inputs:
c     rhs = right hand side vector
c     sprtype = type of spring 
c      0 = no springs
c      1 = Hookean
c      2 = Fraenkel
c      3 = FENE
c      4 = WLC
c     Hspr = spring constant
c     delt = time step 
c     istab = .FALSE. if a look-up table does not already exist,
c     and .TRUE. otherwise; initialized to .FALSE. in program,
c     used for FENE and WLC springs only,
c     modified to .TRUE. during first call to Qfromrhs
c     The resulting spring vector is stored in vector Qvect

      subroutine Qfromrhs(rhs, sprtype, Hspr, delt, Qvect, istab)
      implicit none
      integer sprtype, icoord, ind, numel
      real*8 rhs(3), Hspr, delt, Qvect(3)
      real*8  tabint, rmag, Qmag
      parameter (numel = 1.d7, tabint = 1.d-5)
      real*8 table(0:numel) 
      logical istab
      save table
           
c     icoord = coord index from 1..3
c     ind =  index of element in look-up table
c     numel = number of nonzero elements in look-up table for FENE or WLC springs
c     tabint = intervals by which rmag values in look-up table are separated,
c     i.e., there are 1e7 nonzero elements in table at intervals of 1e-5;
c     elements range from 0 to 100
c     rmag = magnitude of rhs vector

      rmag = sqrt(rhs(1)**2 + rhs(2)**2 + rhs(3)**2)

c     Qmag = magnitude of spring vector Qvect
c     Note that Qvect and rhs are parallel, 
c     hence, Qvect(icoord) = rhs(icoord)*Qmag/rmag
c     table = look-up table with Qmag corresponding to 
c     numel nonzero values of rmag at intervals of tabint

      if(sprtype.eq.0) then 
c      no springs
       Qmag = rmag 
       do 10 icoord = 1, 3
        Qvect(icoord) = rhs(icoord)
 10    continue
      elseif(sprtype.eq.1) then 
c      Hookean
       Qmag = rmag/(1.0+2.0*delt*Hspr)
       do 20 icoord = 1, 3 
        Qvect(icoord) = rhs(icoord)/(1.0+2.0*delt*Hspr)
 20    continue
      elseif(sprtype.eq.2) then
c      Fraenkel
       Qmag = (rmag+2.0*delt*Hspr)/(1.0+2.0*delt*Hspr)
       do 30 icoord = 1, 3  
        Qvect(icoord) = rhs(icoord)*Qmag/rmag
 30    continue
      elseif(sprtype.eq.3) then
c      FENE
c      check if look-up table exists, and if not, create it now
c      and set istab to .TRUE. once the table is created
       if(.NOT.istab) then
        call maketab(numel, tabint, sprtype, Hspr, delt, table)
        istab = .TRUE.
       endif
c      Extract Qmag from the table
c      The index of the largest element of the table less than rmag is [rmag/tabint]
       ind = int(rmag/tabint)
c      if rmag lies within the table, linearly interpolate;
c      otherwise, use the asymp. expansion for large rmag, correct to 1/rmag**2
       if(ind.lt.numel) then
        Qmag = table(ind) + (rmag/tabint - dble(ind))*
     &                     (table(ind+1)-table(ind))
       else
        Qmag = 1.0 - delt*Hspr/rmag 
     &         - 1/(rmag**2)*delt*Hspr*(1.0-0.5*delt*Hspr)
       endif
c      Calculate the elements of Qvect
       do 35 icoord = 1, 3
        Qvect(icoord) = rhs(icoord)*Qmag/rmag
 35    continue 
      elseif(sprtype.eq.4) then
c      WLS     
c      check if look-up table exists, and if not, create it now
c      and set istab to .TRUE. once the table is created
       if(.NOT.istab) then
        call maketab(numel, tabint, sprtype, Hspr, delt, table)
        istab = .TRUE.
       endif
c      Extract Qmag from the table
c      The index of the largest element of the table less than rmag is [rmag/tabint]
       ind = int(rmag/tabint)
c      if rmag lies within the table, linearly interpolate;
c      otherwise, use the asymp. expansion for large rmag, correct to 1/sqrt(rmag)
c      assuming that delt*Hspr/(3*rmag) < 1
       if(ind.lt.numel) then
        Qmag = table(ind) + (rmag/tabint - dble(ind))*
     &                     (table(ind+1)-table(ind))
       else
        Qmag = 1.0 - sqrt(delt*Hspr/(3.0*rmag))
       endif
c      Calculate the elements of Qvect
       do 40 icoord = 1, 3
        Qvect(icoord) = rhs(icoord)*Qmag/rmag
 40    continue
      else
c      unknown spring
       write(*, *) 'unknown spring'
       write(*, *) 'bye'
       stop
      endif

      return
      end
