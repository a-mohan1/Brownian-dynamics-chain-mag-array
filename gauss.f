c     This function returns a standard normal deviate 
c     by the Box-Muller transformation, 
c     using function rannum to generate uniform deviates
c     See Numerical Recipes in Fortran by Press et al., pp. 279-280

      real*8 function gauss(whicRNG, ranflag, jseed)

      implicit none
      integer whicRNG, ranflag, jseed
      real*8 rannum
      integer iset
      real fac, gset, rsq, v1, v2
      save iset, gset
      data iset/0/
      if(iset.eq.0) then
 1      v1 = 2.0*rannum(whicRNG, ranflag, jseed)-1.0
        v2 = 2.0*rannum(whicRNG, ranflag, jseed)-1.0
        rsq = v1**2 + v2**2
        if(rsq.ge.1.0.OR.rsq.eq.0.0) goto 1
        fac = sqrt(-2.0*log(rsq)/rsq)        
        gset = v1*fac
        gauss = v2*fac
        iset = 1
      else
        gauss = gset
        iset = 0
      endif

      return
      end

      

      

