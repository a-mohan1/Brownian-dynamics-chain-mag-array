c     This subroutine calculates the spring force due to one spring
c     taking as inputs:
c     the spring vector (sprvec),
c     the type of spring (sprtype)
c       0 = no spring
c       1 = Hookean
c       2 = Fraenkel
c       3 = FENE
c       4 = WLC (Marko-Siggia)
c     and the non-dimensional spring constant (Hspr)
c     The spring force laws are taken from DPL v.2, p. 21
c     and the MS law is from the Jendrejack papers/ Larson's review
c     The result is stored in the array force
c     If a FENE or a WLC spring is used and the spring is overstretched,
c     force will be set to a zero vector and 
c     overext will be set to .TRUE.

      subroutine calcFspr(sprtype, sprvec, Hspr, force, overext)
      
      implicit none
      integer sprtype, i
      real*8 sprvec(3), Hspr, force(3)
      real*8 Qmag
      logical overext

c     initialize overext to .FALSE.
      overext = .FALSE.
 
c     spring end-to-end distance
      Qmag = sqrt(sprvec(1)**2 + sprvec(2)**2 + sprvec(3)**2)

c     identify spring type and invoke the appropriate force law
 
      if (sprtype.eq.0) then
c     no spring, no force
        do 10 i = 1, 3 
           force(i) = 0.d0
 10     continue
      elseif (sprtype.eq.1) then
c     linear Hookean  
        do 11 i = 1, 3 
           force(i) = Hspr*sprvec(i)
 11     continue
      elseif (sprtype.eq.2) then
c     Fraenkel
        do 12 i = 1, 3 
           force(i) = Hspr*(1-1/Qmag)*sprvec(i)
 12     continue
      elseif (sprtype.eq.3) then
c     FENE
        if(Qmag.ge.1.d0) then
          write(*, *) 'overstretched FENE spring encountered'
          overext = .TRUE.
          do 13 i = 1, 3
             force(i) = 0.d0
 13       continue
        else
          do 14 i = 1, 3   
             force(i) = Hspr*sprvec(i)/(1-Qmag**2)
 14       continue
        endif   
      elseif (sprtype.eq.4) then
c     WLC
        if(Qmag.ge.1.d0) then
          write(*, *) 'overstretched WLS encountered'
          overext = .TRUE.
          do 15 i = 1, 3
             force(i) = 0.d0
 15       continue
        else
          do 16 i = 1, 3   
             force(i) = 2.d0*Hspr/(3.d0*Qmag)*(1/(4.d0*(1-Qmag)**2) -
     &       1/4.d0 + Qmag)*sprvec(i)
 16       continue
        endif      
      else
c     spring force law not specified
         write(*,*) 'unknown spring'
         write(*,*) 'bye'
         stop
      endif

      return
      end
      
