c     This function calculates one element of Blake's Green's tensor
c     taking as inputs
c     j = direction of velocity at 'destination' bead m
c     k = direction of force at 'source' bead n
c     r = vector directed from bead n to bead m
c     rmag = distance between m and n
c     Rim = vector directed from image of n to m
c     Rimmag = distance between m and image of n
c     h = vertical distance of bead n stokeslet from xy plane
c     a = bead radius

      real*8 function blake(j, k, r, rmag, Rim, Rimmag, h, a)
      implicit none
      integer j, k, kdelta
      real*8 r(3), rmag, Rim(3), Rimmag, h, a
      real*8 temp 

      blake = (3.d0*a)/(4.d0*rmag)*
     & ( kdelta(j, k) + r(j)*r(k)/(rmag*rmag) ) -
     & (3.d0*a)/(4.d0*Rimmag)*( kdelta(j, k) +
     & Rim(j)*Rim(k)/(Rimmag*Rimmag) )
 
      if (k.eq.1) then
         temp = (3.d0*a*h)/(2.d0*Rimmag*Rimmag)*(h*kdelta(j, 1)/Rimmag -
     &   3.d0*h*Rim(j)*Rim(1)/(Rimmag**3) + kdelta(j, 3)*Rim(1)/Rimmag - 
     & kdelta(j, 1)*Rim(3)/Rimmag+3.d0*Rim(j)*Rim(3)*Rim(1)/(Rimmag**3))
      elseif (k.eq.2) then
         temp = (3.d0*a*h)/(2.d0*Rimmag*Rimmag)*(h*kdelta(j, 2)/Rimmag -
     &   3.d0*h*Rim(j)*Rim(2)/(Rimmag**3) + kdelta(j, 3)*Rim(2)/Rimmag - 
     & kdelta(j, 2)*Rim(3)/Rimmag+3.d0*Rim(j)*Rim(3)*Rim(2)/(Rimmag**3))
      else
         temp = -3.d0*a*h/(2.d0*Rimmag*Rimmag)*( h/Rimmag*kdelta(j, 3) - 
     & 3.d0*h*Rim(j)*Rim(3)/(Rimmag**3) - Rim(j)/Rimmag + 
     & 3.d0*Rim(j)*Rim(3)*Rim(3)/(Rimmag**3) )
      endif 

      blake = blake + temp
      
      return
      end
