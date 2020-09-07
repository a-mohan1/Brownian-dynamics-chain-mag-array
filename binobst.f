c     This subroutine sorts the obstacles into bins
c     Ref: Allen and Tildesley, Comp. Simulation of Liquids, p. 151
c
c     March 19 2007
c
c     Nobs = # obstacles
c     xobs, yobs = obstacle center x, y coords
c     Nbinsx = # bins along x coord
c     Nbinsy = # bins along y coord
c     Nbins = total # bins (= Nbinsx * Nbinsy)
c     binsizex = bin dimension along x coord
c     binsizey = bin dimension along y coord
c     binhead = array containing one obstacle index of each bin
c     binlist = array containing remaining obstacle indices of each bin 
c     pointed to by the head of each bin, separated by 0s
c
c     Note- origin is at channel lower left hand corner
c     Do the following in the main program (channel dims = L, W): 
c     Nbinsx = int(sqrt(approx #bins * L/W))
c     Nbinsy = int(sqrt(approx #bins * W/L))
c     =>  Nbinsx/Nbinsy = L/W
c     binsizex = L/Nbinsx
c     binsizey = W/Nbinsy
c     => binsizex * binsizey * (Nbinsx*Nbinsy) = L*W
c     Nbins = Nbinsx*Nbinsy (not equal to approx #bins due to truncation by int)

      subroutine binobst(Nobs, xobs, yobs, Nbinsx, Nbinsy, Nbins,
     &                   binsizex, binsizey, binhead, binlist)

      implicit none
      integer Nobs, Nbinsx, Nbinsy, Nbins, binhead(Nbins), binlist(Nobs)
      integer iobs, ibin
      real*8 xobs(Nobs), yobs(Nobs), binsizex, binsizey

c     iobs = counter over obstacles
c     ibin = counter over bins
c     initialize binhead to 0 (binlist is generated in 2nd do loop)
c     iterate over obstacles and generate 2 arrays: binhead and binlist
c     binhead contains the last obstacle index of each bin
c     binlist contains the remaining obstacle indices of each bin
c     pointed to by each head, separated by 0s
c     so the bin indexed I contains the ostacles indexed
c     head(I), list(head(I)),... until list(list..(head(I)..)) becomes 0

      do 100 ibin = 1, Nbins
         binhead(ibin) = 0
 100  continue

      do 200 iobs = 1, Nobs
         ibin = 1 + int(xobs(iobs)/binsizex) + int(yobs(iobs)/binsizey) 
     &               * Nbinsx
         binlist(iobs) = binhead(ibin)
         binhead(ibin) = iobs
 200  continue

      return
      end
 
      
