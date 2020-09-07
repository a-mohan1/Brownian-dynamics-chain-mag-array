c     This subroutine makes a look-up table for calculating
c     Qmag from rmag
c     for FENE and WLC springs
c     Inputs:
c     numel = number of nonzero elements in the table
c     tabint = interval between successive values of rmag 
c     sprtype = spring force law
c      3 = FENE
c      4 = WLS
c     Hspr = spring constant 
c     delt = time step
c     The table is stored in array table

      subroutine maketab(numel, tabint, sprtype, Hspr, delt, 
     &                   table)

      implicit none
      integer numel, sprtype, i
      real*8 tabint, Hspr, delt, table(0:numel), rmag
      real*8 wlcq, feneq
 
c     i = table index
c     rmag = magnitude of rhs vector      
c     Qmag = 0 when rmag = 0
      table(0) = 0.d0
      
c     table contains Qmag values corresponding to numel nonzero values of rmag 
c     spaced at intervals of tabint
c     Function wlcq supplies Qmag for a given rmag for a WLS
c     and feneq supplies the same for a FENE spring

      if(sprtype.eq.4) then
       do 10 i = 1, numel
          rmag = dble(i)*tabint
          table(i) = wlcq(rmag, Hspr, delt)
 10    continue
      elseif(sprtype.eq.3) then
       do 20 i = 1, numel
          rmag = dble(i)*tabint
          table(i) = feneq(rmag, Hspr, delt)
 20    continue
      else
       write(*, *) 'unknown spring'
       write(*, *) 'bye'
       stop
      endif

      return
      end
