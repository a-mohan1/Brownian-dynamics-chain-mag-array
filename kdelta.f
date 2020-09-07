c     This function implements the Kronecker delta

      integer function kdelta(j, k)
      implicit none 
      integer j, k
      if (j.eq.k) then
         kdelta = 1
      else
         kdelta = 0
      endif

      return
      end
   
