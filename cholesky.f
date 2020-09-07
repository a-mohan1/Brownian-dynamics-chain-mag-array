c     This subroutine implements the cholesky decomposition of matrix D
c     taking as inputs
c     D = matrix to be decomposed (must be SPD)
c     ndim = dimension of D
c     istep = current time step
c     and stores the result in lower triangular matrix B

      subroutine cholesky(D, ndim, istep, B)
      implicit none
      integer ndim, istep
      integer i, j, k
      real*8 D(ndim, ndim), B(ndim, ndim)
      real*8 sum
 
      do 10 i = 1, ndim
       do 20 j = i, ndim
        sum = D(i, j) 
        do 30 k = i-1, 1, -1    
         sum = sum - B(i, k)*B(j, k)
 30     continue
        if (i.eq.j) then
         if(sum.le.0) then
          write(*, *) 'non-PD diffusion tensor, time = ', istep
          write(*, *) 'bye'
          stop
         else
          B(i, i) = sqrt(sum)
         endif
        else
         B(j, i) = sum/B(i, i)
         B(i, j) = 0
        endif
 20    continue
 10   continue 

      return 
      end
