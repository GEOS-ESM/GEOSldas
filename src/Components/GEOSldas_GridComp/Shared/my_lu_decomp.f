c
c  routines adapted from Numerical Recipes
c  
c  reichle, 20 Apr 01
c  
c  
c  ***********************************************************************
c  
c  ------- documentation from C function ---------
c
c  LU decomposition after Numerical Recipes 
c  
c  Given a matrix a[1..n][1..n], this routine replaces it by the
c  LU decomposition of a rowwise permutation of itself. a and n are
c  inputs. a is output, arranged as in equation (2.3.14) (see book).
c  indx[1..n] is an output vector that records the row permutation
c  effected by the partial pivoting; This routine is used in combination
c  with lubksb to solve linear equations or invert a matrix.
c  
c  edited 20 Apr 01, reichle
c
c  - eliminated d (for computation of determinant)
c  - eliminated np (can use automatic arrays nowadays...)
c
c  edited 03 Jun 02, reichle
c  
c  - eliminated parameter NMAX, using dynamic allocation
c
c  -----------------------------------------------------
c  
      SUBROUTINE ludcmp(a,n,indx)
      INTEGER n,indx(n)
      REAL a(n,n),TINY
      PARAMETER (TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(n)
      do 12 i=1,n
         aamax=0.
         do 11 j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue
         if (aamax.eq.0.) pause 'singular matrix in ludcmp'
         vv(i)=1./aamax
 12   continue
      do 19 j=1,n
         do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
 13         continue
            a(i,j)=sum
 14      continue
         aamax=0.
         do 16 i=j,n
            sum=a(i,j)
           do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
 15        continue
           a(i,j)=sum
           dum=vv(i)*abs(sum)
           if (dum.ge.aamax) then
              imax=i
              aamax=dum
           endif
 16     continue
        if (j.ne.imax)then
           do 17 k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
 17        continue
           vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
           dum=1./a(j,j)
           do 18 i=j+1,n
              a(i,j)=a(i,j)*dum
 18        continue
        endif
 19   continue
      return
      END      

  
c  ***********************************************************************
c  
c  ------- documentation from C function ---------
c
c  LU backsubstitution after Numerical Recipes 
c 
c  Solves the set of n linear equations A X = B. Here a[1..n][1..n] is
c  input, not as the matrix A but rather as its LU decomposition, 
c  determined by the routine ludcmp. indx[1..n] is input as the permutation
c  vector returned by ludcmp. b[1..n] is input as the right-hand side
c  vector B, and returns the solution vector X. a,n, and indx ar not 
c  modified by this routine and can be left in place for successive calls
c  with different right-hand sides b. This routine takes into account
c  the possibility that b will begin with many zero elements, so it is
c  efficient for use in matrix inversion.
c  
c  edited 20 Apr 01, reichle
c  
c  - eliminated np (can use automatic arrays nowadays...)
c  
c  -----------------------------------------------------
c  
      SUBROUTINE lubksb(a,n,indx,b)
      INTEGER n,indx(n)
      REAL a(n,n),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 11         continue
         else if (sum.ne.0.) then
            ii=i
         endif
         b(i)=sum
 12   continue
      do 14 i=n,1,-1
         sum=b(i)
         do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
 13      continue
         b(i)=sum/a(i,i)
 14   continue
      return
      END
c
c
c ********** EOF ******************************************************

