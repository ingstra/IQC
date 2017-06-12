
MODULE hermite
 integer, parameter :: dp = selected_real_kind(15, 307)

CONTAINS
 RECURSIVE FUNCTION HermitePoly(n) RESULT(hp)

!
! More info
! http://mathworld.wolfram.com/HermitePolynomial.html
!

      REAL(dp) hp(n+1)
      REAL(dp) hp1(n+1),hp2(n+1)
      IF(n .EQ. 0) THEN

         hp(1)=1.0
         RETURN
      END IF
      IF(n .EQ. 1) THEN

         hp(1)=2.0
         hp(2)=0.0
      ELSE

       hp1(1:n+1)=0.0
       hp1(1:n)=2.0*HermitePoly(n-1)

       hp2(1:n+1)=0.0
       hp2(3:)=2.0*(n-1)*HermitePoly(n-2)

       hp=hp1-hp2

      END IF

      END FUNCTION

      FUNCTION evalHermitePoly(ix,n) RESULT(y)
      INTEGER n,ip
      REAL(dp) ix(:),y(size(ix)),h(n+1)

      k=size(ix)

      h=HermitePoly(n)

      y(1:k)=h(n+1)
      ip=1
      DO i=n,1,-1
        DO j=1,k
         y(j)=y(j)+h(i)*ix(j)**ip
        END DO

        ip=ip+1
      END DO
      END FUNCTION

END MODULE hermite

