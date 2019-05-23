
MODULE cubic
    IMPLICIT NONE
    real (KIND = 8), parameter :: U = 0.20d0, valin = 0.85d0, VV = -8.0d0, pi = 4.0d0*datan(1.0d0)
    integer (KIND = 8), parameter :: npt=300
    real(KIND = 8) :: delta0 = 0.0d0,VTI,ezero,t1=1.0d0,t2=0.4d0,ek,deltak,gammak,tot,ptreal
    real(KIND = 8) :: g3,g4,g6,un,valn,V
    integer (KIND = 8) :: check,pontos = 120,eita
END MODULE cubic

MODULE PolynomialRoots
! ---------------------------------------------------------------------------
! PURPOSE - Solve for the roots of a polynomial equation with real
!   coefficients, up to quartic order. Retrun a code indicating the nature
!   of the roots found.

! AUTHORS - Alfred H. Morris, Naval Surface Weapons Center, Dahlgren,VA
!           William L. Davis, Naval Surface Weapons Center, Dahlgren,VA
!           Alan Miller,  CSIRO Mathematical & Information Sciences
!                         CLAYTON, VICTORIA, AUSTRALIA 3169
!                         http://www.mel.dms.csiro.au/~alan
!           Ralph L. Carmichael, Public Domain Aeronautical Software
!                         http://www.pdas.com
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!    ??    1.0 AHM&WLD Original coding                     
! 27Feb97  1.1   AM    Converted to be compatible with ELF90
! 12Jul98  1.2   RLC   Module format; numerous style changes
!  4Jan99  1.3   RLC   Made the tests for zero constant term exactly zero


  IMPLICIT NONE

  INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0), DP=KIND(1.0D0)
  REAL(DP),PARAMETER,PRIVATE:: ZERO=0.0D0, FOURTH=0.25D0, HALF=0.5D0
  REAL(DP),PARAMETER,PRIVATE:: ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, FOUR=4.0D0
  COMPLEX(DP),PARAMETER,PRIVATE:: CZERO=(0.D0,0.D0)

  REAL(DP),PARAMETER,PRIVATE:: EPS=EPSILON(ONE)

  CHARACTER(LEN=*),PARAMETER,PUBLIC:: POLYROOTS_VERSION= "1.3 (4 Jan 1999)"
  INTEGER,PRIVATE:: outputCode
!    =0 degenerate equation
!    =1 one real root
!    =21 two identical real roots
!    =22 two distinct real roots
!    =23 two complex roots
!    =31 multiple real roots
!    =32 one real and two complex roots
!    =33 three distinct real roots
!    =41
!    =42 two real and two complex roots
!    =43
!    =44 four complex roots

  PRIVATE:: CubeRoot
  PUBLIC:: LinearRoot
  PRIVATE:: OneLargeTwoSmall
  PUBLIC:: QuadraticRoots
  PUBLIC:: CubicRoots
  PUBLIC:: QuarticRoots
  PUBLIC:: SolvePolynomial
!----------------------------------------------------------------------------

  INTERFACE Swap
    MODULE PROCEDURE SwapDouble, SwapSingle
  END INTERFACE

CONTAINS

!+
FUNCTION CubeRoot(x) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the Cube Root of a REAL(DP) number. If the argument is
!   negative, then the cube root is also negative.

  REAL(DP),INTENT(IN) :: x
  REAL(DP):: f
!----------------------------------------------------------------------------
  IF (x < ZERO) THEN
    f=-EXP(LOG(-x)/THREE)
  ELSE IF (x > ZERO) THEN
    f=EXP(LOG(x)/THREE)
  ELSE
    f=ZERO
  END IF
  RETURN
END Function CubeRoot   ! ---------------------------------------------------

!+
SUBROUTINE LinearRoot(a, z)
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z 
!     AND STORES THE RESULTS IN Z. It is assumed that a(2) is non-zero.
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  REAL(DP),INTENT(OUT):: z
!----------------------------------------------------------------------------
  IF (a(2)==0.0) THEN
    z=ZERO
  ELSE
    z=-a(1)/a(2)
  END IF
  RETURN
END Subroutine LinearRoot   ! -----------------------------------------------

!+
SUBROUTINE OneLargeTwoSmall(a1,a2,a4,w, z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the roots of a cubic when one root, w, is known to be
!   much larger in magnitude than the other two

  REAL(DP),INTENT(IN):: a1,a2,a4
  REAL(DP),INTENT(IN):: w
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


  REAL(DP),DIMENSION(3):: aq
!----------------------------------------------------------------------------
  aq(1)=a1
  aq(2)=a2+a1/w
  aq(3)=-a4*w
  CALL QuadraticRoots(aq, z)
  z(3)=CMPLX(w,ZERO,DP)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3)=z(2)
  z(2)=z(1)
  z(1)=CMPLX(w,ZERO,DP)
  RETURN
END Subroutine OneLargeTwoSmall   ! -----------------------------------------

!+
SUBROUTINE QuadraticRoots(a, z)
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.

  REAL(DP),INTENT(IN),DIMENSION(:):: a
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


  REAL(DP):: d, r, w, x, y
!----------------------------------------------------------------------------
  IF(a(1)==0.0) THEN     ! EPS is a global module constant (private)
    z(1) = CZERO               ! one root is obviously zero
    z(2) = CMPLX(-a(2)/a(3), ZERO,DP)    ! remainder is a linear eq.
    outputCode=21   ! two identical real roots
    RETURN
  END IF

  d = a(2)*a(2) - FOUR*a(1)*a(3)             ! the discriminant
  IF (ABS(d) <= TWO*eps*a(2)*a(2)) THEN
    z(1) = CMPLX(-HALF*a(2)/a(3), ZERO, DP) ! discriminant is tiny
    z(2) = z(1)
    outputCode=22  ! two distinct real roots
    RETURN
  END IF

  r = SQRT(ABS(d))
  IF (d < ZERO) THEN
    x = -HALF*a(2)/a(3)        ! negative discriminant => roots are complex   
    y = ABS(HALF*r/a(3))
    z(1) = CMPLX(x, y, DP)
    z(2) = CMPLX(x,-y, DP)   ! its conjugate
    outputCode=23                        !  COMPLEX ROOTS
    RETURN
  END IF

  IF (a(2) /= ZERO) THEN              ! see Numerical Recipes, sec. 5.5
    w = -(a(2) + SIGN(r,a(2)))
    z(1) = CMPLX(TWO*a(1)/w,  ZERO, DP)
    z(2) = CMPLX(HALF*w/a(3), ZERO, DP)
    outputCode=22           ! two real roots
    RETURN
  END IF

  x = ABS(HALF*r/a(3))   ! a(2)=0 if you get here
  z(1) = CMPLX( x, ZERO, DP)
  z(2) = CMPLX(-x, ZERO, DP)
  outputCode=22
  RETURN
END Subroutine QuadraticRoots   ! -------------------------------------------

!+
SUBROUTINE CubicRoots(a, z)
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z

  REAL(DP),PARAMETER:: RT3=1.7320508075689D0    ! (Sqrt(3)
  REAL (DP) :: aq(3), arg, c, cf, d, p, p1, q, q1
  REAL(DP):: r, ra, rb, rq, rt
  REAL(DP):: r1, s, sf, sq, sum, t, tol, t1, w
  REAL(DP):: w1, w2, x, x1, x2, x3, y, y1, y2, y3

! NOTE -   It is assumed that a(4) is non-zero. No test is made here.
!----------------------------------------------------------------------------
  IF (a(1)==0.0) THEN
    z(1) = CZERO  ! one root is obviously zero
    CALL QuadraticRoots(a(2:4), z(2:3))   ! remaining 2 roots here
    RETURN
  END IF

  p = a(3)/(THREE*a(4))
  q = a(2)/a(4)
  r = a(1)/a(4)
  tol = FOUR*EPS

  c = ZERO
  t = a(2) - p*a(3)
  IF (ABS(t) > tol*ABS(a(2))) c = t/a(4)

  t = TWO*p*p - q
  IF (ABS(t) <= tol*ABS(q)) t = ZERO
  d = r + p*t
  IF (ABS(d) <= tol*ABS(r)) GO TO 110

!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

  s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
  p1 = a(3)/(THREE*s)
  q1 = a(2)/s
  r1 = a(1)/s

  t1 = q - 2.25D0*p*p
  IF (ABS(t1) <= tol*ABS(q)) t1 = ZERO
  w = FOURTH*r1*r1
  w1 = HALF*p1*r1*t
  w2 = q1*q1*t1/27.0D0

  IF (w1 >= ZERO) THEN
    w = w + w1
    sq = w + w2
  ELSE IF (w2 < ZERO) THEN
    sq = w + (w1 + w2)
  ELSE
    w = w + w2
    sq = w + w1
  END IF

  IF (ABS(sq) <= tol*w) sq = ZERO
  rq = ABS(s/a(4))*SQRT(ABS(sq))
  IF (sq >= ZERO) GO TO 40

!                   ALL ROOTS ARE REAL

  arg = ATAN2(rq, -HALF*d)
  cf = COS(arg/THREE)
  sf = SIN(arg/THREE)
  rt = SQRT(-c/THREE)
  y1 = TWO*rt*cf
  y2 = -rt*(cf + rt3*sf)
  y3 = -(d/y1)/y2

  x1 = y1 - p
  x2 = y2 - p
  x3 = y3 - p

  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)
  IF (ABS(x2) > ABS(x3)) CALL Swap(x2,x3)
  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)

  w = x3

  IF (ABS(x2) < 0.1D0*ABS(x3)) GO TO 70
  IF (ABS(x1) < 0.1D0*ABS(x2)) x1 = - (r/x3)/x2
  z(1) = CMPLX(x1, ZERO,DP)
  z(2) = CMPLX(x2, ZERO,DP)
  z(3) = CMPLX(x3, ZERO,DP)
  RETURN

!                  REAL AND COMPLEX ROOTS

40 ra =CubeRoot(-HALF*d - SIGN(rq,d))
  rb = -c/(THREE*ra)
  t = ra + rb
  w = -p
  x = -p
  IF (ABS(t) <= tol*ABS(ra)) GO TO 41
  w = t - p
  x = -HALF*t - p
  IF (ABS(x) <= tol*ABS(p)) x = ZERO
  41 t = ABS(ra - rb)
  y = HALF*rt3*t
  
  IF (t <= tol*ABS(ra)) GO TO 60
  IF (ABS(x) < ABS(y)) GO TO 50
  s = ABS(x)
  t = y/x
  GO TO 51
50 s = ABS(y)
  t = x/y
51 IF (s < 0.1D0*ABS(w)) GO TO 70
  w1 = w/s
  sum = ONE + t*t
  IF (w1*w1 < 0.01D0*sum) w = - ((r/sum)/s)/s
  z(1) = CMPLX(w, ZERO,DP)
  z(2) = CMPLX(x, y,DP)
  z(3) = CMPLX(x,-y,DP)
  RETURN

!               AT LEAST TWO ROOTS ARE EQUAL

60 IF (ABS(x) < ABS(w)) GO TO 61
  IF (ABS(w) < 0.1D0*ABS(x)) w = - (r/x)/x
  z(1) = CMPLX(w, ZERO,DP)
  z(2) = CMPLX(x, ZERO,DP)
  z(3) = z(2)
  RETURN
  61 IF (ABS(x) < 0.1D0*ABS(w)) GO TO 70
  z(1) = CMPLX(x, ZERO,DP)
  z(2) = z(1)
  z(3) = CMPLX(w, ZERO,DP)
  RETURN

!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

70 aq(1) = a(1)
  aq(2) = a(2) + a(1)/w
  aq(3) = -a(4)*w
  CALL QuadraticRoots(aq, z)
  z(3) = CMPLX(w, ZERO,DP)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3) = z(2)
  z(2) = z(1)
  z(1) = CMPLX(w, ZERO,DP)
  RETURN
!-----------------------------------------------------------------------


!                   CASE WHEN D = 0

110 z(1) = CMPLX(-p, ZERO,DP)
  w = SQRT(ABS(c))
  IF (c < ZERO) GO TO 120
  z(2) = CMPLX(-p, w,DP)
  z(3) = CMPLX(-p,-w,DP)
  RETURN

120 IF (p /= ZERO) GO TO 130
  z(2) = CMPLX(w, ZERO,DP)
  z(3) = CMPLX(-w, ZERO,DP)
  RETURN

130 x = -(p + SIGN(w,p))
  z(3) = CMPLX(x, ZERO,DP)
  t = THREE*a(1)/(a(3)*x)
  IF (ABS(p) > ABS(t)) GO TO 131
  z(2) = CMPLX(t, ZERO,DP)
  RETURN
131 z(2) = z(1)
  z(1) = CMPLX(t, ZERO,DP)
  RETURN
END Subroutine CubicRoots   ! -----------------------------------------------


!+
SUBROUTINE QuarticRoots(a,z)
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!               A(1) + A(2)*Z + ... + A(5)*Z**4

  REAL(DP), INTENT(IN)     :: a(:)
  COMPLEX(DP), INTENT(OUT) :: z(:)

  COMPLEX(DP) :: w
  REAL(DP):: b,b2, c, d, e, h, p, q, r, t
  REAL(DP),DIMENSION(4):: temp
  REAL(DP):: u, v, v1, v2, x, x1, x2, x3, y


! NOTE - It is assumed that a(5) is non-zero. No test is made here

!----------------------------------------------------------------------------

  IF (a(1)==0.0) THEN
    z(1) = CZERO    !  one root is obviously zero
    CALL CubicRoots(a(2:), z(2:))
    RETURN
  END IF


  b = a(4)/(FOUR*a(5))
  c = a(3)/a(5)
  d = a(2)/a(5)
  e = a(1)/a(5)
  b2 = b*b

  p = HALF*(c - 6.0D0*b2)
  q = d - TWO*b*(c - FOUR*b2)
  r = b2*(c - THREE*b2) - b*d + e

! SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
! NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
! THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE
!     Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
! WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
! THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.

  temp(1) = -q*q/64.0D0
  temp(2) = 0.25D0*(p*p - r)
  temp(3) =  p
  temp(4) = ONE
  CALL CubicRoots(temp,z)
  IF (AIMAG(z(2)) /= ZERO) GO TO 60

!         THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
!         REORDER THE ROOTS IN INCREASING ORDER

  x1 = DBLE(z(1))
  x2 = DBLE(z(2))
  x3 = DBLE(z(3))
  IF (x1 > x2) CALL Swap(x1,x2)
  IF (x2 > x3) CALL Swap(x2,x3)
  IF (x1 > x2) CALL Swap(x1,x2)

  u = ZERO
  IF (x3 > ZERO) u = SQRT(x3)
  IF (x2 <= ZERO) GO TO 41
  IF (x1 >= ZERO) GO TO 30
  IF (ABS(x1) > x2) GO TO 40
  x1 = ZERO

30 x1 = SQRT(x1)
  x2 = SQRT(x2)
  IF (q > ZERO) x1 = -x1
  temp(1) = (( x1 + x2) + u) - b
  temp(2) = ((-x1 - x2) + u) - b
  temp(3) = (( x1 - x2) - u) - b
  temp(4) = ((-x1 + x2) - u) - b
  CALL SelectSort(temp)
  IF (ABS(temp(1)) >= 0.1D0*ABS(temp(4))) GO TO 31
  t = temp(2)*temp(3)*temp(4)
  IF (t /= ZERO) temp(1) = e/t
31 z(1) = CMPLX(temp(1), ZERO,DP)
  z(2) = CMPLX(temp(2), ZERO,DP)
  z(3) = CMPLX(temp(3), ZERO,DP)
  z(4) = CMPLX(temp(4), ZERO,DP)
  RETURN

40 v1 = SQRT(ABS(x1))
v2 = ZERO
GO TO 50
41 v1 = SQRT(ABS(x1))
v2 = SQRT(ABS(x2))
IF (q < ZERO) u = -u

50 x = -u - b
y = v1 - v2
z(1) = CMPLX(x, y,DP)
z(2) = CMPLX(x,-y,DP)
x =  u - b
y = v1 + v2
z(3) = CMPLX(x, y,DP)
z(4) = CMPLX(x,-y,DP)
RETURN

!                THE RESOLVENT CUBIC HAS COMPLEX ROOTS

60 t = DBLE(z(1))
x = ZERO
IF (t < ZERO) THEN
  GO TO 61
ELSE IF (t == ZERO) THEN
  GO TO 70
ELSE
  GO TO 62
END IF
61 h = ABS(DBLE(z(2))) + ABS(AIMAG(z(2)))
IF (ABS(t) <= h) GO TO 70
GO TO 80
62 x = SQRT(t)
IF (q > ZERO) x = -x

70 w = SQRT(z(2))
  u = TWO*DBLE(w)
  v = TWO*ABS(AIMAG(w))
  t =  x - b
  x1 = t + u
  x2 = t - u
  IF (ABS(x1) <= ABS(x2)) GO TO 71
  t = x1
  x1 = x2
  x2 = t
71 u = -x - b
  h = u*u + v*v
  IF (x1*x1 < 0.01D0*MIN(x2*x2,h)) x1 = e/(x2*h)
  z(1) = CMPLX(x1, ZERO,DP)
  z(2) = CMPLX(x2, ZERO,DP)
  z(3) = CMPLX(u, v,DP)
  z(4) = CMPLX(u,-v,DP)
  RETURN

80 v = SQRT(ABS(t))
  z(1) = CMPLX(-b, v,DP)
  z(2) = CMPLX(-b,-v,DP)
  z(3) = z(1)
  z(4) = z(2)
  RETURN

END Subroutine QuarticRoots

!+
SUBROUTINE SelectSort(a)
! ---------------------------------------------------------------------------
! PURPOSE - Reorder the elements of in increasing order.
  REAL(DP),INTENT(IN OUT),DIMENSION(:):: a

  INTEGER:: j
  INTEGER,DIMENSION(1):: k
! NOTE - This is a n**2 method. It should only be used for small arrays. <25
!----------------------------------------------------------------------------
  DO j=1,SIZE(a)-1
    k=MINLOC(a(j:))
    IF (j /= k(1)) CALL Swap(a(k(1)),a(j))
  END DO
  RETURN
END Subroutine SelectSort   ! -----------------------------------------------

!+
SUBROUTINE SolvePolynomial(quarticCoeff, cubicCoeff, quadraticCoeff, &
  linearCoeff, constantCoeff, code, root1,root2,root3,root4)
! ---------------------------------------------------------------------------
  REAL(DP),INTENT(IN):: quarticCoeff
  REAL(DP),INTENT(IN):: cubicCoeff, quadraticCoeff
  REAL(DP),INTENT(IN):: linearCoeff, constantCoeff
  INTEGER,INTENT(OUT):: code
  COMPLEX(DP),INTENT(OUT):: root1,root2,root3,root4
  REAL(DP),DIMENSION(5):: a
  COMPLEX(DP),DIMENSION(5):: z
!----------------------------------------------------------------------------
  a(1)=constantCoeff
  a(2)=linearCoeff
  a(3)=quadraticCoeff
  a(4)=cubicCoeff
  a(5)=quarticCoeff

  IF (quarticCoeff /= ZERO) THEN
    CALL QuarticRoots(a,z)  
  ELSE IF (cubicCoeff /= ZERO) THEN
    CALL CubicRoots(a,z)
  ELSE IF (quadraticCoeff /= ZERO) THEN
    CALL QuadraticRoots(a,z)
  ELSE IF (linearCoeff /= ZERO) THEN
    z(1)=CMPLX(-constantCoeff/linearCoeff, 0, DP)
    outputCode=1
  ELSE
    outputCode=0    !  { no roots }
  END IF

  code=outputCode
  IF (outputCode > 0) root1=z(1)
  IF (outputCode > 1) root2=z(2)
  IF (outputCode > 23) root3=z(3)
  IF (outputCode > 99) root4=z(4)
  RETURN
END Subroutine SolvePolynomial   ! ------------------------------------------

!+
SUBROUTINE SwapDouble(a,b)
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(DP),INTENT(IN OUT):: a,b
  REAL(DP):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine SwapDouble   ! -----------------------------------------------

!
SUBROUTINE SwapSingle(a,b)
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(SP),INTENT(IN OUT):: a,b
  REAL(SP):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine SwapSingle   ! -----------------------------------------------


END Module PolynomialRoots   ! ==





! C                                                                                                !
! C     CALCULO NUMERICO DA REFERENCIA PHYSICA B 404 (2009) 3102                                   !
! C------------------------------------------------------------------------------------------------!
! C     PROGRAMA PARA CALCULAR O PARAMETRO DE ORDEM SUPERCONDUTOR                                  !
! C     USANDO APROXIMACAO HUBBARD I                                                               !
! C                                                                                                !
! C                                                                                                !
! C                                                                                                !                                                                             !
! C------------------------------------------------------------------------------------------------!









program fullmodel_v2
  use cubic
  IMPLICIT NONE
  real (KIND = 8) :: kx,ky,kz,aw1,aw2,aw3,boltz,kbt,pq_write,enek 
    real (KIND = 8) :: limx,limy,qreal,e1k,e2k,e3k,ireal,zzreal,nword,eitareal
    integer (KIND = 8) :: qq,ii=npt,zz,iti
    V = VV
    VTI = valn*V*2.0d0
    ezero = -VTI/5.65d0
    un = U*valn
    ireal = real(ii,8)
    OPEN(45,file='w1nt85_U02V8.dat')
    OPEN(46,file='w2nt85_U02V8.dat')
    OPEN(47,file='w3nt85_U02V8.dat')
    OPEN(48,file='w4nt85_U02V8.dat')
    OPEN(49,file='w5nt85_U02V8.dat')
    OPEN(50,file='w6nt85_U02V8.dat')
    print*, U,V
    read*,
    !OPEN(45,file='qw1_V0_U4_del0_nt0.85.dat')
    !OPEN(46,file='qw2_V0_U4_del0_nt0.85.dat')
    !OPEN(47,file='qw3_V0_U4_del0_nt0.85.dat')
    !OPEN(48,file='qw4_V0_U4_del0_nt0.85.dat')
    !OPEN(49,file='qw5_V0_U4_del0_nt0.85.dat')
    !OPEN(50,file='qw6_V0_U4_del0_nt0.85.dat')    
    !delta0 = 0.0d0
      iti = 1
    SELECT CASE(iti)

    CASE(1)
    79        format(1x,F7.2,1x,F9.3,1X,I5.2,1X)
        
        limx = pi/ireal
        limy = pi/ireal
        kx = 0.0d0
        ky = 0.0d0
        kz = 0.0d0
        do qq=0,ii
            qreal=real(qq,8)
            kx = limx*qreal
            ky = limy*qreal
            call Bands(kx,ky,e1k,e2k,e3k,delta0,valin)
            print*, e1k,e2k,e3k
            write(45,79) dsqrt((kx**2.0d0) + (ky**2.0d0)),e1k,qq
            write(46,79) dsqrt((kx**2.0d0) + (ky**2.0d0)),e2k,qq
            write(47,79) dsqrt((kx**2.0d0) + (ky**2.0d0)),e3k,qq
            write(48,79) dsqrt((kx**2.0d0) + (ky**2.0d0)),-e1k,qq
            write(49,79) dsqrt((kx**2.0d0) + (ky**2.0d0)),-e2k,qq
            write(50,79) dsqrt((kx**2.0d0) + (ky**2.0d0)),-e3k,qq
            
        end do
        
        
        do qq=0,ii
            qreal=real(qq,8)
            ky = pi-(limy*qreal)
            print*, e1k,e2k,e3k
            call Bands(kx,ky,e1k,e2k,e3k,delta0,valin)
            write(45,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(limy*qreal),e1k,ii+qq
            write(46,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(limy*qreal),e2k,qq+ii
            write(47,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(limy*qreal),e3k,qq+ii
            write(48,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(limy*qreal),-e1k,qq+ii
            write(49,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(limy*qreal),-e2k,qq+ii
            write(50,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(limy*qreal),-e3k,qq+ii
            
        end do 
        do qq=0,ii
            qreal=real(qq,8)
            kx = pi-(limx*qreal)     
            print*, e1k,e2k,e3k
            call Bands(kx,ky,e1k,e2k,e3k,delta0,valin)
            write(45,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(pi)+(limx*qreal),e1k,ii+qq
            write(46,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(pi)+(limx*qreal),e2k,qq+ii
            write(47,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(pi)+(limx*qreal),e3k,qq+ii
            write(48,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(pi)+(limx*qreal),-e1k,qq+ii
            write(49,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(pi)+(limx*qreal),-e2k,qq+ii
            write(50,79) dsqrt((pi**2.0d0) + (pi**2.0d0))+(pi)+(limx*qreal),-e3k,qq+ii
      end do

  CASE(2) 
    print*, "par_ordem"
    call par_ordem(pontos,valin)
  
  CASE(3)
    78     format(1x,F20.15,1x,F20.15,1x,F20.15) 
    OPEN(60,file='pq_n0.85_V8U01kbt2var.dat')
    print*, "pquio"
     do eita = 1, pontos 
      ireal = real(pontos,8)
      eitareal = real(eita,8)
      V = 10.0d0*eitareal/ireal 
     do qq = 1, pontos
        
        qreal = real(qq,8)
        enek = 1.00d0*qreal/ireal 
        kbt = 0.50d0 
        !kbt = 10.0d0*qreal/ireal
        boltz = 1.0d0/kbt
      
        ! do qq = 1, nnpt
    !    qqreal = real(qq,8)
    !    betaizi = 10.0d0*qqreal/ppreal
    !    betaizi = 1.0d0/betaizi 
    !call pquimico(pontos,delta0,boltz,pq_write,enek,nword)
    call bissec_quim(-10000.0d0,2000000.0d0,pq_write,delta0,boltz,enek)


    write(60,78) V, enek, pq_write 
  end do 
end do 
  END SELECT
end program fullmodel_v2



SUBROUTINE Bands(akx,aky,w1,w2,w3,deltaz,nk)
  use cubic 
  IMPLICIT NONE 
  real (KIND = 8) :: ap,aq,ar,a1,a2,d1,ac0,ac1,pwrd,uti,nk
  real (KIND = 8) :: w1,w2,w3,akx,aky,akz,az,d2,dx,l2_a,l2_b,w4,l4 
  real (KIND = 8) :: l1,l2,l3,l0,ab1,ab2,ab3,deta=1d-12,ax12,ax13,ax23,ax32
  real (KIND = 8) :: vt,l1_a,l1_b,l1_c,g3_up,g4_up,g6_up,deltaz,w11,w22,w33
    VTI = 2.0d0*V*nk
    ezero = -VTI/5.65
    ek = -2.0d0*t1*(dcos(akx) + dcos(aky)) + 4.0d0*t2*dcos(akx)*dcos(aky) + ezero
    gammak = dabs(dcos(akx) - dcos(aky))

    deltak = 2.0d0*deltaz*gammak
    valn = nk
    un = U*nk
    
    !print*, V,U,valn

    !dx = 4.0d0*deltak*deltak
 !   ab1 = VTI + U
 !   ab2 = VTI + U*valn 
 !   ab3 = VTI + U*(1.0d0-valn)
    !l0 = -((ek*ek*(un**2.0d0)*(ab3**2.0d0)) + 4.0d0*deltak*deltak*(-2.0d0*un*VTI*(VTI*VTI + VTI*un) - VTI**4.0d0 - un*un*ab2*ab2))  
    !l1 = ab1*ab1*ab2*ab2 +  ek*ek*(un*un + ab3*ab3) + 2.0d0*ek*(ab1*ab3*VTI + un*un*ab2) - 4.0d0*deltak*deltak*ab2*ab2
    !l2 = -(ab1*ab1 + ab2*ab2 + ek*ek + 2.0d0*ek*ab2)
    l3 = 1.0d0
    l4 = 0.0d0


 l0= -4.0d0*deltak**2.0d0*VTI**4.0d0 - 8.0d0*deltak**2.0d0*VTI**3.0d0*U - 4.0d0*deltak**2.0d0*VTI**2.0d0*U**2.0d0& 
     -  8.0d0*valn**2.0d0*deltak**2.0d0*VTI**2.0d0*U**2.0d0 - valn**2.0d0*ek**2.0d0*VTI**2.0d0*U**2.0d0& 
     -  8.0d0*valn**2.0d0*deltak**2.0d0*VTI*U**3.0d0 - 2.0d0*valn**2.0d0*ek**2.0d0*VTI*U**3.0d0&
     +  2.0d0*valn**3.0d0*ek**2.0d0*VTI*U**3.0d0 - 4.0d0*valn**4.0d0*deltak**2.0d0*U**4.0d0 - valn**2.0d0*ek**2.0d0*U**4.0d0 & 
     +  2.0d0*valn**3.0d0*ek**2.0d0*U**4.0d0 - valn**4.0d0*ek**2.0d0*U**4.0d0
     
        
     l1=4.0d0*deltak**2.0d0*VTI**2.0d0 + ek**2.0d0*VTI**2.0d0 + 2.0d0*ek*VTI**3.0d0 + VTI**4.0d0 & 
     +  8.0d0*valn*deltak**2.0d0*VTI*U + 2.0d0*ek**2.0d0*VTI*U - 2.0d0*valn*ek**2.0d0*VTI*U + 4.0d0*ek*VTI**2.0d0*U & 
     -  2.0d0*valn*ek*VTI**2.0d0*U + 2.0d0*VTI**3.0d0*U + 2.0d0*valn*VTI**3.0d0*U + 4.0d0*valn**2.0d0*deltak**2.0d0*U**2.0d0 &
     +  ek**2.0d0*U**2.0d0 - 2.0d0*valn*ek**2.0d0*U**2.0d0 + 2.0d0*valn**2.0d0*ek**2.0d0*U**2.0d0 + 2.0d0*ek*VTI*U**2.0d0 & 
     -  2.0d0*valn*ek*VTI*U**2.0d0 + 2.0d0*valn**2.0d0*ek*VTI*U**2.0d0 + VTI**2.0d0*U**2.0d0 + 4.0d0*valn*VTI**2.0d0*U**2.0d0 & 
     +  valn**2.0d0*VTI**2.0d0*U**2.0d0 + 2.0d0*valn**3.0d0*ek*U**3.0d0 + 2.0d0*valn*VTI*U**3.0d0 + 2.0d0*valn**2.0d0*VTI*U**3.0d0 & 
     +  valn**2.0d0*U**4.0d0
    
      
     l2 = -ek**2.0d0 - 2.0d0*ek*VTI - 2.0d0*VTI**2.0d0 - 2.0d0*valn*ek*U - 2.0d0*VTI*U - 2.0d0*valn*VTI*U - U**2.0d0 
     l2 = l2 - (valn**2.0d0)*(U**2.0d0)

    !print*, l0,l1,l2,l3,l4
    !read*,  
    call QuarticSolver(l0,l1,l2,l3,l4,w11,w22,w33,w4)

    if(dabs(w11) .lt. 1d-5) then
      w11 = 0.0D0
    end if 
    if(dabs(w22) .lt. 1d-5) then
      w22 = 0.0D0
    end if
    if(dabs(w33) .lt. 1d-5) then
      w33 = 0.0D0
    end if  



      w1 = dsqrt(w11)
      w2 = dsqrt(w22)
      w3 = dsqrt(w33)   
     
  !    print*, w1
  !    print*, w2
  !    print*, w3
  !    read*,
   
         if (isnan(w1))then
           print*,'w1,w11 =',w1,w11 
           read(*,*)
        endif             
       if (isnan(w3))then
           print*,'w3,w22 =',w2,w22 
           read(*,*)
        endif             

       if (isnan(w3))then
           print*,'w3,w33 =',w3,w33 
           read(*,*)
        endif             

! if (dabs(w1) .lt. 1.0d-8) then
!        w1 = 0.0d0
!    end if
!    if (dabs(w2) .lt. 1.0d-8) then
!        w2 = 0.0d0
!    end if
!    if (dabs(w3) .lt. 1.0d-8) then
!        w3 = 0.0d0
!    end if
    !read*,
    !ax12 = (w1*w1) - (w2*w2)
    !ax13 = (w1*w1) - (w3*w3)
    !ax23 = (w2*w2) - (w3*w3)
    !ax32 = (w3*w3) - (w2*w2)
    !print*, w1,w2,w3, 'w'
    !print*, ax12, ax23,ax13, 'ax'
    !read*,

    !     (-l0 + l1*w1 - l2*w1**2 + l3*w1**3 - l4*w1**4)/
    !     (2.0d0*w1*(w1 - w2)*(w1 + w2)*(w1 - w3)*(w1 + w3))!*(w + w1)  
    !  -  (l0 + l1*w1 + l2*w1**2 + l3*w1**3 + l4*w1**4)/
    !     (2.0d0*w1*(w1 - w2)*(w1 + w2)*(w1 - w3)*(w1 + w3))!*(w - w1) 
    !  -  (l0 - l1*w2 + l2*w2**2 - l3*w2**3 + l4*w2**4)/
    !     (2.0d0*(w1 - w2)*w2*(w1 + w2)*(w2 - w3)*(w2 + w3))!*(w + w2) + 
    !  -  (l0 + l1*w2 + l2*w2**2 + l3*w2**3 + l4*w2**4)/
    !     (2.0d0*w2*(-w1 + w2)*(w1 + w2)*(w2 - w3)*(w2 + w3))!*(w - w2) + 
    !  -  (l0 - l1*w3 + l2*w3**2 - l3*w3**3 + l4*w3**4)/
    !     (2.0d0*(w1 - w3)*w3*(w1 + w3)*(-w2 + w3)*(w2 + w3))!*(w + w3) + 
    !  -  (l0 + l1*w3 + l2*w3**2 + l3*w3**3 + l4*w3**4)/
    !     (2.0d0*w3*(-w1 + w3)*(w1 + w3)*(-w2 + w3)*(w2 + w3))!*(w - w3)

    g3 = (2.0d0*w1*(w1 - w2)*(w1 + w2)*(w1 - w3)*(w1 + w3))
    g4 = (2.0d0*(w1 - w2)*w2*(w1 + w2)*(w2 - w3)*(w2 + w3))
    g6 = (2.0d0*(w1 - w3)*w3*(w1 + w3)*(-w2 + w3)*(w2 + w3))
   ! print*, g3,g4,g6 


    !g3_up = 2.0d0*w1*w1*(ax13*ax12)*(ax13*ax12) + deta*deta 
    !g3 = w1*(ax13*ax12)/g3_up 
    !g4_up = 2.0d0*w2*w2*(ax12*ax23)*(ax12*ax23) + deta*deta
    !g4 = w2*(ax12*ax23)/g4_up 
    !g6_up = 2.0d0*w3*w3*(ax13*ax32)*(ax13*ax32) + deta*deta
    !g6 = -w3*(ax13*ax32)/g6_up 
    if (dabs(g3) .lt. deta) then
        g3 = 0.0d0
    end if
    if (dabs(g4) .lt. deta) then
        g6 = 0.0d0
    end if
    if (dabs(g6) .lt. deta) then
        g6 = 0.0d0
    end if



END SUBROUTINE Bands





SUBROUTINE QuarticSolver(a0,a1,a2,a3,a4,z1,z2,z3,z4)
! ---------------------------------------------------------------------------
! PURPOSE - Test the functioningg of the procedures in the module
!    PolynomialRoots and compare to those in the original module
!    by Alan Miller
! AUTHORS - Alan Miller,  CSIRO Mathematical & Information Sciences
!                         CLAYTON, VICTORIA, AUSTRALIA 3169
!                         http://www.mel.dms.csiro.au/~alan
!           Ralph L. Carmichael, Public Domain Aeronautical Software
!                         http://www.pdas.com
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 27Feb97  1.0   AM    Original coding
! 16Jul98  1.1   RLC   Numerous style changes

  USE PolynomialRoots
  IMPLICIT NONE

  INTEGER,PARAMETER:: DP=KIND(1.0D0)
  CHARACTER(LEN=*),PARAMETER:: FMT = "(2ES20.12)"

!  CHARACTER(LEN=80):: buffer
  INTEGER:: i
  REAL(dp),DIMENSION(0:4):: a
  COMPLEX(dp),DIMENSION(4):: z
  REAL(dp),intent(in) :: a0,a1,a2,a3,a4 
  REAL(dp),intent(out) :: z1,z2,z3,z4
!----------------------------------------------------------------------------
!  WRITE(*,*) "Solve quadratic, cubic or quartic eqns. with real coefficients"
!  WRITE(*,*) " "
!  WRITE(*,*) "  A*x**4 + B*x**3 + C*x**2 + D*x + E = 0"
!  WRITE(*,*) " "

!  WRITE(*,*) "Enter the coefficient of the constant term:"
!  READ(*,'(A)') buffer
!  READ(buffer,'(F80.0)') a(0)

!  WRITE(*,*) "Enter the coefficient of the linear term:"
!  READ(*,'(A)') buffer
!  READ(buffer,'(F80.0)') a(1)

!  WRITE(*,*) "Enter the coefficient of the quadratic term:"
!  READ(*,'(A)') buffer
!  READ(buffer,'(F80.0)') a(2)

!  WRITE(*,*) "Enter the coefficient of the cubic term:"
!  READ(*,'(A)') buffer
!  READ(buffer,'(F80.0)') a(3)

!  WRITE(*,*) "Enter the coefficient of the quartic term:"
!  READ(*,'(A)') buffer
!  READ(buffer,'(F80.0)') a(4)
  a(0)=a0
  a(1)=a1
  a(2)=a2
  a(3)=a3
  a(4)=a4
 
 
 


  IF (a(4) /= 0.0) THEN
    CALL QuarticRoots(a, z)
!    WRITE(*,*) " Roots: REAL PART   IMAGINARY PART"
!    WRITE(*,FMT) (DBLE(z(i)), AIMAG(z(i)), i=1,4)
     z1=DBLE(z(1))
     z2=DBLE(z(2))
     z3=DBLE(z(3))
     z4=DBLE(z(4))

  ELSE IF(a(3) /= 0.0) THEN
    CALL CubicRoots(a, z)
!    WRITE(*,*) " Roots: REAL PART   IMAGINARY PART"
!    WRITE(*,FMT) (DBLE(z(i)), AIMAG(z(i)), i=1,3)
     z1=DBLE(z(1))
     z2=DBLE(z(2))
     z3=DBLE(z(3))
     z4=0.0d0


  ELSE IF(a(2) /= 0.0) THEN
    CALL QuadraticRoots(a, z)
!    WRITE(*,*) " Roots: REAL PART   IMAGINARY PART"
!    WRITE(*,FMT) (DBLE(z(i)), AIMAG(z(i)), i=1,2)
     z1=DBLE(z(1))
     z2=DBLE(z(2))
     z3=0.0d0
     z4=0.0d0


  ELSE IF(a(1) /= 0.0) THEN
    z(1)=CMPLX(-a(0)/a(1), 0.0, DP)      
!    WRITE(*,*) " Roots: REAL PART   IMAGINARY PART"
!    WRITE(*,FMT) (DBLE(z(i)), AIMAG(z(i)), i=1,1)
     z1=DBLE(z(1))
     z2=0.0d0
     z3=0.0d0
     z4=0.0d0

  ELSE
    WRITE(*,*) "The equation is degenerate."
  END IF
return
END SUBROUTINE QuarticSolver


subroutine par_ordem(pts,ene)
    use cubic
    IMPLICIT NONE 
    real(KIND = 8) :: vx,vy,limax,limay,izireal,ggreal,hhreal,deltra,delta,temp,ktemp
    real(KIND = 8) :: soma_par,deltaold,err_delta = 1d-8,deldif,stotal,ttreal,q0_a,q0_b,q0,q2,q4,q1
    real(KIND = 8) :: p3=0.0d0,p4=0.0d0,p6=0.0d0,tg1,tg2,tg3,lam0,lam2,lam4,lamtot,fermi
    real(KIND = 8) :: d1,d2,d3,d4,d5,d6,d1a,d1b,d2a,d2b,d3a,d3b,d4a,d4b,d5a,d5b,d6a,d6b
    real(KIND = 8) :: f1,f2,f3,f4,f5,f6,df,pqmu,ene,q3 
    integer (KIND = 8) :: hh,gg,izi,pts,tt!,check
    OPEN(51,file='idelta_n0.85_V2U6.dat')
    78     format(1x,F20.15,1x,F20.15,1x) 
    izi = pts
    ptreal = real(izi*izi,8)
    izireal = real(izi,8)
    limax = pi/izireal
    limay = pi/izireal
    valn = ene 
    VTI = 2.0d0*V*valn
    un = U*valn 
    ezero = -VTI/5.65d0
    soma_par = 0.0d0 
    !deltaold = 0.0001d0
    !q4 = 1.0d0
   ! q2 = (U + VTI)**2.0d0 + 2.0d0*ek*(VTI + U) + un*(2.0d0*VTI + U*(1.0d0 - valn))
    !q2 = -(U + VTI)**2.0d0 + un*un - U*un - VTI*U
    !q0_a = ek*un*(VTI+U)*(VTI+U*(1.0d0 - valn)) + 2.0d0*un*VTI*(VTI + U*(1.0d0 - valn)) 
    !q0_b = -un*(VTI + U*(1.0d0 - valn))*(VTI + un)
    !q0_a = -un*ek*(VTI+u)*(VTI+U*(1-valn)) + (VTI*U*un)*(1-valn)*(VTI+un) + un*un*un*(VTI + U*(1-valn))
    !q0 = q0_a 
   ! print*, q4,q2,q0
   ! read*,


   pqmu = 0.0d0
    do tt = 0,izi
        ttreal = real(tt,8)
        ktemp = (10.0d0/izireal)*(1.0d0+ttreal)
        if (ktemp .eq. 0.0d0) then
        temp = 0.0d0
        else  
        temp = 1.0d0/ktemp
        end if
        check = 0
        deltaold = 0.0001d0
    do while (check .eq. 0)
        do hh = 0, izi
            hhreal = real(hh,8)  
            vy = -pi + 2.0d0*limay*hhreal
            !print *, vy
            do gg = 0,izi
                ggreal = real(gg,8)
                vx = -pi + 2.0d0*limax*ggreal
                call Bands(vx,vy,p3,p4,p6,deltaold,ene)
                ek = -2.0d0*t1*(dcos(vx) + dcos(vy)) + 4.0d0*t2*dcos(vx)*dcos(vy) + ezero
                q4 = 1.0d0
                q2 = -(U + VTI)**2.0d0 + un*un - U*un - VTI*U
                q0_a = -un*ek*(VTI+u)*(VTI+U*(1.0d0-valn)) + (VTI*U*un)*(1.0d0-valn)*(VTI+un) + un*un*un*(VTI + U*(1.0d0-valn))
                q0 = q0_a 
                q3 = 0.0d0 
                lam0 = q0*(tg1(temp,p3) + tg2(temp,p4) + tg3(temp,p6))
                lam2 = q2*(tg1(temp,p3)*p3*p3 + tg2(temp,p4)*p4*p4 + tg3(temp,p6)*p6*p6)
                lam4 = q4*(tg1(temp,p3)*(p3**4.0d0) + tg2(temp,p4)*(p4**4.0d0) + tg3(temp,p6)*(p6**4.0d0))
                lamtot = lam0 + lam2 + lam4
                
                !d1a = (q0 + p3*q1 + (p3**2.0d0)*q2 + (p3**4.0d0)*q4)
                !d1b = 2.0d0*p3*(p3 - p4)*(p3 + p4)*(p3 - p6)*(p3 + p6)  
                !d1 = (d1a*d1b)/(d1b*d1b + err_delta*err_delta)
                !3-

               !  d2a = (-q0 + p3*q1 - (p3**2.0d0)*q2 - (p3**4.0d0)*q4)
               !  d2b = (2.0d0*p3*(p3 - p4)*(p3 + p4)*(p3 - p6)*(p3 + p6))  
               !  d2 = (d2a*d2b)/(d2b*d2b + err_delta*err_delta)
               !  !3+

               !  d3a = (q0 + p4*q1 + (p4**2.0d0)*q2 + (p4**4.0d0)*q4)
               !  d3b = (2.0d0*p4*(-p3 + p4)*(p3 + p4)*(p4 - p6)*(p4 + p6))  
               !  d3 = (d3a*d3b)*(d3b*d3b + err_delta*err_delta)
               !  !4-

               !  d4a = (q0 - p4*q1 + (p4**2.0d0)*q2 + (p4**4.0d0)*q4)
               !  d4b = (2.0d0*(p3 - p4)*p4*(p3 + p4)*(p4 - p6)*(p4 + p6))  
               !  d4 = (d4a*d4b)/(d4b*d4b + err_delta*err_delta)
               !  !4+

               !  d5a = (q0 + p6*q1 + (p6**2.0d0)*q2 + (p6**4.0d0)*q4)
               !  d5b = (2.0d0*p6*(-p3 + p6)*(p3 + p6)*(-p4 + p6)*(p4 + p6)) 
               !  d5 = (d5a*d5b)/(d5b*d5b + err_delta*err_delta)
               !  !6-

               !  d6a = (q0 - p6*q1 + p6**2*q2 + p6**4*q4)
               !  d6b = (2.0d0*(p3 - p6)*p6*(p3 + p6)*(-p4 + p6)*(p4 + p6))
               !  d6 = (d6a*d6b)/(d6a*d6b + err_delta*err_delta)
               !  !6+  

               !  f1 = d1*fermi(temp,pqmu,p3)
               !  f2 = d2*fermi(temp,pqmu,-p3)
               !  f3 = d3*fermi(temp,pqmu,p4)
               !  f4 = d4*fermi(temp,pqmu,-p4)
               !  f5 = d5*fermi(temp,pqmu,p6)
               !  f6 = d6*fermi(temp,pqmu,-p6)
               !  df = f1+f2+f3+f4+f5+f6
               ! !print*, lamtot, 'lam'
          !     read*,
                deltra = delta(vx,vy,temp,deltaold,lamtot)
                !print*, deltra
                soma_par = deltra + soma_par
                !soma_par = 0.5d0*soma_par
            end do 
        end do

        soma_par = soma_par/ptreal
        soma_par = 0.5d0*soma_par + 0.5d0*deltaold
        deldif = dabs(soma_par - deltaold)
        if (deldif .lt. err_delta) then
            print*, deldif, 'deldif'
            print*, ktemp, 'ktemp'
           ! print*, soma_par, 'sp'
            print*, deltaold,'old'
            print*, soma_par, 'soma_par'
            !print *, stotal, 'stotal'
            write(51,78) ktemp,soma_par
            call flush(51)
            check = 1
            !read*,
        else 
            deltaold = soma_par
        end if 
        soma_par = 0.0d0
    end do

end do 

END SUBROUTINE par_ordem



subroutine pquimico(nnpt,deltaizi,betaizi,pq_sub,ocup,ncorrel)
  use cubic
  IMPLICIT NONE 
  real (KIND = 8) :: pq,pq1,pq2,pq3,pq4,pq5,pq6,bw1,bw2,bw3,bw4,bw5,bw6,corri=1d-10
  real (KIND = 8) :: k0,k1,k2,k3,k4,k5,k6,m1,m2,m3,m4,m5,m6,fermi,dhelta,pdel,k2_a,k2_b,pq_sub
  real (KIND = 8) :: pq1_a,pq1_b,pq2_a,pq2_b,pq3_a,pq3_b,pq4_a,pq4_b,pq5_a,pq5_b,pq6_a,pq6_b
  real (KIND = 8) :: psix, psiy,psx,psy,qqreal,jjreal,llreal,deltaizi,betaizi,ppreal,pqchute,pqsoma=0.0d0
  real (KIND = 8) :: uni,VVTI, ocup,ncorrel
  integer (KIND = 8) :: jj,ll,nnpt,qq,confir
      
    !  OPEN(60,file='pq_n0.85_V80U8.dat')
    !78     format(1x,F20.15,1x,F20.15,1x) 
      ppreal = real(nnpt,8)
      ! gammak = dabs(dcos(psx) - dcos(psy))
      ! pdel = 4.0d0*deltaizi*deltaizi*gammak*gammak
      ! ek = -2.0d0*t1*(dcos(psx) + dcos(psy)) + 4.0d0*t2*dcos(psx)*dcos(psy) + ezero
      ! k0 = -pdel*VTI*(VTI+U)
      ! k1 = (pdel*VTI*(VTI + U)) + (ek*un*un*(VTI + U*(1-valn))**2.0d0)
      ! k2_a = pdel*(VTI+un) + (2.0d0*VTI + U*(1.0d0+valn))*(VTI+U*(1.0d0-valn)) 
      ! k2_b = - un*((VTI+U)**2.0d0)*(VTI+un) - ek*(un*un + (VTI+U*(1.0d0-valn))**2.0d0)
      ! k2 = k2_a + k2_b
      ! k3 = -(VTI+un)*(2.0d0*VTI + U(1.0d0+valn)) + un*(VTI + U*(1.0d0-valn) + (VTI+U)*(VTI + un))
      !psy = -pi + 2.0d0*(psiy*jjreal)/ppreal 
      !psx = -pi + 2.0d0*(psix*llreal)/ppreal
     ! pqchute = 0.0d0
      !pq = pqchute
      ncorrel = 0.0d0
      valn = ocup
      un = U*valn  
      VTI = 2.0d0*V*valn
        
    ! do qq = 1, nnpt
    !    qqreal = real(qq,8)
    !    betaizi = 10.0d0*qqreal/ppreal
    !    betaizi = 1.0d0/betaizi 
        confir = 0
        pqchute = pq_sub
       ! do while (confir .ne. 1)
        !pqsoma = 0.0d0
        do jj = 0,nnpt
          jjreal = real(jj,8)
          psx = -pi + 2.0d0*(pi*jjreal)/ppreal
        !  print*, psx
        !  read*,
            do ll = 0,nnpt
              llreal = real(ll,8)
              psy = -pi + 2.0d0*(pi*llreal)/ppreal
              gammak = dabs(dcos(psx) - dcos(psy))
              pdel = 16.0d0*deltaizi*deltaizi*gammak*gammak
              ek = -2.0d0*t1*(dcos(psx) + dcos(psy)) + 4.0d0*t2*dcos(psx)*dcos(psy) + ezero
              k0 = -pdel*VTI*(VTI+U)
              k1 = (pdel*VTI*(VTI + U)) + (ek*un*un*(VTI + U*(1.0d0-valn))**2.0d0)
              k2_a = pdel*(VTI+un) + (2.0d0*VTI + U*(1.0d0+valn))*(VTI+U*(1.0d0-valn)) 
              k2_b = - un*((VTI+U)**2.0d0)*(VTI+un) - ek*(un*un + (VTI+U*(1.0d0-valn))**2.0d0)
              k2 = k2_a + k2_b
              k3 = -(VTI+un)*(2.0d0*VTI + U*(1.0d0+valn)) + un*(VTI + U*(1.0d0-valn) + (VTI+U)*(VTI + un))
              k4 = VTI + ek + un
              k5 = 1.0d0

              call Bands(psx,psy,bw2,bw4,bw6,deltaizi,ocup)
            !  print*, bw2,bw4,bw6
            !  print*, k0,k1,k2,k3,k4,k5
            !  read*,
              pq1_a = (-k0 + k1*bw2 - k2*bw2**2.0d0 + k3*bw2**3.0d0 - k4*bw2**4.0d0 + k5*bw2**5.0d0)
              pq1_b = (2.0d0*bw2*(bw2 - bw4)*(bw2 + bw4)*(bw2 - bw6)*(bw2 + bw6))  
              pq1 = (pq1_a*pq1_b)/(pq1_b*pq1_b + corri)

              pq2_a = (k0 + k1*bw2 + k2*bw2**2.0d0 + k3*bw2**3.0d0 + k4*bw2**4.0d0 + k5*bw2**5.0d0)
              pq2_b = (2.0d0*bw2*(bw2 - bw4)*(bw2 + bw4)*(bw2 - bw6)*(bw2 + bw6))  
              pq2 = (pq2_a*pq2_b)/(pq2_b*pq2_b + corri)

              pq3_a = (k0 - k1*bw4 + k2*bw4**2.0d0 - k3*bw4**3.0d0 + k4*bw4**4.0d0 - k5*bw4**5.0d0)
              pq3_b = (2.0d0*(bw2 - bw4)*bw4*(bw2 + bw4)*(bw4 - bw6)*(bw4 + bw6))  
              pq3 = (pq3_a*pq3_b)/(pq3_b*pq3_b + corri)

              pq4_a = (k0 + k1*bw4 + k2*bw4**2.0d0 + k3*bw4**3.0d0 + k4*bw4**4.0d0 + k5*bw4**5.0d0)
              pq4_b = (2.0d0*bw4*(-bw2 + bw4)*(bw2 + bw4)*(bw4 - bw6)*(bw4 + bw6))  
              pq4 = (pq4_a*pq4_b)/(pq4_b*pq4_b + corri)

              pq5_a = (k0 - k1*bw6 + k2*bw6**2.0d0 - k3*bw6**3.0d0 + k4*bw6**4.0d0 - k5*bw6**5.0d0)
              pq5_b = (2.0d0*(bw2 - bw6)*bw6*(bw2 + bw6)*(-bw4 + bw6)*(bw4 + bw6))  
              pq5 = (pq5_a*pq5_b)/(pq5_b*pq5_b + corri)

              pq6_a = (k0 + k1*bw6 + k2*bw6**2.0d0 + k3*bw6**3.0d0 + k4*bw6**4.0d0 + k5*bw6**5.0d0)
              pq6_b = (2.0d0*bw6*(-bw2 + bw6)*(bw2 + bw6)*(-bw4 + bw6)*(bw4 + bw6))
              pq6 = (pq6_a*pq6_b)/(pq6_b*pq6_b + corri)

              m1 = pq1*fermi(betaizi,pqchute,-bw2)
              m2 = pq2*fermi(betaizi,pqchute,bw2)
              m3 = pq3*fermi(betaizi,pqchute,-bw4)
              m4 = pq4*fermi(betaizi,pqchute,bw4)
              m5 = pq5*fermi(betaizi,pqchute,-bw6)
              m6 = pq6*fermi(betaizi,pqchute,bw6)
              
              pq = m1+m2+m3+m4+m5+m6
              !print*, pq
              !read*,
              pqsoma = pqsoma+pq
            !  print*, pqsoma
            !  read*,
            end do 
              
          end do 
          ncorrel = pqsoma/(ppreal*ppreal)
          pqsoma = 0.0d0
          !print*, ncorrel, 'ncorr'

         ! print*, pqsoma
       !   pqsoma = 0.5d0*pqsoma + 0.5d0*pqchute
       !   dhelta = dabs(pqsoma - pqchute)
       !   if (dhelta .lt. 1d-9) then
       !    ! write (60,78) 1.0d0/betaizi,pqsoma
        !    pq_sub = pqsoma
        !    print *, pqsoma,1.0d0/betaizi,ocup 
        !    confir = 1
        !  else 
        !  !print*, pqsoma, 'erro'
        !  pqchute = pqsoma 
        !  end if
        !  pqsoma = 0.0d0
         ! pqchute = 0.0d0
   ! end do 
  !  end do   
end 


subroutine bissec_quim(liminf,limsup,pm_h,del_gap,kabete,nigg)
  use cubic 
  implicit none
  real(KIND=8) :: fa,fb,fp,p,erro_1,liminf,limsup,bolz,epson,bolz_1,dummya,dummyb,qlog
  real(KIND=8) :: kabete, del_gap,nigg,xa,xb,pm_h,qdif,qdifa,qdifb,dma,dmb 
  integer :: dummy = 0,ilog,ipi 
  xa = liminf
  xb = limsup
  !open (unit=55,file ='biss.dat')
!pquimico(nnpt,deltaizi,betaizi,pq_sub,ocup,ncorrel)
  dummy = 0
  erro_1 = 1.0d0
  epson = 1.0d0
 ! print*, nigg
  ipi = 0
  p = (xa + xb)/2.0d0
  call pquimico(pontos,del_gap,kabete,xa,nigg,fa)
  call pquimico(pontos,del_gap,kabete,xb,nigg,fb)
  call pquimico(pontos,del_gap,kabete,p,nigg,fp)
 ! print*, nigg, xa,xb,p, 'coord'
  !print*, fa,fb,fp, 'imag'
  qdifa = fa - nigg
  qdifb = fb - nigg
  bolz = (qdifa/dabs(qdifa))*(qdifb/dabs(qdifb)) 
  !print*, qdifa,qdifb, bolz, 'err'
  !bolz = (fa/abs(fa))*(fb/abs(fb))
  qlog = (log10(xb - xa) - log10(1e-9))/log10(2.0d0)
!  do while (abs(erro_1) .gt. 1d-15 .and. abs(epson) .gt. 1d-15) 
  do while (ipi < qlog)
  !do while(dummy .ne. 1 .or. ipi .lt. qlog )
              !  write(55,'(f12.6,6X,f12.6,6X,f12.6)') xa,xb,fp
  qdifa = fa - nigg
  qdifb = fb - nigg
  bolz = (qdifa/dabs(qdifa))*(qdifb/dabs(qdifb))
 ! print*, qdifa,qdifb, 'qdif'
  !read*,
!bolz = (fa/abs(fa))*(fb/abs(fb))
 !  print*, bolz,fa,fb
   !read*, 
  !  print *,'xa = ',xa, 'xb =', xb, 'f(p) =', fp
    dummya = xa + 1d-10
    if(bolz .eq. 1.0d0) then
      xa = p
      !fa = l(xa)
      call pquimico(pontos,del_gap,kabete,xa,nigg,fa)

      !print*, xa,fa
      !read*,
      qdifa = fa - nigg
  qdifb = fb - nigg
  bolz_1 = (qdifa/dabs(qdifa))*(qdifb/dabs(qdifb))
      !bolz_1 = (fa/abs(fa))*(fb/abs(fb))
      if (bolz .ne. bolz_1) then
        xb = xa
        xa = dummya
        call pquimico(pontos,del_gap,kabete,xa,nigg,fa)
        call pquimico(pontos,del_gap,kabete,xb,nigg,fb)
        !fb = l(xb)
        !fa = l(xa)
      end if
    end if
    


    dummyb = xb - 1d-10
    
    if(bolz .eq. -1.0d0) then     
      xb = p      
      !fb = l(xb)
      call pquimico(pontos,del_gap,kabete,xb,nigg,fb)
      qdifa = fa - nigg
  qdifb = fb - nigg
  bolz_1 = (qdifa/dabs(qdifa))*(qdifb/dabs(qdifb))
      !bolz_1 = (fa/abs(fa))*(fb/abs(fb))
    

      if (bolz .ne. bolz_1) then
        xa = xb
        xb = dummyb
        call pquimico(pontos,del_gap,kabete,xa,nigg,fa)
        call pquimico(pontos,del_gap,kabete,xb,nigg,fb)
        !fa = l(xa)
        !fb = l(xb)
      end if
    end if
 !   print*, xa,fa,nigg,'a'
 !   print*, xb,fb,nigg,'b'
  !  print*, p,fp,'p'
    p = (xa+xb)/2.0d0
    !fp = l(p)
    call pquimico(pontos,del_gap,kabete,p,nigg,fp)
    qdif = dabs(dabs(nigg) - dabs(fp))
    ! print*, xa,xb
    !if(abs(fp) .lt. 1e-9) then!#3 
    if(qdif .lt. 1d-8 ) then 
   !   print*, p,fp,nigg, 'p'
      dummy = 1
      
    end if
    ipi = ipi+1
 !   erro_1 = abs(abs(fb) - abs(fa))
 !   epson = abs(abs(xb) - abs(xa))
    
   ! print*, qdif
    
  end do
  !print*,'A raiz tem valor p = ',p
  print*, p,fp,nigg, 'p'
  pm_h = p
end




real (KIND = 8) function delta(rkx,rky,beta,deltax,lambida)
    use cubic 
    IMPLICIT NONE 
    real (KIND = 8) :: ksi,arg,alfa1,alfa2,diff,aux,delta_par
    real (KIND = 8) :: deltax,arg1,beta,rkx,rky,deltakx,lambida
    gammak = dabs(dcos(rkx) - dcos(rky))
    delta_par =  -lambida*V*deltax*2.0d0*gammak
    delta = delta_par
    !return
end 

real (KIND = 8) function tg1(btz,omega3)
    use cubic
    IMPLICIT NONE 
    real (KIND = 8) :: btz, ele3, omega3
    tg1 = -g3*(dtanh(btz*omega3/2.0d0))
    !print*, tg1
    !print*, omega4
    !print *, dtanh(100.0d0),omega3*btz, 'q4'
end 

real (KIND = 8) function tg2(btz2,omega4)
    use cubic
    IMPLICIT NONE 
    real (KIND = 8) :: btz2, ele4, omega4
    tg2 = g4*(dtanh(btz2*omega4/2.0d0))
   ! print *, (dtanh(btz2*omega4/2.0d0)), 'qe'
end 

real (KIND = 8) function tg3(btz3,omega6)
    use cubic
    IMPLICIT NONE 
    real (KIND = 8) :: btz3, ele6, omega6
    tg3 = g6*(dtanh(btz3*omega6/2.0d0))
   ! print *, (dtanh(btz3*omega6/2.0d0)), 'q2'
end 

real (KIND = 8) function fermi(cabete,muq,omega)
    use cubic
    IMPLICIT NONE 
    real (KIND = 8) :: cabete, muq, omega,difarg
    difarg = -cabete*(omega - muq)
    if (dabs(difarg) .le. 1d-6) then 
      difarg = 0.0d0  
    else if (dabs(difarg) .gt. 1d6 .and. dabs(difarg)/difarg .lt. 0.0d0) then 
      fermi = 0.5d0 
    else 
      fermi = 1.0d0/(dexp(difarg) + 1.0d0)
    end if
end 