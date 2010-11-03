! TAKEN FROM
!
!   http://ftp.cac.psu.edu/pub/ger/fortran/hdk/ginv.f90
!
!
!    This file is part of ForMPS
!
!    ForMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ForMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ForMPS.  If not, see <http://www.gnu.org/licenses/>.

module CACM66InverseCalculator

 program TestGinv
! This file: http://ftp.aset.psu.edu/pub/ger/fortran/hdk/ginv.f90
!
! Compiles/runs with F=Compiler
! use QUADRUPLE_PRECISION
  implicit none

! type (QUAD), Dimension(2,2) :: A,U
! type (QUAD), Dimenision(2) :: AFLAG,ATEMP
  integer, parameter :: DP = selected_real_kind(15)
  real(kind=DP), dimension (3,2) :: A,U, C
  real(kind=DP), dimension (2,2) :: D
  real(kind=DP), dimension (3) :: AFLAG,ATEMP
      integer :: NR,NC, I,J,K

! Simple test case.
      NR=3
      NC=2
      A(1,1)=2.0
      A(2,2)=2.0
      A(1,2)=1.0
      A(2,1)=5.0
      A(3,1)=2.0
      A(3,2)=1.0
      C = A
      print *, "Input Matrix:"
      DO I=1,NR
       print *, (A(I,J),J=1,NC)
      END DO
      print *," "
      print *, "Generalized Inverse should be:"
      print *, " -1.0    1.0   -1.0"
      print *, "  2.5   -2.0    2.5"

      CALL GINV2(A,U,AFLAG,ATEMP,NR,NC)
      print *, " "
      print *, "Generalized Inverse is:"
      DO I=1,NC
       print *, (A(J,I),J=1,NR)
      END DO

!  Show that GINV(A)*A = I
      print *, " "
      print *, "GINV(A)*A=I of Order ",NC,":"
      DO I=1,NC
        DO J=1,NC
          D(I,J) = 0.0
          DO K=1,NR
            D(I,J)=D(I,J) + A(K,I)*C(K,J)
          END DO
        END DO
        PRINT *, (D(I,J), J=1,NC)
      END DO

      contains

   subroutine GINV2 (A,U,AFLAG,ATEMP,NR,NC)

! A Simple Algorithm for Computing the Generalized Inverse of a Matrix
!  by B. Rust, W. R. Burrus and C. Schneeberger
!  CACM 9(5):381-387 (May, 1966)
!
! This routine calculstes the Generalized Inverse of input matrix, A,
! and stores the transpose of it in matrix, A.
! NR -> Number of rows of matrix,  A
! NC -> Number of columns of matrix,  A
! U  -> a bookkeeping matrix.
! AFLAG and ATEMP are temporary working vectors
!
! Notes: If the columns of A are independent, then the Generalized
!        Inverse of A is the Least Squares Inverse of A. That is,
!        GINV can be used to compute Least Squares Regression
!        Coefficients.
!
!        If the matrix A is square with independent columns, then
!        the Generalized Inverse of A is the Inverse of A.
!
!
! type (QUAD), Dimension(:,:) :: A,U
! type (QUAD), Dimenision (:) :: AFLAG,ATEMP
! type (QUAD) :: FAC, TOL, DOT1, DOT2
  real(kind=DP), dimension(:,:), intent(in out) :: A
  real(kind=DP), dimension(:,:), intent(out) :: U
  real(kind=DP), dimension(:), intent(out) :: AFLAG, ATEMP
  real(kind=DP) :: FAC, TOL, DOT1, DOT2
  integer, intent(in) :: NR,NC
  integer  :: I,J,K,L,JM1

      DO  I = 1,NC
        DO J = 1,NC
         U (I,J) = 0.0
        END DO
        U(I,I) = 1.0
      END DO
      FAC = DOT(NR,A,1,1)
      FAC= 1.0_DP/SQRT(FAC)
      DO I = 1,NR
        A(I,1) = A(I,1) * FAC
      END DO
      DO I = 1,NC
        U (I,1) = U(I,1)*FAC
      END DO
      AFLAG(1) = 1.0
!
! Dependent column tolerance, TOL
!
!     N = 27
!     TOL = (10.0 * 0.5**N)**2
      TOL=10.0*EPSILON(FAC)
      DO J = 2,NC
        DOT1 = DOT(NR,A,J,J)
        JM1=J-1
         DO L=1,2
          DO K=1,JM1
            ATEMP(K) = DOT(NR,A,J,K)
         END DO
         DO K=1,JM1
            DO I = 1,NR
              A(I,J) = A(I,J)-ATEMP(K)*A(I,K)*AFLAG(K)
            END DO
            DO I = 1,NC
              U(I,J) = U(I,J)-ATEMP(K)*U(I,K)
            END DO
         END DO
         END DO
        DOT2 = DOT(NR,A,J,J)
        IF((DOT2/DOT1) <= TOL) THEN
          DO I=1,JM1
            ATEMP (I)=0.0
            DO  K=1,I
             ATEMP(I) = ATEMP(I) + U(K,I)*U(K,J)
            END DO
          END DO
          DO I = 1,NR
            A( I,J)=0.0
            DO K=I,JM1
              A(I,J) = A(I,J) - A (I,K)*ATEMP(K)*AFLAG(K)
            END DO
          END DO
          AFLAG(J) = 0.0
          FAC = DOT(NC,U,J,J)
          FAC= 1.0_DP/SQRT(FAC)
        ELSE
          AFLAG(J) = 1.0
          FAC=1.0_DP/SQRT(DOT2)
        ENDIF
        DO I = 1,NR
          A(I,J) = A(I,J)*FAC
        END DO
        DO I = 1,NC
          U(I,J) = U(I,J)*FAC
        END DO
      END DO
      DO J=1,NC
        DO I=1,NR
        FAC = 0.0
        DO K = J,NC
          FAC=FAC+A(I,K)*U(J,K)
        END DO
          A(I,J) = FAC
        END DO
      END DO
      RETURN
   end subroutine Ginv2
      FUNCTION DOT (NR,A,JC,KC) result (PROD)
!  Computes the inner product of columns JC and KC
!  of matrix, A.

! type (QUAD), Dimension(:,:) :: A
! type (QUAD) :: DOT
  real(kind=DP), dimension(:,:), intent(in) :: A
  real(kind=DP) :: PROD
  integer :: I
  integer, intent(in) :: NR,JC,KC
      PROD=0.0
      DO I = 1,NR
        PROD = PROD + A(I,JC)*A(I,KC)
      END DO
      RETURN
   end function DOT
 end program TestGinv

end module CACM66InverseCalculator
