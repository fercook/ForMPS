!!   Copyright 2010 Fernando M. Cucchietti
!
!    This file is part of FortranMPS
!
!    FortranMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    FortranMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    Module Matrix_Helper


  use ErrorHandling
  use Constants

  implicit none


contains

   function vecmat(vector,matrix) result(this)
     real(8),intent(IN) :: vector(:)
     complex(8),intent(IN) :: matrix(:,:)
     complex(8) :: this(size(matrix,1),size(matrix,2))
     integer :: LengthOfVector,LeftDimension,RightDimension
     integer :: i,j

     LengthOfVector=size(vector,1)
     LeftDimension=size(matrix,1)
     RightDimension=size(matrix,2)

     this=zero
     do i=1,min(LeftDimension,LengthOfVector)
        this(i,:)=vector(i)*matrix(i,:)
     enddo
    
   end function vecmat

   function matvec(matrix,vector) result(this)
     real(8),intent(IN) :: vector(:)
     complex(8),intent(IN) :: matrix(:,:)
     complex(8) :: this(size(matrix,1),size(matrix,2))
     integer :: LengthOfVector,LeftDimension,RightDimension
     integer :: i,j

     LengthOfVector=size(vector,1)
     LeftDimension=size(matrix,1)
     RightDimension=size(matrix,2)

     this=zero
     do i=1,min(RightDimension,LengthOfVector)
        this(:,i)=vector(i)*matrix(:,i)
     enddo
    
   end function matvec

   !Simplified interface for LAPACK's ZGESDD routine
   integer function SingularValueDecomposition(matrix,U,Sigma,vTransposed) result(ErrorCode)
     complex(8),intent(IN) :: matrix(:,:)
     complex(8),intent(OUT) :: U(:,:),vTransposed(:,:)
     real(8),intent(OUT) :: Sigma(:)
     integer :: LeftDimension,RightDimension
     !Lapack ugly variables
     integer :: Lwork,LRWork,LIWork,info
     complex(8),allocatable :: Work(:)
     real(8),allocatable :: RWork(:)
     integer(8),allocatable :: IWork(:)
     character,parameter :: Jobz='S' !Always get the minimum only, hopefully the rest of the matrix is zeroed out

     LeftDimension=size(matrix,1); RightDimension=size(matrix,2)

     !Checks
     if( (size(U,1).ne.LeftDimension).or.(size(U,2).ne.LeftDimension).or. &
          & (size(vTransposed,1).ne.RightDimension).or.(size(vTransposed,2).ne.RightDimension).or. &
          & (size(Sigma).ne.Min(LeftDimension,RightDimension)) ) then
        call ThrowException('SingularValueDecomposition','Dimensions of matrices do not match',ErrorCode,CriticalError)
        return
     endif        

     !Recommended values of memory allocation from LAPACK documentation
     LWork=(Min(LeftDimension,RightDimension)*(Min(LeftDimension,RightDimension)+2)+Max(LeftDimension,RightDimension))
     LRWork=5*Min(LeftDimension,RightDimension)*(Min(LeftDimension,RightDimension)+1)
     LIWork=8*Min(LeftDimension,RightDimension)

     allocate(Work(LWork),RWork(LRWork),IWork(LIWork),STAT=ErrorCode)
     If (ErrorCode.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Could not allocate memory',ErrorCode,CriticalError)
        return
     endif
     !For some reason I need to call LAPACK with LWork=-1 first
     !And find out the optimum work storage, otherwise it returns an error
     LWork=-1
     call ZGESDD(JOBZ, LeftDimension, RightDimension, matrix, LeftDimension, Sigma, U, LeftDimension, vTransposed, & 
          & RightDimension,WORK,LWORK,RWORK,IWORK,ErrorCode )
     If (ErrorCode.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Lapack returned error in ZGESDD',ErrorCode,CriticalError)
        return
     endif
     !And now call with right value of LWork
     LWork=Int(Work(1))
     deallocate(Work)
     Allocate(Work(LWork))
     call ZGESDD(JOBZ, LeftDimension, RightDimension, matrix, LeftDimension, Sigma, U, LeftDimension, vTransposed, & 
          & RightDimension,WORK,LWORK,RWORK,IWORK,ErrorCode )
     If (ErrorCode.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Lapack returned error in ZGESDD',ErrorCode,CriticalError)
        return
     endif

     !Clean up
     deallocate(Work,RWork,IWork,STAT=ErrorCode)
     If (ErrorCode.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Problems in deallocation',ErrorCode,CriticalError)
        return
     endif

     ErrorCode=Normal
     
   end function SingularValueDecomposition
     
   real function Difference_btw_Matrices(matrix1, matrix2) result(diff)
     complex(8) :: matrix1(:,:),matrix2(:,:)
     integer :: n,m,alpha,beta

     alpha=size(matrix1,1)
     beta=size(matrix1,2)
     diff=0.0d0
     if(alpha.eq.size(matrix2,1).and.beta.eq.size(matrix2,2)) then
        do n=1,alpha
           do m=1,beta
              diff=diff+(abs(matrix1(n,m)-matrix2(n,m)))**2
           enddo
        enddo
     else
        call ThrowException('Difference_btw_Matrices','Matrices of different shape',NoErrorCode,CriticalError)
     endif
     diff=sqrt(diff)
     return 

   end function Difference_btw_Matrices

end module Matrix_Helper




!This code left here because I don't trust Git :)
!######################################################################################
!!$   subroutine mymatmul(A,B,C,indexL,indexC,indexR,mode)
!!$     complex(8) :: A(:,:),B(:,:),C(:,:)
!!$     integer :: indexL,indexC,indexR
!!$     character*1 :: mode
!!$     integer :: I,J,K,L
!!$     complex(8) TEMP
!!$     ! mode = 'N' is normal multiplication C = A * B + C
!!$     ! mode = 'A' is with A dagged, C = A^+ * B + C
!!$     ! mode = 'B' is with B dagged, C = A * B^+ + C
!!$     if (mode.eq.'N'.or.mode.eq.'n') then
!!$        !C = A * B + C
!!$        DO J = 1,indexR
!!$           DO L = 1,indexC
!!$              IF (B(L,J).NE.ZERO) THEN
!!$                 TEMP = B(L,J)
!!$                 DO I = 1,indexL
!!$                    C(I,J) = C(I,J) + A(I,L)*TEMP
!!$                 enddo
!!$              END IF
!!$           enddo
!!$        enddo
!!$     else if (mode.eq.'A'.or.mode.eq.'a') then
!!$        ! C = A^+ * B + C
!!$        DO J = 1, indexR
!!$           DO I = 1,indexL
!!$              TEMP = ZERO
!!$              DO L = 1,indexC
!!$                 TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
!!$              enddo
!!$             C(I,J) = TEMP + C(I,J)
!!$           enddo
!!$        enddo
!!$     else if (mode.eq.'B'.or.mode.eq.'b') then
!!$        ! C = A * B^+ + C
!!$        DO J = 1,indexR
!!$           DO L = 1,indexC
!!$              IF (B(J,L).NE.ZERO) THEN
!!$                 TEMP = DCONJG(B(J,L))
!!$                 DO I = 1,indexL
!!$                    C(I,J) = C(I,J) + A(I,L)*TEMP
!!$                 enddo
!!$              END IF
!!$           enddo
!!$        enddo
!!$
!!$     endif
!!$   end subroutine mymatmul
