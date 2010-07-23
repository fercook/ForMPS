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
!    along with FortranMPS.  If not, see <http://www.gnu.org/licenses/>.


module Operator_Class

  use ErrorHandling
  use Constants
  use Tensor_Class

  implicit none
  private

  public :: PauliSigma

  integer,parameter :: SPINONEHALF = 2
  integer,parameter :: OPERATOR_BASIS_SIZE = 3
  integer,parameter :: I_Operator=0,X_Operator=1,Y_Operator=2,Z_Operator=3

!###############################
!#####  The class main object
!###############################
  type,public,extends(Tensor2) :: SpinOperator
     private
  end type SpinOperator


  complex(8),parameter,dimension(SPINONEHALF,SPINONEHALF):: PauliSigmaXMatrix = Reshape([zero,one,one,zero],  [SPINONEHALF,SPINONEHALF])
  complex(8),parameter,dimension(SPINONEHALF,SPINONEHALF):: PauliSigmaYMatrix = Reshape([zero,II,-II,zero],   [SPINONEHALF,SPINONEHALF])
  complex(8),parameter,dimension(SPINONEHALF,SPINONEHALF):: PauliSigmaZMatrix = Reshape([one,zero,zero,-one], [SPINONEHALF,SPINONEHALF])
  complex(8),parameter,dimension(SPINONEHALF,SPINONEHALF):: IdentityOperatorMatrix = Reshape([one,zero,zero,one], [SPINONEHALF,SPINONEHALF])

contains

    function PauliSigma(whichOperator) result(this)
        integer,intent(IN) :: whichOperator
        type(SpinOperator) :: this

        select case (whichOperator)
            case(I_Operator)
                this=new_Tensor(IdentityOperatorMatrix)
            case(X_Operator)
                this=new_Tensor(PauliSigmaXMatrix)
            case(Y_Operator)
                this=new_Tensor(PauliSigmaYMatrix)
            case(Z_Operator)
                this=new_Tensor(PauliSigmaZMatrix)
            case default
                call ThrowException('PauliSigma','Unknown operator requested',whichOperator,CriticalError)
        end select

    end function PauliSigma

 end module Operator_Class
