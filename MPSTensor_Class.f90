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


!!  This module contains a class of objects with three legs as used by MPS algorithms.
!!  One leg is "special", it is the physical leg of the tensor.
!!  In notation it is nice to have the spin as the first index, but in the implementation 
!!  it is the third, as it makes it better for performance (this hasn't been checked, just intuition)

module MPSTensor_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class

  implicit none

  private

  public :: new_MPSTensor,LeftCanonize,RightCanonize
  public :: MPSTensor_times_matrix, matrix_times_MPSTensor !operator(.times.)

!###############################
!#####  The class main object
!###############################
  type,public,extends(Tensor3) :: MPSTensor
     private
     integer :: spin,DLeft,DRight
     integer :: Canonized=NO
   contains
     procedure,public :: getDRight => DRight_MPSTensor
     procedure,public :: getDLeft => DLeft_MPSTensor
     procedure,public :: getSpin => Spin_MPSTensor
     procedure,public :: CollapseSpinWithBond => Collapse_Spin_With_Bond_Dimension
!     procedure,public :: delete => delete_MPSTensor
!     procedure,public :: LCanonize => Left_Canonize_MPSTensor
!     procedure,public :: RCanonize => Right_Canonize_MPSTensor
     procedure,public :: ApplyOperator => Apply_Operator_To_Spin_Dimension
  end type MPSTensor

!###############################
!#####  Operators and methods
!###############################
  interface new_MPSTensor
     module procedure new_MPSTensor_Random,new_MPSTensor_fromMPSTensor, & !new_MPSTensor_withData, &
          & new_MPSTensor_withConstant, new_MPSTensor_with_SplitData
  end interface

  interface assignment (=)
     module procedure new_MPSTensor_fromAssignment
  end interface

  interface operator (.times.)
     module procedure MPSTensor_times_matrix, matrix_times_MPSTensor
  end interface

  interface LeftCanonize
    module procedure Left_Canonize_MPSTensor
  end interface

  interface RightCanonize
    module procedure Right_Canonize_MPSTensor
  end interface

  interface SplitSpinFromBond
    module procedure Split_Spin_From_Bond_Dimension
  end interface

!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################

 contains

!######################################################################################
!#####                           Creation operators
!######################################################################################
   function new_MPSTensor_Random (spin,DLeft,DRight) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     type(MPSTensor) :: this

     this=new_Tensor(DLeft,DRight,spin)
     !initialize internal variables
     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight

   end function new_MPSTensor_Random

!##################################################################
!   function new_MPSTensor_withData (originalData) result (this)
!     complex(8),intent(in) :: originalData(:,:,:)
!     type(MPSTensor) this
!
!     this=new_Tensor(originalData)
!     !initialize internal variables
!     this%spin=size(originalData,3)
!     this%DLeft=size(originalData,1)
!     this%DRight=size(originalData,2)
!
!   end function new_MPSTensor_withData

!##################################################################
   function new_MPSTensor_with_SplitData (originalDataUP,originalDataDOWN) result (this)
     complex(8),intent(in) :: originalDataUP(:,:),originalDataDOWN(:,:)
     complex(8),allocatable :: JoinedData(:,:,:)
     type(MPSTensor) this
     integer :: upDims(2),downDims(2)

     upDims=shape(originalDataUP)
     downDims=shape(originalDataDOWN)

     if (upDims.equalvector.downDims) then
        allocate (JoinedData(upDims(1),upDims(2),2))
        JoinedData(:,:,1)=originalDataUP
        JoinedData(:,:,2)=originalDataDOWN

        this=new_Tensor(JoinedData)
        !initialize internal variables
        this%spin=2
        this%DLeft=size(originalDataUP,1)
        this%DRight=size(originalDataUP,2)

      else
        call ThrowException('new_MPSTensor_with_SplitData','Up and down data have different size',upDims(1)-downDims(2),CriticalError)
      endif
   end function new_MPSTensor_with_SplitData

!##################################################################
   function new_MPSTensor_withConstant (spin,DLeft,DRight,constant) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     complex(8),intent(in) :: constant
     type(MPSTensor) :: this

     this=new_Tensor(DLeft,DRight,spin,constant)
     !initialize internal variables
     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight

   end function new_MPSTensor_withConstant

!##################################################################
   function new_MPSTensor_fromMPSTensor (tensor) result (this)
     class(MPSTensor),intent(in) :: tensor
     type(MPSTensor) this

     this=new_Tensor(tensor)
     !initialize internal variables
     this%spin=tensor%spin
     this%DLeft=tensor%DLeft
     this%DRight=tensor%DRight

   end function new_MPSTensor_fromMPSTensor

!##################################################################
   subroutine new_MPSTensor_fromAssignment(lhs,rhs)
     class(MPSTensor),intent(out) :: lhs
     type(MPSTensor),intent(in) :: rhs

     lhs=new_Tensor(rhs)

     lhs%spin=rhs%spin
     lhs%DLeft=rhs%DLeft
     lhs%DRight=rhs%DRight

   end subroutine new_MPSTensor_fromAssignment


!##################################################################
!###########       Accessor methods
!##################################################################
   integer function spin_MPSTensor(this) result(s)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
 
    if(.not.(this%IsInitialized())) then
        call ThrowException('Spin','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        s=this%spin
     endif

   end function spin_MPSTensor
!##################################################################

   integer function DLeft_MPSTensor(this) result(DL)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DLeft','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DL=this%DLeft
     endif

   end function DLeft_MPSTensor
!##################################################################
   integer function DRight_MPSTensor(this) result(DR)

     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DR=this%DRight
     endif

   end function DRight_MPSTensor

!##################################################################
!#######################        Products by things
!##################################################################

   function Apply_Operator_To_Spin_Dimension(this,anOperator) result(aTensor)
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     type(MPSTensor) :: aTensor
     class(SpinOperator),intent(IN) :: anOperator
     integer :: opDims(2),tensorDims(3)

     opDims=anOperator%GetDimensions()
     if(this%IsInitialized()) then
        if(opDims(1).eq.this%spin.and.opDims(2).eq.this%spin) then
           aTensor = new_Tensor(this*anOperator)
           tensorDims=aTensor%GetDimensions()
           aTensor%spin=tensorDims(3)
           aTensor%DLeft=tensorDims(1)
           aTensor%DRight=tensorDims(2)
        else
           call ThrowException('Apply_Operator_To_Spin_Dimension','Operator is not of the rigt size',opDims(1)-opDims(2),CriticalError)
        endif
     else
        call ThrowException('Apply_Operator_To_Spin_Dimension','Tensor not initialized',NoErrorCode,CriticalError)
     endif

   end function Apply_Operator_To_Spin_Dimension
!##################################################################
    function MPSTensor_times_matrix(aTensor,aMatrix) result(theResult)
        class(MPSTensor),intent(IN) :: aTensor
        class(Tensor2),intent(IN) :: aMatrix
        type(MPSTensor) :: theResult
        integer :: newDims(3)

        theResult=TensorTranspose( TensorTranspose(aTensor,[1,3,2])*aMatrix , [1,3,2] )
        newDims=theResult%GetDimensions()

        theResult%spin=newDims(3)
        theResult%DLeft=newDims(1)
        theResult%DRight=newDims(2)

    end function MPSTensor_times_matrix
!##################################################################
    function matrix_times_MPSTensor(aMatrix,aTensor) result(theResult)
        class(MPSTensor),intent(IN) :: aTensor
        class(Tensor2),intent(IN) :: aMatrix
        type(MPSTensor) :: theResult
        integer :: newDims(3)

        theResult=aMatrix*aTensor
        newDims=theResult%GetDimensions()

        theResult%spin=newDims(3)
        theResult%DLeft=newDims(1)
        theResult%DRight=newDims(2)

    end function matrix_times_MPSTensor


!##################################################################
!##################################################################
!##################################################################
!##################################################################

    function Collapse_Spin_With_Bond_Dimension(this,whichDimension) result(collapsedTensor)
        class(MPSTensor),intent(IN) :: this
        integer,intent(IN) :: whichDimension(1)
        type(Tensor2) :: collapsedTensor

        if(this%IsInitialized()) then
            select case (whichDimension(1))
                case (FIRST(1))
                    collapsedTensor=this%JoinIndices(THIRDANDFIRST,SECOND)
                case (SECOND(1))
                    collapsedTensor=this%JoinIndices(THIRDANDSECOND,FIRST)
                case default
                    call ThrowException('Collapse_Spin_With_Bond_Dimension','Dimension must be FIRST or SECOND',whichDimension(1),CriticalError)
            end select
        else
            call ThrowException('Collapse_Spin_With_Bond_Dimension','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function Collapse_Spin_With_Bond_Dimension

!##################################################################

    function Split_Spin_From_Bond_Dimension(this,whichDimension,spinSize) result(splitTensor)
        class(Tensor2),intent(IN) :: this
        integer,intent(IN) :: whichDimension(1),spinSize
        type(MPSTensor) :: splitTensor
        integer :: newDims(3)

        if(this%IsInitialized()) then
            select case (whichDimension(1))
                case (FIRST(1))
                    splitTensor=TensorTranspose( SplitIndexOf(this,FIRST,spinSize), [3,1,2] )
                    newDims=splitTensor%GetDimensions()
                    splitTensor%spin=newDims(3)
                    splitTensor%DLeft=newDims(1)
                    splitTensor%DRight=newDims(2)
                case (SECOND(1))
                    splitTensor=TensorTranspose( SplitIndexOf(this,SECOND,spinSize), [1,3,2] )
                    newDims=splitTensor%GetDimensions()
                    splitTensor%spin=newDims(3)
                    splitTensor%DLeft=newDims(1)
                    splitTensor%DRight=newDims(2)
                case default
                    call ThrowException('Split_Spin_From_Bond_Dimension','Dimension must be FIRST or SECOND',whichDimension(1),CriticalError)
            end select
        else
            call ThrowException('Split_Spin_From_Bond_Dimension','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function Split_Spin_From_Bond_Dimension

!##################################################################
!##################################################################
! Right Site Canonization -- Returns the matrix that needs to be multiplied
! to the adjacent site on the RIGHT
!##################################################################
!##################################################################

  function Right_Canonize_MPSTensor(this) result(matrix)
    class(MPSTensor),intent(INOUT) :: this !
    type(Tensor2) :: matrix
    type(Tensor2) :: U,Sigma,vTransposed,collapsedTensor
    integer :: error,newUDims(2)

    if(.not.this%IsInitialized()) then
       call ThrowException('Left_Canonize_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    collapsedTensor=this%CollapseSpinWithBond(FIRST)

    if (WasThereError()) then
       call ThrowException('Left_Canonize_MPSTensor','Could not collapse the tensor',NoErrorCode,CriticalError)
       return
    endif

    call SingularValueDecomposition(CollapsedTensor,U,Sigma,vTransposed,error)

    if (error.ne.Normal) then
       call ThrowException('Left_Canonize_MPSTensor','Could not split the matrix',NoErrorCode,CriticalError)
       return
    endif

    !We need to trim the matrix to get rid of useless dimensions
    newUDims=[ this%spin*this%DLeft, min(this%spin * this%DLeft, this%Dright)]
    this=SplitSpinFromBond( TensorPad(U,newUDims), FIRST, this%spin )

    !matrix is reshaped to fit the product with the tensor on the right
    matrix=TensorPad(sigma*Conjugate(vTransposed), [ newUDims(2), this%DRight ] )

  end function Right_Canonize_MPSTensor

!##################################################################
!##################################################################
! Left Site Canonization -- Returns the matrix that needs to be multiplied
! to the adjacent site on the LEFT
!##################################################################
!##################################################################

  function Left_Canonize_MPSTensor(this) result(matrix)
    class(MPSTensor),intent(INOUT) :: this !
    type(Tensor2) :: matrix
    type(Tensor2) :: U,Sigma,vTransposed,collapsedTensor
    integer :: error,newVDims(2)

    if(.not.this%IsInitialized()) then
       call ThrowException('Left_Canonize_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    collapsedTensor=this%CollapseSpinWithBond(SECOND)
    if (WasThereError()) then
       call ThrowException('Left_Canonize_MPSTensor','Could not collapse the tensor',NoErrorCode,CriticalError)
       return
    endif

    call SingularValueDecomposition(CollapsedTensor,U,Sigma,vTransposed,error)
    if (error.ne.Normal) then
       call ThrowException('Left_Canonize_MPSTensor','Could not split the matrix',NoErrorCode,CriticalError)
       return
    endif

    newVDims=[ min(this%spin * this%DRight, this%DLeft), this%spin*this%DRight]
    this=SplitSpinFromBond( TensorPad(Conjugate(vTransposed),newVDims), FIRST, this%spin )

    !matrix is reshaped to fit the product with the tensor on the right
    matrix=TensorPad( U*sigma, [ this%DLeft, newVDims(1) ] )

  end function Left_Canonize_MPSTensor

!#######################################################################################
!#######################################################################################


 end module MPSTensor_Class
