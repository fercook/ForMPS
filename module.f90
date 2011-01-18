

module Tensor_Class

  implicit none

  complex(8),parameter :: II=(0.0d0,1.0d0)

  type,private :: Tensor
  	private
  	integer :: Initialized=.false.
  contains
    procedure,public :: getDimensions => getDimensions_Of_Tensor
  end type Tensor

  type,public,extends(Tensor) :: Tensor2
  	private
  	complex(8),allocatable :: data(:,:)
  contains
  end type Tensor2

  type,public,extends(Tensor) :: Tensor3
  	private
  	complex(8),allocatable :: data(:,:,:)
  contains
    procedure,public :: Foo => FunctionForTensor3
  end type Tensor3

!###############################
!#####  Operators and methods
!###############################

  interface new_Tensor
     module procedure new_Tensor2,new_Tensor3
  end interface

  interface assignment (=)
     module procedure  new_Tensor3_fromAssignment
  end interface

 contains

!######################################################################################
!#####                           Creation operators
!######################################################################################

   function new_Tensor2 (dim1,dim2) result (this)
     integer,intent(in) :: dim1,dim2
     type(Tensor2) :: this
     real(8) :: randomtensorR(dim1,dim2),randomtensorC(dim1,dim2)

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor2

   function new_Tensor3 (dim1,dim2,dim3) result (this)
     integer,intent(in) :: dim1,dim2,dim3
     type(Tensor3) :: this
     real(8) :: randomtensorR(dim1,dim2,dim3),randomtensorC(dim1,dim2,dim3)

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor3


 !##################################################################

   subroutine new_Tensor3_fromAssignment(lhs,rhs)
     class(Tensor3),intent(out) :: lhs
     type(Tensor3),intent(in) :: rhs

     if(lhs%Initialized) deallocate(lhs%data)
     allocate(lhs%data(size(rhs%data,1),size(rhs%data,2),size(rhs%data,3)))
     lhs%data=rhs%data
     lhs%Initialized=.true.

   end subroutine new_Tensor3_fromAssignment

!##################################################################
!#################   Polymorphic function         ################
!##################################################################

    function getDimensions_Of_Tensor(this) result(Dims)
        class(Tensor),intent(IN) :: this
!       type(Tensor1) :: Dims
        integer,allocatable :: Dims(:)

        select type (Typed_this => this)
            class is (Tensor2)
                allocate(Dims(2))
                Dims=shape(Typed_this%data)
            class is (Tensor3)
                allocate(Dims(3))
                Dims=shape(Typed_this%data)
            class is (Tensor)
                print *,'Dimensions not defined'
                return
        end select

   end function getDimensions_Of_Tensor

!!##################################################################
!!############3####   Type specific functions    ###################
!!##################################################################

function FunctionForTensor3(this,firstindex,secondindex) result (aTensor)
	integer,intent(IN) :: firstindex(:),secondindex(:)
	class(tensor3),intent(IN) :: this
	type(tensor2) :: aTensor

		aTensor=new_Tensor(10,10)

end function FunctionForTensor3


end module Tensor_Class
