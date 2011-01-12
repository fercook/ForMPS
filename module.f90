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


module Tensor_Class

  use ErrorHandling
  use Constants
!  use lapack95
!  use NNLSInverseCalculator
!  use CACM66InverseCalculator

  implicit none
! Need to learn how to make operators public.

  integer,parameter :: Max_Combined_Dimension = 10000000

  type,private :: Tensor
  	private
  	integer :: Initialized=.false.
  contains
  	procedure,public :: IsInitialized => Is_Tensor_init !Commented out because of Ifort bug
  	procedure,public :: print => print_Tensor
	procedure,public :: PrintDimensions => Print_Tensor_Dimensions
    procedure,public :: getDimensions => getDimensions_Of_Tensor
    procedure,public :: Norm => Norm_Of_Tensor
  end type Tensor

  type,public,extends(Tensor) :: Tensor1
  	private
  	complex(8),allocatable :: data(:)
  contains
    procedure,public :: Slice => Take_Slice_Of_Tensor1
    final :: delete_Tensor1
  end type Tensor1

  type,public,extends(Tensor) :: Tensor2
  	private
  	complex(8),allocatable :: data(:,:)
  contains
    procedure,public :: SVD => SingularValueDecompositionTensor2
    procedure,public :: SplitIndex => SplitIndexOfTensor2
    procedure,public :: Slice => Take_Slice_Of_Tensor2
    procedure,public :: dagger => ConjugateTranspose2
    procedure,public :: trace => Tensor2Trace
    procedure,public :: CompactFromLeft => Mirror_Compact_Left_With_Tensor3
    procedure,public :: CompactFromRight => Mirror_Compact_Right_With_Tensor3
    final :: delete_Tensor2
  end type Tensor2

  type,public,extends(Tensor) :: Tensor3
  	private
  	complex(8),allocatable :: data(:,:,:)
  contains
    procedure,public :: JoinIndices => JoinIndicesOfTensor3
    procedure,public :: SplitIndex => SplitIndexOfTensor3
    procedure,public :: CompactFromBelow => Compact_From_Below_With_Tensor4
    procedure,public :: Slice => Take_Slice_Of_Tensor3
    procedure,public :: PartialTrace => Tensor3Trace
    procedure,public :: Unfold => UnfoldTensor3
    final :: delete_Tensor3
  end type Tensor3

  type,public,extends(Tensor) :: Tensor4
    private
    complex(8),allocatable :: data(:,:,:,:)
  contains
    procedure,public :: JoinIndices => JoinIndicesOfTensor4
    procedure,public :: Slice => Take_Slice_Of_Tensor4
    procedure,public :: PartialTrace => Tensor4Trace
    procedure,public :: Unfold => UnfoldTensor4
    final :: delete_Tensor4
  end type Tensor4



!###############################
!#####  Operators and methods
!###############################

  interface new_Tensor
     module procedure new_Tensor1_Random,new_Tensor1_fromTensor1,new_Tensor1_fromData,new_Tensor1_withConstant, &
		& new_Tensor2_Random,new_Tensor2_fromTensor2,new_Tensor2_fromData,new_Tensor2_withConstant, &
		& new_Tensor3_Random,new_Tensor3_fromTensor3,new_Tensor3_fromData,new_Tensor3_withConstant, &
		& new_Tensor4_Random,new_Tensor4_fromTensor4,new_Tensor4_fromData,new_Tensor4_withConstant
  end interface

  interface operator (*)
     module procedure &
     	  & number_times_Tensor1,number_times_Tensor2,number_times_Tensor3,number_times_Tensor4, &
     	  & Tensor2_matmul_Tensor2, Tensor2_matmul_Tensor1, Tensor1_matmul_Tensor2, &
     	  & Tensor1_dotProduct_Tensor1, &
     	  & Tensor3_matmul_Tensor3, Tensor2_matmul_Tensor3, Tensor3_matmul_Tensor2
  end interface

!  interface operator ()
!     module procedure &
!     &   Tensor4_doubletimes_Tensor4,Tensor2_doubletimes_Tensor2
!  end interface

  interface operator (.x.)
     module procedure &
          & Tensor2_matmul_Tensor2, Tensor2_matmul_Tensor1, Tensor1_matmul_Tensor2, &
          & Tensor1_dotProduct_Tensor1, Tensor3_matmul_Tensor3, Tensor2_matmul_Tensor3, &
          & Tensor3_matmul_Tensor2
  end interface

  interface operator (.xxx.)
     module procedure &
     &   Tensor3_tripletimes_Tensor3
  end interface

  interface operator (.xplus.)
     module procedure &
     &   MultAndCollapse_Tensor3_Tensor4, MultAndCollapse_Tensor3_Tensor3, &
     &   MultAndCollapse_Tensor4_Tensor4
  end interface

  interface operator (+)
     module procedure add_Tensor1,add_Tensor2,add_Tensor3,add_Tensor4
  end interface

  interface operator (-)
     module procedure subtract_Tensor1,subtract_Tensor2,subtract_Tensor3,subtract_Tensor4
  end interface

  interface MultAndCollapse
    module procedure MultAndCollapse_Tensor3_Tensor4, &
        &  MultAndCollapse_Tensor3_Tensor3, MultAndCollapse_Tensor4_Tensor4
  end interface

  interface assignment (=)
     module procedure new_Tensor1_fromAssignment, new_Tensor2_fromAssignment, &
          & new_Tensor3_fromAssignment, new_Tensor4_fromAssignment
  end interface

  interface operator (.diff.)
     module procedure Difference_btw_Tensors
  end interface

  interface operator (.absdiff.)
     module procedure Difference_btw_Tensors_WithAbsoluteValue
  end interface

  interface operator (.equaldims.)
     module procedure  Tensors_are_of_equal_Shape
  end interface

  interface operator (.equaltype.)
     module procedure Tensors_are_of_equal_Type
  end interface

  interface operator (.unfold.)
     module procedure UnfoldTensor3,UnfoldTensor4
  end interface

  interface TensorSlice
    module procedure Take_Slice_Of_Tensor1,Take_Slice_Of_Tensor2,Take_Slice_Of_Tensor3, &
        & Take_Slice_Of_Tensor4
  end interface

  interface JoinIndicesOf
  	module procedure JoinIndicesOfTensor3,JoinIndicesOfTensor4 !,JoinTwoIndicesOfTensor4
  end interface

  interface SplitIndexOf
    module procedure SplitIndexOfTensor2,SplitIndexOfTensor3
  end interface

  interface TensorPad
    module procedure Pad_Tensor2
  end interface

  interface Conjugate
    module procedure ConjugateTensor1,ConjugateTensor2,ConjugateTensor3,ConjugateTensor4
  end interface

  interface TensorTranspose
    module procedure TensorTranspose2,TensorTranspose3,TensorTranspose4
  end interface

  interface ConjugateTranspose
    module procedure ConjugateTranspose2,ConjugateTranspose3,ConjugateTranspose4
  end interface

  interface TensorTrace
    module procedure Tensor2Trace,Tensor3Trace,Tensor4Trace
  end interface

  interface CompactLeft
    module procedure Compact_Tensor3_From_Left_With_Tensor2,Mirror_Compact_Left_With_Tensor3
  end interface

  interface CompactRight
    module procedure Compact_Tensor3_From_Right_With_Tensor2,Mirror_Compact_Right_With_Tensor3
  end interface

  interface CompactBelow
    module procedure Compact_From_Below_With_Tensor4
  end interface


  interface SingularValueDecomposition
    module procedure SingularValueDecompositionTensor2
  end interface

  interface SolveLinearProblem
    module procedure SolveLinearProblem_LAPACK_matrix, SolveLinearProblem_LAPACK_vector
    !SolveLinearProblem_Netlib,SolveLinearProblem_NNLSNasa,SolveLinearProblem_CACM66
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
   function new_Tensor1_Random (dim1) result (this)
     integer,intent(in) :: dim1
     type(Tensor1) :: this
     real(8) :: randomtensorR(dim1),randomtensorC(dim1)

     if(dim1.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor1_Random','Dimension is larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1) then
        call ThrowException('new_Tensor1_Random','Dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor1_Random

   function new_Tensor2_Random (dim1,dim2) result (this)
     integer,intent(in) :: dim1,dim2
     type(Tensor2) :: this
     real(8) :: randomtensorR(dim1,dim2),randomtensorC(dim1,dim2)

     if(dim1*dim2.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor2_Random','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1) then
        call ThrowException('new_Tensor2_Random','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor2_Random

   function new_Tensor3_Random (dim1,dim2,dim3) result (this)
     integer,intent(in) :: dim1,dim2,dim3
     type(Tensor3) :: this
     real(8) :: randomtensorR(dim1,dim2,dim3),randomtensorC(dim1,dim2,dim3)

     if(dim1*dim2*dim3.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor3_Random','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1) then
        call ThrowException('new_Tensor3_Random','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor3_Random

   function new_Tensor4_Random (dim1,dim2,dim3,dim4) result (this)
     integer,intent(in) :: dim1,dim2,dim3,dim4
     type(Tensor4) :: this
     real(8) :: randomtensorR(dim1,dim2,dim3,dim4),randomtensorC(dim1,dim2,dim3,dim4)

     if(dim1*dim2*dim3*dim4.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor4_Random','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1.or.dim4.lt.1) then
        call ThrowException('new_Tensor4_Random','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3,dim4))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor4_Random


!##################################################################
   function new_Tensor1_fromData (originalData) result (this)
     complex(8),intent(in) :: originalData(:)
     integer :: dim1
     type(Tensor1) this

     dim1=size(originalData,1)

     if(dim1.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor1_fromData','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1) then
        call ThrowException('new_Tensor1_fromData','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1))

     this%data=originalData
     this%Initialized=.true.

   end function new_Tensor1_fromData

   function new_Tensor2_fromData (originalData) result (this)
     complex(8),intent(in) :: originalData(:,:)
     integer :: dim1,dim2
     type(Tensor2) this

	dim1=size(originalData,1)
	dim2=size(originalData,2)

     if(dim1*dim2.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor2_fromData','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1) then
        call ThrowException('new_Tensor3_fromData','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2))

     this%data=originalData
     this%Initialized=.true.

   end function new_Tensor2_fromData

   function new_Tensor3_fromData (originalData) result (this)
     complex(8),intent(in) :: originalData(:,:,:)
     integer :: dim1,dim2,dim3
     type(Tensor3) this

	dim1=size(originalData,1)
	dim2=size(originalData,2)
	dim3=size(originalData,3)
     if(dim1*dim2*dim3.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor3_fromData','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1) then
        call ThrowException('new_Tensor3_fromData','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))
     this%data=originalData
     this%Initialized=.true.

   end function new_Tensor3_fromData

   function new_Tensor4_fromData (originalData) result (this)
     complex(8),intent(in) :: originalData(:,:,:,:)
     integer :: dim1,dim2,dim3,dim4
     type(Tensor4) this

    dim1=size(originalData,1)
    dim2=size(originalData,2)
    dim3=size(originalData,3)
    dim4=size(originalData,4)
     if(dim1*dim2*dim3*dim4.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor4_fromData','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1.or.dim4.lt.1) then
        call ThrowException('new_Tensor4_fromData','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3,dim4))
     this%data=originalData
     this%Initialized=.true.

   end function new_Tensor4_fromData

!##################################################################
   function new_Tensor1_withConstant (dim1,constant) result (this)
     integer,intent(in) :: dim1
     complex(8),intent(in) :: constant
     type(Tensor1) this

     if(dim1.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor1_withConstant','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1) then
        call ThrowException('new_Tensor1_withConstant','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1))
     this%data=constant
     this%Initialized=.true.

   end function new_Tensor1_withConstant

   function new_Tensor2_withConstant (dim1,dim2,constant) result (this)
     integer,intent(in) :: dim1,dim2
     complex(8),intent(in) :: constant
     type(Tensor2) this

     if(dim1*dim2.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor2_withConstant','Dimensions are larger than maximum',dim1*dim2,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1) then
        call ThrowException('new_Tensor2_withConstant','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2))
     this%data=constant
     this%Initialized=.true.

   end function new_Tensor2_withConstant

   function new_Tensor3_withConstant (dim1,dim2,dim3,constant) result (this)
     integer,intent(in) :: dim1,dim2,dim3
     complex(8),intent(in) :: constant
     type(Tensor3) this

     if(dim1*dim2*dim3.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor3_withConstant','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1) then
        call ThrowException('new_Tensor3_withConstant','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))
     this%data=constant
     this%Initialized=.true.

   end function new_Tensor3_withConstant

   function new_Tensor4_withConstant (dim1,dim2,dim3,dim4,constant) result (this)
     integer,intent(in) :: dim1,dim2,dim3,dim4
     complex(8),intent(in) :: constant
     type(Tensor4) this

     if(dim1*dim2*dim3*dim4.gt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor4_withConstant','Dimensions are larger than maximum',dim1*dim2*dim3*dim4,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1.or.dim4.lt.1) then
        call ThrowException('new_Tensor4_withConstant','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3,dim4))
     this%data=constant
     this%Initialized=.true.

   end function new_Tensor4_withConstant

!##################################################################
   function new_Tensor1_fromTensor1 (tensor) result (this)
     class(Tensor1),intent(in) :: tensor
     type(Tensor1) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor1_fromTensor1

   function new_Tensor2_fromTensor2 (tensor) result (this)
     class(Tensor2),intent(in) :: tensor
     type(Tensor2) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor2_fromTensor2

   function new_Tensor3_fromTensor3 (tensor) result (this)
     class(Tensor3),intent(in) :: tensor
     type(Tensor3) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2),size(tensor%data,3)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor3_fromTensor3

   function new_Tensor4_fromTensor4 (tensor) result (this)
     class(Tensor4),intent(in) :: tensor
     type(Tensor4) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2),size(tensor%data,3),size(tensor%data,4)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor4_fromTensor4



 !##################################################################

   subroutine new_Tensor1_fromAssignment(lhs,rhs)
     class(Tensor1),intent(out) :: lhs
     type(Tensor1),intent(in) :: rhs

     if(lhs%Initialized) deallocate(lhs%data)
     allocate(lhs%data(size(rhs%data,1)))
     lhs%data=rhs%data
     lhs%Initialized=.true.

   end subroutine new_Tensor1_fromAssignment

   subroutine new_Tensor2_fromAssignment(lhs,rhs)
     class(Tensor2),intent(out) :: lhs
     type(Tensor2),intent(in) :: rhs

     if(lhs%Initialized) deallocate(lhs%data)
     allocate(lhs%data(size(rhs%data,1),size(rhs%data,2)))
     lhs%data=rhs%data
     lhs%Initialized=.true.

   end subroutine new_Tensor2_fromAssignment

   subroutine new_Tensor3_fromAssignment(lhs,rhs)
     class(Tensor3),intent(out) :: lhs
     type(Tensor3),intent(in) :: rhs

     if(lhs%Initialized) deallocate(lhs%data)
     allocate(lhs%data(size(rhs%data,1),size(rhs%data,2),size(rhs%data,3)))
     lhs%data=rhs%data
     lhs%Initialized=.true.

   end subroutine new_Tensor3_fromAssignment

   subroutine new_Tensor4_fromAssignment(lhs,rhs)
     class(Tensor4),intent(out) :: lhs
     type(Tensor4),intent(in) :: rhs

     if(lhs%Initialized) deallocate(lhs%data)
     allocate(lhs%data(size(rhs%data,1),size(rhs%data,2),size(rhs%data,3),size(rhs%data,4)))
     lhs%data=rhs%data
     lhs%Initialized=.true.

   end subroutine new_Tensor4_fromAssignment



 !######################################    delete
   subroutine delete_Tensor1 (this)
     type(Tensor1) :: this
	 if (allocated(this%data)) then
	   deallocate(this%data)
     endif
     this%Initialized=.false.
   end subroutine delete_Tensor1
!##################################################################
   subroutine delete_Tensor2 (this)
     type(Tensor2) :: this
     if (allocated(this%data)) then
       deallocate(this%data)
     endif
     this%Initialized=.false.
   end subroutine delete_Tensor2
!##################################################################
   subroutine delete_Tensor3 (this)
     type(Tensor3) :: this
     if (allocated(this%data)) then
       deallocate(this%data)
     endif
     this%Initialized=.false.
   end subroutine delete_Tensor3
!##################################################################
   subroutine delete_Tensor4 (this)
     type(Tensor4) :: this
     if (allocated(this%data)) then
       deallocate(this%data)
     endif
     this%Initialized=.false.
   end subroutine delete_Tensor4
!##################################################################
!##################################################################
!##################################################################

logical function Is_Tensor_init(this) result(AmIInitialized)
    class(Tensor) :: this

    AmIInitialized=this%Initialized

end function Is_Tensor_Init

integer function InitializationCheck(this) result(error)
    class(Tensor),intent(IN) :: this

    if (.not.this%Initialized) then
       error=CriticalError
       call ThrowException('Internal Routine ','Uninitialized tensor',NoErrorCode,error)
    else
       error=Normal
    endif

  end function InitializationCheck


!######################################     print
   subroutine Print_Tensor(this,message,error)
     class(Tensor),intent(IN) :: this
     integer i,j,k
     character*(*),optional :: message
     integer,optional :: error

     If(present(error)) error = Warning

     if(.not.(this%Initialized)) then
        call ThrowException('PrintTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     endif

     if(present(message)) print *,message

	 select type (Typed_this => this)
	 	class is (Tensor1)
        	print *,'Vector data:'
	 		print *,Typed_this%data
	 	class is (Tensor2)
        	print *,'Matrix data:'
	 		print *,Typed_this%data
	 	class is (Tensor3)
        	print *,'3-Tensor data:'
	 		print *,Typed_this%data
        class is (Tensor4)
            print *,'4-Tensor data:'
            print *,Typed_this%data
	 	class is (Tensor)
	 		print *,'No data in raw tensor'
	 end select

     If(present(error)) error=Normal

   end subroutine Print_Tensor


   subroutine Print_Tensor_Dimensions(this,message)
     class(Tensor),intent(IN) :: this  !!<<TYPE>>!!
     character*(*),optional :: message
     integer i,j,k

     if(.not.(this%Initialized)) then
        call ThrowException('Print Tensor','Tensor not initialized',NoErrorCode,Warning)
        return
     endif

     if (present(message)) write(*,'(A)'),message

	 select type (Typed_this => this)
	 	class is (Tensor1)
            write(*,'("Vector Dimension:",I6)'), size(Typed_this%data,1)
	 	class is (Tensor2)
            write(*,'("Matrix Dimensions:",I4," x",I4)'),size(Typed_this%data,1),size(Typed_this%data,2)
	 	class is (Tensor3)
            write(*,'("3-Tensor Dimensions:",I3," x",I3," x",I4)'),size(Typed_this%data,1),size(Typed_this%data,2), &
                  & size(Typed_this%data,3)
        class is (Tensor4)
            write(*,'("4-Tensor Dimensions:",I3," x",I3," x",I3," x",I3)'),size(Typed_this%data,1),size(Typed_this%data,2), &
                  & size(Typed_this%data,3),size(Typed_this%data,4)
	 	class is (Tensor)
	 		print *,'No data in raw tensor'
	 end select

   end subroutine Print_Tensor_Dimensions
!##################################################################


	function getDimensions_Of_Tensor(this) result(Dims)
		class(Tensor),intent(IN) :: this
!		type(Tensor1) :: Dims
		integer,allocatable :: Dims(:)

	    if(.not.(this%Initialized)) then
    	    call ThrowException('getDimensions','Tensor not initialized',NoErrorCode,CriticalError)
        	return
	    endif

	 	select type (Typed_this => this)
	 		class is (Tensor1)
	 		    allocate(Dims(1))
	 		    Dims=shape(Typed_this%data)
	 		class is (Tensor2)
                allocate(Dims(2))
                Dims=shape(Typed_this%data)
	 		class is (Tensor3)
                allocate(Dims(3))
                Dims=shape(Typed_this%data)
            class is (Tensor4)
                allocate(Dims(4))
                Dims=shape(Typed_this%data)
	 		class is (Tensor)
	    	    call ThrowException('getDimensions','Dimensions not defined',NoErrorCode,CriticalError)
    	    	return
	 	end select

   end function getDimensions_Of_Tensor
!##################################################################

   real(8) function Norm_Of_Tensor(this)
     class(Tensor),intent(IN) :: this

	 if(.not.(this%Initialized)) then
   		call ThrowException('Norm_Of_Tensor','Tensor not initialized',NoErrorCode,CriticalError)
        return
	 endif

     Norm_Of_Tensor=0.0d0
	 select type (Typed_this => this)
	 	class is (Tensor1)
	 		Norm_Of_Tensor=sum(abs(Typed_this%data))
	 	class is (Tensor2)
	 		Norm_Of_Tensor=sum(abs(Typed_this%data))
	 	class is (Tensor3)
	 		Norm_Of_Tensor=sum(abs(Typed_this%data))
        class is (Tensor4)
            Norm_Of_Tensor=sum(abs(Typed_this%data))
	 	class is (Tensor)
	    	call ThrowException('Norm_Of_Tensor','Norm is not defined',NoErrorCode,CriticalError)
    	    return
	 end select

     return
   end function Norm_Of_Tensor

!
!! FUTURE POLYMORPHIC CODE
!   function number_times_Tensor(constant, aTensor) result(this)
!     complex(8),intent(IN) :: constant
!     class(Tensor),intent(IN) :: aTensor
!     class(Tensor),allocatable :: this
!
!     if(.not.aTensor%Initialized) then
!        call ThrowException('Number_times_Tensor','Tensor not initialized',NoErrorCode,CriticalError)
!	 else
!	 	 allocate(this,SOURCE=aTensor)
!     endif
!	 return
!
!   end function Number_times_Tensor
!

   function number_times_Tensor1(constant, aTensor) result(this)
     complex(8),intent(IN) :: constant
     class(Tensor1),intent(IN) :: aTensor
     type(Tensor1) :: this

     if(.not.aTensor%Initialized) then
        call ThrowException('Number_times_Tensor','Tensor not initialized',NoErrorCode,CriticalError)
	 else
	 	 this=new_Tensor(constant*aTensor%data)
     endif
	 return

   end function Number_times_Tensor1

   function number_times_Tensor2(constant, aTensor) result(this)
     complex(8),intent(IN) :: constant
     class(Tensor2),intent(IN) :: aTensor
     type(Tensor2) :: this

     if(.not.aTensor%Initialized) then
        call ThrowException('Number_times_Tensor','Tensor not initialized',NoErrorCode,CriticalError)
	 else
	 	 this=new_Tensor(constant*aTensor%data)
     endif
	 return

   end function Number_times_Tensor2

   function number_times_Tensor3(constant, aTensor) result(this)
     complex(8),intent(IN) :: constant
     class(Tensor3),intent(IN) :: aTensor
     type(Tensor3) :: this

     if(.not.aTensor%Initialized) then
        call ThrowException('Number_times_Tensor','Tensor not initialized',NoErrorCode,CriticalError)
	 else
	 	 this=new_Tensor(constant*aTensor%data)
     endif
	 return

   end function Number_times_Tensor3

   function number_times_Tensor4(constant, aTensor) result(this)
     complex(8),intent(IN) :: constant
     class(Tensor4),intent(IN) :: aTensor
     type(Tensor4) :: this

     if(.not.aTensor%Initialized) then
        call ThrowException('Number_times_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
     else
         this=new_Tensor(constant*aTensor%data)
     endif
     return

   end function Number_times_Tensor4


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function add_Tensor1(tensorA,tensorB) result(this)
        class(tensor1),intent(IN) :: tensorA,tensorB
        type(Tensor1) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data+tensorB%data)
        else
           call ThrowException('add_Tensor1','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function add_Tensor1

    function add_Tensor2(tensorA,tensorB) result(this)
        class(tensor2),intent(IN) :: tensorA,tensorB
        type(Tensor2) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data+tensorB%data)
        else
           call ThrowException('add_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function add_Tensor2

    function add_Tensor3(tensorA,tensorB) result(this)
        class(tensor3),intent(IN) :: tensorA,tensorB
        type(Tensor3) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data+tensorB%data)
        else
           call ThrowException('add_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function add_Tensor3

    function add_Tensor4(tensorA,tensorB) result(this)
        class(tensor4),intent(IN) :: tensorA,tensorB
        type(Tensor4) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data+tensorB%data)
        else
           call ThrowException('add_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function add_Tensor4


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function subtract_Tensor1(tensorA,tensorB) result(this)
        class(tensor1),intent(IN) :: tensorA,tensorB
        type(Tensor1) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data-tensorB%data)
        else
           call ThrowException('subtract_Tensor1','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function subtract_Tensor1

    function subtract_Tensor2(tensorA,tensorB) result(this)
        class(tensor2),intent(IN) :: tensorA,tensorB
        type(Tensor2) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data-tensorB%data)
        else
           call ThrowException('subtract_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function subtract_Tensor2

    function subtract_Tensor3(tensorA,tensorB) result(this)
        class(tensor3),intent(IN) :: tensorA,tensorB
        type(Tensor3) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data-tensorB%data)
        else
           call ThrowException('subtract_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function subtract_Tensor3

    function subtract_Tensor4(tensorA,tensorB) result(this)
        class(tensor4),intent(IN) :: tensorA,tensorB
        type(Tensor4) :: this
        if(TensorA%Initialized.and.TensorB%Initialized.and.(tensorA.equaldims.tensorB)) then
           this=new_Tensor(tensorA%data-tensorB%data)
        else
           call ThrowException('subtract_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif
        return
    end function subtract_Tensor4

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function Tensor2_matmul_Tensor2(tensorA,tensorB) result(this)
      class(Tensor2),intent(IN) :: tensorA,tensorB
      type(Tensor2) :: this

      if(TensorA%Initialized.and.TensorB%Initialized) then
         if (size(TensorA%data,2).eq.size(TensorB%data,1)) then
             this=new_Tensor(matmul(tensorA%data,tensorB%data))
         else
             call ThrowException('Tensor2_times_Tensor2','Tensors do not conform shape',NoErrorCode,CriticalError)
         endif
      else
         call ThrowException('Tensor2_times_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
      endif
      return
   end function Tensor2_matmul_Tensor2

   function Tensor1_matmul_Tensor2(tensorA,tensorB) result(this)
      class(Tensor1),intent(IN) :: tensorA
      class(Tensor2),intent(IN) :: tensorB
      type(Tensor1) :: this

      if(TensorA%Initialized.and.TensorB%Initialized) then
         this=new_Tensor(matmul(tensorA%data,tensorB%data))
      else
         call ThrowException('Tensor1_times_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
      endif
      return
   end function Tensor1_matmul_Tensor2

   function Tensor2_matmul_Tensor1(tensorA,tensorB) result(this)
      class(Tensor2),intent(IN) :: tensorA
      class(Tensor1),intent(IN) :: tensorB
      type(Tensor1) :: this

      if(TensorA%Initialized.and.TensorB%Initialized) then
         this=new_Tensor(matmul(tensorA%data,tensorB%data))
      else
         call ThrowException('Tensor2_times_Tensor1','Tensor not initialized',NoErrorCode,CriticalError)
      endif
      return
   end function Tensor2_matmul_Tensor1


   function Tensor1_dotProduct_Tensor1(tensorA,tensorB) result(this)
      class(Tensor1),intent(IN) :: tensorA,tensorB
      real(8) :: this

      if(TensorA%Initialized.and.TensorB%Initialized) then
         this=dot_product(tensorA%data,tensorB%data)
      else
         call ThrowException('Tensor1_times_Tensor1','Tensor not initialized',NoErrorCode,CriticalError)
      endif
      return
   end function Tensor1_dotProduct_Tensor1

    function Tensor3_matmul_Tensor3(tensorA,tensorB) result(this)
        class(Tensor3),intent(IN) :: tensorA,tensorB
        type(tensor4) :: this
        integer :: sumIndex,left1,left2,right1,right2
        integer :: dims_Of_A(3),dims_Of_B(3),new_Dims(4)
        complex(8) :: ctemp
        complex(8),allocatable :: anArray(:,:,:,:)

        if(TensorA%Initialized.and.TensorB%Initialized) then
            dims_Of_A=shape(tensorA%data)
            dims_Of_B=shape(tensorB%data)
            if(dims_Of_A(3).ne.dims_Of_B(1)) then
                call ThrowException('Tensor3_times_Tensor3','Tensor indexes have different size',NoErrorCode,CriticalError)
                return
            endif
            new_dims(1:2)=dims_Of_A(1:2)
            new_dims(3:4)=dims_Of_B(2:3)
            allocate(anArray(new_dims(1),new_dims(2),new_dims(3),new_dims(4)))
           !Structure inspired by BLAS3
            do right2=1,new_dims(4)
              do right1=1,new_dims(3)
                do sumIndex=1,dims_Of_A(3)
                  if(tensorB%data(sumIndex,right1,right2).ne.ZERO) then
                      ctemp=tensorB%data(sumIndex,right1,right2)
                      do left2=1,new_dims(2)
                        do left1=1,new_dims(1)
                            anArray(left1,left2,right1,right2)=anArray(left1,left2,right1,right2)+tensorA%data(left1,left2,sumIndex)*ctemp
                        enddo
                      enddo
                  endif
                enddo
              enddo
            enddo

            this=new_Tensor(anArray)
            deallocate(anArray)

          else
            call ThrowException('Tensor3_times_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
          endif

          return

       end function Tensor3_matmul_Tensor3

  function Tensor2_matmul_Tensor3(tensorA,tensorB) result(this)
        class(Tensor2),intent(IN) :: tensorA
        class(Tensor3),intent(IN) :: tensorB
        type(tensor3) :: this
        integer :: sumIndex,left1,right1,right2
        integer :: dims_Of_A(2),dims_Of_B(3),new_Dims(3)
        complex(8) :: ctemp
        complex(8),allocatable :: anArray(:,:,:)

        if(TensorA%Initialized.and.TensorB%Initialized) then
            dims_Of_A=shape(tensorA%data)
            dims_Of_B=shape(tensorB%data)
            if(dims_Of_A(2).ne.dims_Of_B(1)) then
                call ThrowException('Tensor2_times_Tensor3','Tensor indexes have different size',NoErrorCode,CriticalError)
                return
            endif
            new_dims(1)=dims_Of_A(1)
            new_dims(2:3)=dims_Of_B(2:3)
            allocate(anArray(new_dims(1),new_dims(2),new_dims(3)))
            anArray=ZERO
           !Structure inspired by BLAS3
            do right2=1,new_dims(3)
              do right1=1,new_dims(2)
                do sumIndex=1,dims_Of_A(2)
                  if(tensorB%data(sumIndex,right1,right2).ne.ZERO) then
                      ctemp=tensorB%data(sumIndex,right1,right2)
                      do left1=1,new_dims(1)
                          anArray(left1,right1,right2)=anArray(left1,right1,right2)+tensorA%data(left1,sumIndex)*ctemp
                      enddo
                  endif
                enddo
              enddo
            enddo

            this=new_Tensor(anArray)
            deallocate(anArray)

          else
            call ThrowException('Tensor2_times_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
          endif

          return

       end function Tensor2_matmul_Tensor3

   function Tensor3_matmul_Tensor2(tensorA,tensorB) result(this)
        class(Tensor3),intent(IN) :: tensorA
        class(Tensor2),intent(IN) :: tensorB
        type(tensor3) :: this
        integer :: sumIndex,left1,left2,right1
        integer :: dims_Of_A(3),dims_Of_B(2),new_Dims(3)
        complex(8) :: ctemp
        complex(8),allocatable :: anArray(:,:,:)

        if(TensorA%Initialized.and.TensorB%Initialized) then
            dims_Of_A=shape(tensorA%data)
            dims_Of_B=shape(tensorB%data)
            if(dims_Of_A(3).ne.dims_Of_B(1)) then
                call ThrowException('Tensor3_times_Tensor2','Tensor indexes have different size',NoErrorCode,CriticalError)
                return
            endif
            new_dims(1:2)=dims_Of_A(1:2)
            new_dims(3)=dims_Of_B(2)
            allocate(anArray(new_dims(1),new_dims(2),new_dims(3)))
            anArray=ZERO
           !Structure inspired by BLAS3
            do right1=1,new_dims(3)
              do sumIndex=1,dims_Of_A(3)
                if(tensorB%data(sumIndex,right1).ne.ZERO) then
                    ctemp=tensorB%data(sumIndex,right1)
                    do left2=1,new_dims(2)
                      do left1=1,new_dims(1)
                          anArray(left1,left2,right1)=anArray(left1,left2,right1)+tensorA%data(left1,left2,sumIndex)*ctemp
                      enddo
                    enddo
                endif
              enddo
            enddo

            this=new_Tensor(anArray)
            deallocate(anArray)

          else
            call ThrowException('Tensor3_times_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
          endif

          return

       end function Tensor3_matmul_Tensor2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function Tensor4_doubletimes_Tensor4 ( tensorA,tensorB) result(this)
        class(Tensor4),intent(IN) :: tensorA,tensorB
        type(tensor4) :: this
        integer :: sumIndex1,sumIndex2,left1,left2,right1,right2
        integer :: dims_Of_A(4),dims_Of_B(4),new_Dims(4)
        complex(8) :: ctemp
        complex(8),allocatable :: anArray(:,:,:,:)

        if(TensorA%Initialized.and.TensorB%Initialized) then
            dims_Of_A=shape(tensorA%data)
            dims_Of_B=shape(tensorB%data)
            if(dims_Of_A(3).ne.dims_Of_B(1).or.dims_Of_A(4).ne.dims_Of_B(2)) then
                call ThrowException('Tensor4_Doubletimes_Tensor4','Tensor indexes have different size',NoErrorCode,CriticalError)
                return
            endif
            new_dims(1:2)=dims_Of_A(1:2)
            new_dims(3:4)=dims_Of_B(3:4)
            allocate(anArray(new_dims(1),new_dims(2),new_dims(3),new_dims(4)))
            anArray=ZERO !Crucial to initialize to zero
            !Structure inspired by BLAS3
            do right2=1,new_dims(4)
              do right1=1,new_dims(3)
                do sumIndex2=1,dims_Of_A(4)
                  do sumIndex1=1,dims_Of_A(3)
                    if(tensorB%data(sumIndex1,sumIndex2,right1,right2).ne.ZERO) then
                      ctemp=tensorB%data(sumIndex1,sumIndex2,right1,right2)
                      do left2=1,new_dims(2)
                        do left1=1,new_dims(1)
                            anArray(left1,left2,right1,right2)=anArray(left1,left2,right1,right2)+tensorA%data(left1,left2,sumIndex1,sumIndex2)*ctemp
                        enddo
                      enddo
                    endif
                  enddo
                enddo
              enddo
            enddo

            this=new_Tensor(anArray)
            deallocate(anArray)

        else
            call ThrowException('Tensor4_doubletimes_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif

        return

    end function Tensor4_doubletimes_Tensor4

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function Tensor2_doubletimes_Tensor2 ( tensorA,tensorB) result(theResult)
        class(Tensor2),intent(IN) :: tensorA,tensorB
        complex(8) :: theResult

        if(TensorA%Initialized.and.TensorB%Initialized) then
            theResult=TensorTrace(tensorA*tensorB)
        else
            call ThrowException('Tensor2_doubletimes_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function Tensor2_doubletimes_Tensor2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function Tensor3_tripletimes_Tensor3 ( tensorA,tensorB) result(theResult)
        class(Tensor3),intent(IN) :: tensorA,tensorB
        complex(8) :: theResult
        integer :: a,b,c,dimsOfA(3),dimsOfB(3)

        if(TensorA%Initialized.and.TensorB%Initialized) then
            dimsofA=shape(tensorA%data)
            dimsofB=shape(tensorB%data)
            if (dimsofA.equalvector.dimsOfB) then
                theResult=ZERO
                do c=1,dimsOfA(3)
                  do b=1,dimsOfA(3)
                    do a=1,dimsOfA(3)
                      theResult=theResult+tensorA%data(a,b,c)*tensorB%data(a,b,c)
                    enddo
                  enddo
                enddo
            else
                call ThrowException('Tensor3_tripletimes_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
            endif
        else
            call ThrowException('Tensor3_tripletimes_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function Tensor3_tripletimes_Tensor3

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function Mirror_Compact_Left_With_Tensor3(this,aTensor,IndexToCompact) result(theResult)
       class(Tensor2),intent(IN) :: this
       integer :: IndexToCompact(1)
       class(Tensor3),intent(IN) :: aTensor
       type(Tensor2) :: theResult

       theResult=CompactLeft(this,aTensor,aTensor,IndexToCompact)

   end function Mirror_Compact_Left_With_Tensor3

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function Compact_Tensor3_From_Left_With_Tensor2(LeftTensor,upTensor,downTensor,IndexToCompact) result(theResult)
      class(Tensor3),intent(IN) :: upTensor,downTensor
      integer :: IndexToCompact(1)
      class(Tensor2),intent(IN) :: LeftTensor
      type(Tensor2) :: theResult
      complex(8),allocatable :: aMatrix(:,:)
      integer :: upDims(3), downDims(3), leftDims(2), n

      if(.not.(upTensor%Initialized .and. downTensor%Initialized .and. leftTensor%Initialized)) then
        call ThrowException('Compact_From_Left','Tensor not initialized',NoErrorCode,CriticalError)
        return
      endif

      upDims=shape(upTensor%data)
      downDims=shape(downTensor%data)
      leftDims=shape(leftTensor%data)

      if ( upDims(IndexToCompact(1)).ne.downDims(IndexToCompact(1)) ) then
        call ThrowException('Compact_From_Left','Index to Compact does not have equal dimensions',IndexToCompact(1),CriticalError)
        return
      endif

      If(IndexToCompact.equalvector.FIRST) then
          allocate ( aMatrix(downDims(3),upDims(3)) )
          aMatrix=ZERO
          do n=1,upDims(IndexToCompact(1))
            aMatrix=aMatrix+ matmul( Transpose(downTensor%data(n,:,:)) ,matmul( LeftTensor%data, upTensor%data(n,:,:) ) )
          enddo
      else if (IndexToCompact.equalvector.SECOND) then
          allocate ( aMatrix(downDims(3),upDims(3)) )
          aMatrix=ZERO
          do n=1,upDims(IndexToCompact(1))
            aMatrix=aMatrix+ matmul( Transpose(downTensor%data(:,n,:)) ,matmul( LeftTensor%data, upTensor%data(:,n,:) ) )
          enddo
      else if (IndexToCompact.equalvector.THIRD) then
          allocate ( aMatrix(downDims(2),upDims(2)) )
          aMatrix=ZERO
          do n=1,upDims(IndexToCompact(1))
            aMatrix=aMatrix+ matmul( Transpose(downTensor%data(:,:,n)) ,matmul( LeftTensor%data, upTensor%data(:,:,n) ) )
          enddo
      else
        call ThrowException('Compact_From_Left','Index to Compact is not 1,2, or 3',IndexToCompact(1),CriticalError)
        return
      endif

      theResult=new_Tensor(aMatrix)

   end function Compact_Tensor3_From_Left_With_Tensor2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function Mirror_Compact_Right_With_Tensor3(this,aTensor,IndexToCompact) result(theResult)
       class(Tensor2),intent(IN) :: this
       integer :: IndexToCompact(1)
       class(Tensor3),intent(IN) :: aTensor
       type(Tensor2) :: theResult

       theResult=CompactRight(this,aTensor,aTensor,IndexToCompact)

   end function Mirror_Compact_Right_With_Tensor3

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function Compact_Tensor3_From_Right_With_Tensor2(RightTensor,upTensor,downTensor,IndexToCompact) result(theResult)
      class(Tensor3),intent(IN) :: upTensor,downTensor
      integer :: IndexToCompact(1)
      class(Tensor2),intent(IN) :: RightTensor
      type(Tensor2) :: theResult
      complex(8),allocatable :: aMatrix(:,:)
      integer :: upDims(3), downDims(3), rightDims(2), n

      if(.not.(upTensor%Initialized .and. downTensor%Initialized .and. rightTensor%Initialized)) then
        call ThrowException('Compact_From_Right','Tensor not initialized',NoErrorCode,CriticalError)
        return
      endif

      upDims=shape(upTensor%data)
      downDims=shape(downTensor%data)
      rightDims=shape(rightTensor%data)

      if ( upDims(IndexToCompact(1)).ne.downDims(IndexToCompact(1)) ) then
        call ThrowException('Compact_From_Right','Index to Compact does not have equal dimensions',IndexToCompact(1),CriticalError)
        return
      endif

      If(IndexToCompact.equalvector.FIRST) then
          allocate ( aMatrix(upDims(2),downDims(2)) )
          aMatrix=ZERO
          do n=1,upDims(IndexToCompact(1))
            aMatrix=aMatrix+ matmul( matmul( upTensor%data(n,:,:), RightTensor%data ) , Transpose(downTensor%data(n,:,:)) )
          enddo
      else if (IndexToCompact.equalvector.SECOND) then
          allocate ( aMatrix(upDims(1),downDims(1)) )
          aMatrix=ZERO
          do n=1,upDims(IndexToCompact(1))
            aMatrix=aMatrix+ matmul( matmul( upTensor%data(:,n,:), RightTensor%data ) , Transpose(downTensor%data(:,n,:)) )
          enddo
      else if (IndexToCompact.equalvector.THIRD) then
          allocate ( aMatrix(upDims(1),downDims(1)) )
          aMatrix=ZERO
          do n=1,upDims(IndexToCompact(1))
            aMatrix=aMatrix+ matmul( matmul( upTensor%data(:,:,n), RightTensor%data ) , Transpose(downTensor%data(:,:,n)) )
          enddo
      else
        call ThrowException('Compact_From_Right','Index to Compact is not 1,2, or 3',IndexToCompact(1),CriticalError)
        return
      endif

      theResult=new_Tensor(aMatrix)

   end function Compact_Tensor3_From_Right_With_Tensor2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function Compact_From_Below_With_Tensor4(this,bound_index_of_3,aTensor4,bound_index_of_4,free_index_of_4) result(theResult)
     class(Tensor3),intent(IN) :: this
     class(Tensor4),intent(IN) :: aTensor4
     integer,intent(IN) :: bound_index_of_3(1),bound_index_of_4(1),free_index_of_4(1)
     type(Tensor3) :: theResult
     type(Tensor3) :: thisTransposed
     type(Tensor4) :: Tensor4Transposed
     integer :: permutation_Of_3(3),permutation_Of_4(4)
     integer :: dims_Of_4(4),dims_Of_3(3)
     integer :: leftDim,rightDim,downDim
     integer :: index,n
     integer :: freeDims3(2),freeDims4(2)

     if(.not.(this%Initialized)) then
        call ThrowException('Compact_From_Below_4','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

    dims_Of_3=shape(this%data)
    dims_Of_4=shape(aTensor4%data)

     if(dims_Of_3(bound_index_of_3(1)).ne.dims_Of_4(bound_index_of_4(1))) then
        call ThrowException('Compact_From_Below_4','Contracted index is not equal on tensors', &
           & dims_Of_3(bound_index_of_3(1))-dims_Of_4(bound_index_of_4(1)),CriticalError)
        return
     endif

    !Skip the bound indices to find the left and right dimensions
    index=1
    do n=1,3
      if (n.ne.bound_index_of_3(1)) then
         freeDims3(index)=n
         index=index+1
      endif
    enddo
    index=1
    do n=1,4
      if (n.ne.bound_index_of_4(1).and.n.ne.free_index_of_4(1)) then
         freeDims4(index)=n
         index=index+1
      endif
    enddo

    !Prepare permutation vectors for transposition
    permutation_Of_3(bound_index_of_3(1))=1
    permutation_Of_3(freeDims3(1))=2
    permutation_Of_3(freeDims3(2))=3
    permutation_Of_4(bound_index_of_4(1))=1
    permutation_Of_4(freeDims4(1))=2
    permutation_Of_4(freeDims4(2))=3
    permutation_Of_4(free_index_of_4(1))=4

    !Transpose
    thisTransposed=TensorTranspose(this,permutation_Of_3)
    Tensor4Transposed=TensorTranspose(aTensor4,permutation_Of_4)

    !Prepare permutation for the final result
    permutation_Of_3(1)=freeDims3(1)
    permutation_Of_3(2)=freeDims3(2)
    permutation_Of_3(3)=bound_index_of_3(1)

    theResult=TensorTranspose(MultAndCollapse(thisTransposed,Tensor4Transposed),permutation_Of_3)

   end function Compact_From_Below_With_Tensor4


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function MultAndCollapse_Tensor3_Tensor4(aTensor3,aTensor4) result(this)
    class(tensor3),intent(IN) :: aTensor3
    class(tensor4),intent(IN) :: aTensor4
    type(Tensor3) :: this
    integer :: dims_Of_4(4),dims_Of_3(3),new_dims(3)
    integer :: sumIndex,a,b,alpha,beta,freeIndex,leftIndex,rightIndex
    complex(8) :: ctemp
    complex(8),allocatable :: anArray(:,:,:)

    if(.not.(aTensor3%Initialized).or..not.(aTensor4%Initialized)) then
        call ThrowException('MultAndCollapse34','Tensor not initialized',NoErrorCode,CriticalError)
        return
    endif

    dims_Of_3=shape(aTensor3%data)
    dims_Of_4=shape(aTensor4%data)

    if(dims_Of_3(1).ne.dims_Of_4(1)) then
        call ThrowException('MultAndCollapse34','Contracted index is not equal on tensors', &
           & dims_Of_3(1)-dims_Of_4(1),CriticalError)
        return
    endif

    new_dims(1)=dims_Of_3(2)*dims_Of_4(2)
    new_dims(2)=dims_Of_3(3)*dims_Of_4(3)
    new_dims(3)=dims_Of_4(4)
    allocate(anArray(new_dims(1),new_dims(2),new_dims(3)))
    anArray=ZERO

    !Structure of loop inspired from BLAS level 3
    do freeIndex=1,dims_Of_4(4)
    do beta=1,dims_Of_4(3)
      do b=1,dims_Of_3(3)
        rightIndex=b+(beta-1)*dims_Of_3(3)
          do alpha=1,dims_Of_4(2)
            do a=1,dims_Of_3(2)
              leftIndex=a+(alpha-1)*dims_Of_3(2)
              ctemp=ZERO
              do sumIndex=1,dims_Of_3(1)
                 ctemp=ctemp+aTensor3%data(sumIndex,a,b)*aTensor4%data(sumIndex,alpha,beta,freeIndex)
              enddo
              anArray(leftIndex,rightIndex,freeIndex)=ctemp+anArray(leftIndex,rightIndex,freeIndex)
            enddo
          enddo
        enddo
      enddo
    enddo

    this=new_Tensor(anArray)

    end function MultAndCollapse_Tensor3_Tensor4
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function MultAndCollapse_Tensor3_Tensor3(TensorA,TensorB) result(this)
    class(tensor3),intent(IN) :: TensorA,TensorB
    type(Tensor3) :: this
    integer :: dims_Of_A(3),dims_Of_B(3),new_dims(3)
    integer :: sumIndex,a,b,beta,gama,newIndex
    complex(8) :: ctemp

    if(.not.(TensorA%Initialized).or..not.(TensorB%Initialized)) then
        call ThrowException('MultAndCollapse33','Tensor not initialized',NoErrorCode,CriticalError)
        return
    endif

    dims_Of_A=shape(TensorA%data)
    dims_Of_B=shape(TensorB%data)

    if(dims_Of_A(3).ne.dims_Of_B(1)) then
        call ThrowException('MultAndCollapse33','Contracted index is not equal on tensors', &
           & dims_Of_A(3)-dims_Of_B(1),CriticalError)
        return
    endif

    new_dims(1)=dims_Of_A(1)
    new_dims(2)=dims_Of_A(2)*dims_Of_B(2)
    new_dims(3)=dims_Of_B(3)
    this=new_Tensor(new_dims(1),new_dims(2),new_dims(3),ZERO)

    do gama=1,dims_Of_B(3)
      do beta=1,dims_Of_B(2)
        do sumIndex=1,dims_Of_B(1)
          ctemp=tensorB%data(sumIndex,beta,gama)
          do b=1,dims_Of_A(2)
            newIndex=b+(beta-1)*dims_Of_A(2)
            do a=1,dims_Of_A(1)
              this%data(a,newIndex,gama)=this%data(a,newIndex,gama)+tensorA%data(a,b,sumIndex)*ctemp
            enddo
          enddo
        enddo
      enddo
    enddo

    end function MultAndCollapse_Tensor3_Tensor3

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function MultAndCollapse_Tensor4_Tensor4 (TensorA,TensorB) result(this)
    class(tensor4),intent(IN) :: TensorA,TensorB
    type(Tensor4) :: this
    integer :: dims_Of_A(4),dims_Of_B(4),new_dims(4)
    integer :: sumIndex,a,b,c,d,alpha,beta,gama,delta,newIndex1,newIndex2
    complex(8) :: ctemp

    if(.not.(TensorA%Initialized).or..not.(TensorB%Initialized)) then
        call ThrowException('MultAndCollapse44','Tensor not initialized',NoErrorCode,CriticalError)
        return
    endif

    dims_Of_A=shape(TensorA%data)
    dims_Of_B=shape(TensorB%data)

    if(dims_Of_A(2).ne.dims_Of_B(1)) then
        call ThrowException('MultAndCollapse44','Contracted index is not equal on tensors', &
           & dims_Of_A(2)-dims_Of_B(1),CriticalError)
        return
    endif

    new_dims(1)=dims_Of_A(1)
    new_dims(2)=dims_Of_B(2)
    new_dims(3)=dims_Of_A(3)*dims_Of_B(3)
    new_dims(4)=dims_Of_A(4)*dims_Of_B(4)
    this=new_Tensor(new_dims(1),new_dims(2),new_dims(3),new_dims(4),ZERO)

    do delta=1,dims_Of_B(4)
     do gama=1,dims_Of_B(3)
      do beta=1,dims_Of_B(2)
          do d=1,dims_Of_A(4)
           newIndex2=d+(delta-1)*dims_Of_A(4)
           do c=1,dims_Of_A(3)
             newIndex1=c+(gama-1)*dims_Of_A(3)
             do sumIndex=1,dims_Of_A(2)
               ctemp=tensorB%data(sumIndex,beta,gama,delta)
               do a=1,dims_Of_A(1)
                this%data(a,beta,newIndex1,newIndex2)=this%data(a,beta,newIndex1,newIndex2)+tensorA%data(a,sumIndex,c,d)*ctemp
               enddo
             enddo
           enddo
          enddo
       enddo
     enddo
    enddo

    end function MultAndCollapse_Tensor4_Tensor4

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!!##################################################################
!!##################################################################
!!##################################################################

   real(8) function Difference_btw_Tensors(tensorA, tensorB) result(diff)
     class(Tensor),intent(IN) :: tensorA,tensorB

     diff=0.0d0

	if(.not.(tensorA%Initialized .and. tensorB%Initialized)) then
		call ThrowException('Difference_btw_Tensors','Tensor not initialized',NoErrorCode,CriticalError)
        return
	endif
	if( .not. same_type_as(tensorA,tensorB) ) then
		call ThrowException('Difference_btw_Tensors','Tensor not of same type',NoErrorCode,CriticalError)
        return
	endif
	if (tensorA.equaldims.tensorB) then
	!The following ugly structure is the best I have to cast the tensors into
	!their corresponding type
	select type (Typed_A => tensorA)
	    class is (Tensor1)
			select type (Typed_B => tensorB)
				class is (Tensor1)
					diff=sum(Typed_A%data-Typed_B%data)
				class default
			    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    			    return
			end select
	 	class is (Tensor2)
			select type (Typed_B => tensorB)
				class is (Tensor2)
					diff=sum(Typed_A%data-Typed_B%data)
				class default
			    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    			    return
			end select
	 	class is (Tensor3)
			select type (Typed_B => tensorB)
				class is (Tensor3)
					diff=sum(Typed_A%data-Typed_B%data)
				class default
			    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    			    return
			end select
        class is (Tensor4)
            select type (Typed_B => tensorB)
                class is (Tensor4)
                    diff=sum(Typed_A%data-Typed_B%data)
                class default
                    call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
                    return
            end select
	 	class is (Tensor)
	    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    	    return
	 end select
    else !Not equal shape
        call ThrowException('Difference_btw_Tensors','Tensor not of same shape',NoErrorCode,CriticalError)
    endif

   end function Difference_btw_Tensors


   real function Difference_btw_Tensors_WithAbsoluteValue(tensorA, tensorB) result(diff)
     class(Tensor),intent(IN) :: tensorA,tensorB

     diff=0.0d0

	if(.not.(tensorA%Initialized .and. tensorB%Initialized)) then
		call ThrowException('Difference_btw_Tensors','Tensor not initialized',NoErrorCode,CriticalError)
        return
	endif
	if( .not. same_type_as(tensorA,tensorB) ) then
		call ThrowException('Difference_btw_Tensors','Tensor not of same type',NoErrorCode,CriticalError)
        return
	endif

	if (tensorA.equaldims.tensorB) then
	!The following ugly structure is the best I have to cast the tensors into
	!their corresponding type
	select type (Typed_A => tensorA)
	    class is (Tensor1)
			select type (Typed_B => tensorB)
				class is (Tensor1)
					diff=sum(abs(Typed_A%data-Typed_B%data))
				class default
			    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    			    return
			end select
	 	class is (Tensor2)
			select type (Typed_B => tensorB)
				class is (Tensor2)
					diff=sum(abs(Typed_A%data-Typed_B%data))
				class default
			    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    			    return
			end select
	 	class is (Tensor3)
			select type (Typed_B => tensorB)
				class is (Tensor3)
					diff=sum(abs(Typed_A%data-Typed_B%data))
				class default
			    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    			    return
			end select
        class is (Tensor4)
            select type (Typed_B => tensorB)
                class is (Tensor4)
                    diff=sum(abs(Typed_A%data-Typed_B%data))
                class default
                    call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
                    return
            end select
	 	class is (Tensor)
	    	call ThrowException('Difference_btw_Tensors','Unknown error',NoErrorCode,CriticalError)
    	    return
	 end select
    else !Not equal shape
        call ThrowException('Difference_btw_Tensors','Tensor not of same shape',NoErrorCode,CriticalError)
    endif

   end function Difference_btw_Tensors_WithAbsoluteValue

!##################################################################

   logical function Tensors_are_of_equal_Shape(tensorA,tensorB) result(equals)
     class(Tensor),intent(IN) :: tensorA,tensorB

     if(.not.(tensorA%Initialized .and. tensorB%Initialized)) then
        call ThrowException('MPSTensors_are_of_equal_Shape','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif
	if( .not. same_type_as(tensorA,tensorB) ) then
		call ThrowException('Tensors_are_of_equal_Shape','Tensor not of same type',NoErrorCode,CriticalError)
        return
	endif
	select type (Typed_A => tensorA)
      class is (Tensor1)
         select type (Typed_B => tensorB)
             class is (Tensor1)
                equals=sum(abs(shape(typed_A%data)-shape(typed_B%data))).eq.ZERO
             class default
                call ThrowException('Tensors_are_equal','Unknown error',NoErrorCode,CriticalError)
                return
         end select
      class is (Tensor2)
         select type (Typed_B => tensorB)
             class is (Tensor2)
                equals=sum(abs(shape(typed_A%data)-shape(typed_B%data))).eq.ZERO
             class default
                call ThrowException('Tensors_are_equal','Unknown error',NoErrorCode,CriticalError)
                return
         end select
      class is (Tensor3)
         select type (Typed_B => tensorB)
             class is (Tensor3)
                equals=sum(abs(shape(typed_A%data)-shape(typed_B%data))).eq.ZERO
             class default
                call ThrowException('Tensors_are_equal','Unknown error',NoErrorCode,CriticalError)
                return
         end select
      class is (Tensor4)
         select type (Typed_B => tensorB)
             class is (Tensor4)
                equals=sum(abs(shape(typed_A%data)-shape(typed_B%data))).eq.ZERO
             class default
                call ThrowException('Tensors_are_equal','Unknown error',NoErrorCode,CriticalError)
                return
         end select
      end select

   end function Tensors_are_of_equal_Shape

  logical function Tensors_are_of_equal_Type(tensorA,tensorB) result(equals)
     class(Tensor),intent(IN) :: tensorA,tensorB

     if(.not.(tensorA%Initialized .and. tensorB%Initialized)) then
        call ThrowException('MPSTensors_are_of_equal_Shape','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif
	if( .not. same_type_as(tensorA,tensorB) ) then
		call ThrowException('Tensors_are_of_equal_Shape','Tensor not of same type',NoErrorCode,CriticalError)
        return
	endif
	equals=same_type_as(tensorA,tensorB)
	return

  end function Tensors_are_of_equal_Type

!##################################################################

!POSIBLE CALLING FORMS:
!aTensor%JoinIndices(FirstAndSecond,Third)
!aTensor%JoinIndices(Second,FirstAndThird)
!aTensor%JoinIndices(ThirdAndSecond,First)
!aTensor%JoinIndices(ThirdAndFirst,Second)
!aTensor%JoinIndices(SecondAndFirst,Third)
!
!Tensor4 will be  *  combinations
!1 2 34         2 * 6
!1 3 24         2 * 6
!1 4 23         2 * 6
!2 3 14         2 * 6
!2 4 13         2 * 6
!3 4 12         2 * 6
!12 34         2 * 4
!13 24         2 * 4
!14 23         2 * 4
!1 234         2 * 6
!2 134         2 * 6
!3 124         2 * 6
!1234         1
! For a total of 44*3 combinations, well beyond a reasonable If-else chain


function JoinIndicesOfTensor3(this,firstindex,secondindex) result (aTensor)
	integer,intent(IN) :: firstindex(:),secondindex(:)
	class(tensor3),intent(IN) :: this
	type(tensor2) :: aTensor
	integer :: i,dim(3)
	logical :: isNormalOrder
	complex(8),allocatable :: anArray(:,:)

     if(.not.(this%Initialized)) then
        call ThrowException('JoinIndicesOfTensor3','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

	dim=shape(this%data)

	if ((firstindex .equalvector. FIRST ).or. (secondindex .equalvector. FIRST)) then
		allocate(anArray(dim(1),dim(2)*dim(3)))
		if ((secondindex .equalvector. SECONDANDTHIRD ).or. (firstindex .equalvector. SECONDANDTHIRD) ) then
			do i=1,dim(1)
				anArray(i,:)=reshape( this%data(i,:,:), [dim(2)*dim(3)])
			enddo
		else if ((secondindex .equalvector. THIRDANDSECOND) .or. (firstindex .equalvector. THIRDANDSECOND)) then
			do i=1,dim(1)
				anArray(i,:)=reshape( transpose(this%data(i,:,:)), [dim(2)*dim(3)])
			enddo
		else
	        call ThrowException('JoinIndicesOfTensor3','Incorrect joint index',NoErrorCode,CriticalError)
    	    return
		endif
        if (firstindex .equalvector. FIRST) then
            isNormalOrder=.true.
        else
            isNormalOrder=.false.
        endif
	else if((firstindex .equalvector. SECOND) .or. (secondindex .equalvector. SECOND)) then
		allocate(anArray(dim(2),dim(1)*dim(3)))
		if ((secondindex .equalvector. FIRSTANDTHIRD) .or. (firstindex .equalvector. FIRSTANDTHIRD) ) then
			do i=1,dim(2)
				anArray(i,:)=reshape( this%data(:,i,:), [dim(1)*dim(3)])
			enddo
		else if ((secondindex .equalvector. THIRDANDFIRST) .or. (firstindex .equalvector. THIRDANDFIRST)) then
			do i=1,dim(2)
				anArray(i,:)=reshape( transpose(this%data(:,i,:)), [dim(1)*dim(3)])
			enddo
		else
	        call ThrowException('JoinIndicesOfTensor3','Incorrect joint index',NoErrorCode,CriticalError)
    	    return
		endif
        if (firstindex .equalvector. SECOND) then
            isNormalOrder=.true.
        else
            isNormalOrder=.false.
        endif
	else if((firstindex .equalvector. THIRD) .or. (secondindex .equalvector. THIRD)) then
		allocate(anArray(dim(3),dim(1)*dim(2)))
		if ((secondindex .equalvector. FIRSTANDSECOND) .or.( firstindex .equalvector. FIRSTANDSECOND) ) then
			do i=1,dim(3)
				anArray(i,:)=reshape( this%data(:,:,i), [dim(1)*dim(2)])
			enddo
		else if ((secondindex .equalvector. SECONDANDFIRST) .or. (firstindex .equalvector. SECONDANDFIRST)) then
			do i=1,dim(3)
				anArray(i,:)=reshape( transpose(this%data(:,:,i)), [dim(1)*dim(2)])
			enddo
		else
	        call ThrowException('JoinIndicesOfTensor3','Incorrect joint index',NoErrorCode,CriticalError)
    	    return
		endif
        if (firstindex .equalvector. THIRD) then
            isNormalOrder=.true.
        else
            isNormalOrder=.false.
        endif

	else
        call ThrowException('JoinIndicesOfTensor3','Incorrect single index',NoErrorCode,CriticalError)
    	return
	endif

	!Finally create the tensor checking for the reverse order or normal order
	if (isNormalOrder) then
		aTensor=new_Tensor(anArray)
	else
		aTensor=new_Tensor(Transpose(anArray))
	endif

end function JoinIndicesOfTensor3

!##################################################################

function JoinIndicesOfTensor4(this,firstindex,secondindex) result (aTensor)
    integer,intent(IN) :: firstindex(:),secondindex(:)
    class(tensor4),intent(IN) :: this
    type(tensor2) :: aTensor
    integer :: dim(4),indexVector(4)
    integer :: x1,x2,x3,x4
    integer :: leftDim,rightDim,leftIndex,rightIndex
    complex(8),allocatable :: anArray(:,:)

     if(.not.(this%Initialized)) then
        call ThrowException('JoinIndicesOfTensor4','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

    !TODO::  ERROR CHECKING of INPUT DIMENSIONS
    dim=shape(this%data)

    leftDim=dim(firstindex(1))*dim(firstindex(2))
    RightDim=dim(secondindex(1))*dim(secondindex(2))

    allocate(anArray(leftDim,RightDim))

    do x4=1,dim(4)
      indexVector(4)=x4
      do x3=1,dim(3)
        indexVector(3)=x3
        do x2=1,dim(2)
          indexVector(2)=x2
          do x1=1,dim(1)
            indexVector(1)=x1
            LeftIndex=indexVector(firstindex(1))+(indexVector(firstindex(2))-1)*dim(firstindex(1))
            RightIndex=indexVector(secondindex(1))+(indexVector(secondindex(2))-1)*dim(secondindex(1))
            anArray(leftIndex,rightIndex)=this%data(x1,x2,x3,x4)
          enddo
        enddo
      enddo
    enddo

    aTensor=new_Tensor(anArray)

end function JoinIndicesOfTensor4

!##################################################################

!##################################################################


!##################################################################

function Take_Slice_Of_Tensor4(this,range1,range2,range3,range4) result (aTensor)
    integer,intent(IN) :: range1(2),range2(2),range3(2),range4(2)
    class(tensor4),intent(IN) :: this
    type(tensor4) :: aTensor
    integer :: dims(4)

     if(this%Initialized) then
        dims=shape(this%data)
        if(range1(1).ge.1.and.range2(1).ge.1.and.range3(1).ge.1.and.range4(1).ge.1.and. &
          & range1(2).le.dims(1).and.range2(2).le.dims(2).and.range3(2).le.dims(3).and.range4(2).le.dims(4) ) then
            aTensor=new_Tensor(this%data(range1(1):range1(2),range2(1):range2(2),range3(1):range3(2),range4(1):range4(2)))
	  	else
	  	    call ThrowException('Take_Slice_Of_Tensor4','Index is inappropriate',noErrorCode,CriticalError)
	  	endif
     else
        call ThrowException('Take_Slice_Of_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Take_Slice_Of_Tensor4

!##################################################################

function Take_Slice_Of_Tensor3(this,range1,range2,range3) result (aTensor)
    integer,intent(IN) :: range1(2),range2(2),range3(2)
    class(tensor3),intent(IN) :: this
    type(tensor3) :: aTensor
    integer :: dims(3)

     if(this%Initialized) then
         dims=shape(this%data)
         if(range1(1).ge.1.and.range2(1).ge.1.and.range3(1).ge.1.and. &
          & range1(2).le.dims(1).and.range2(2).le.dims(2).and.range3(2).le.dims(3)) then
                aTensor=new_Tensor(this%data(range1(1):range1(2),range2(1):range2(2),range3(1):range3(2)))
        else
            call ThrowException('Take_Slice_Of_Tensor3','Index is inappropriate',NoErrorCode,CriticalError)
        endif
     else
        call ThrowException('Take_Slice_Of_Tensor3','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Take_Slice_Of_Tensor3

!##################################################################

function Take_Slice_Of_Tensor2(this,range1,range2) result (aTensor)
    integer,intent(IN) :: range1(2),range2(2)
    class(tensor2),intent(IN) :: this
    type(tensor2) :: aTensor
    integer :: dims(2)

     if(this%Initialized) then
         dims=shape(this%data)
         if(range1(1).ge.1.and.range2(1).ge.1.and. &
          & range1(2).le.dims(1).and.range2(2).le.dims(2)) then
                aTensor=new_Tensor(this%data(range1(1):range1(2),range2(1):range2(2)))
        else
            call ThrowException('Take_Slice_Of_Tensor2','Index is inappropriate',NoErrorCode,CriticalError)
        endif
     else
        call ThrowException('Take_Slice_Of_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Take_Slice_Of_Tensor2

!##################################################################

function Take_Slice_Of_Tensor1(this,range1) result (aTensor)
    integer,intent(IN) :: range1(2)
    class(tensor1),intent(IN) :: this
    type(tensor1) :: aTensor
    integer :: dims(1)

     if(this%Initialized) then
         dims=shape(this%data)
         if(range1(1).ge.1.and.range1(2).le.dims(1)) then
                aTensor=new_Tensor(this%data(range1(1):range1(2)))
        else
            call ThrowException('Take_Slice_Of_Tensor1','Index is inappropriate',NoErrorCode,CriticalError)
        endif
     else
        call ThrowException('Take_Slice_Of_Tensor1','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Take_Slice_Of_Tensor1

!##################################################################
!##################################################################

function SplitIndexOfTensor2(this,WhichIndex,Partition) result (aTensor)
    integer,intent(IN) :: WhichIndex(:),Partition
    class(tensor2),intent(IN) :: this
    type(tensor3) :: aTensor
    integer :: i,NewDims(3),OrigDims(2),error

     if(.not.(this%Initialized)) then
        call ThrowException('SplitIndexOfTensor2','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

    OrigDims=shape(this%data)

    !Ugly notation with (1) follows because select case needs numbers

    !First check for error in multiple
    select case (WhichIndex(1))
      case (first(1))
        error=mod( OrigDims(1),Partition )
      case (second(1))
        error=mod( OrigDims(2),Partition )
      case default
        call ThrowException('SplitIndexOfTensor2','Index is inappropriate',WhichIndex(1),CriticalError)
        return
    end select
    if(error.ne.0.or.Partition.le.0.or.Partition.gt.maxval(OrigDims)) then
        call ThrowException('SplitIndexOfTensor2','Requested partition seems incorrect',Partition,CriticalError)
        return
     endif

    !Now for real, I need to repeat this code or I have to repeat the error checking
    select case (WhichIndex(1))
      case (first(1))
        NewDims=[ Partition, OrigDims(1)/Partition, OrigDims(2) ]
      case (second(1))
        NewDims=[ OrigDims(1), Partition, OrigDims(2)/Partition ]
    end select

    aTensor=new_Tensor(reshape(this%data,NewDims))

end function SplitIndexOfTensor2

!##################################################################

function SplitIndexOfTensor3(this,WhichIndex,Partition) result (aTensor)
    integer,intent(IN) :: WhichIndex(:),Partition
    class(tensor3),intent(IN) :: this
    type(tensor4) :: aTensor
    integer :: i,NewDims(4),OrigDims(3),error

     if(.not.(this%Initialized)) then
        call ThrowException('SplitIndexOfTensor3','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

    OrigDims=shape(this%data)

    !Ugly notation with (1) follows because select case needs numbers

    !First check for error in multiple
    select case (WhichIndex(1))
      case (first(1))
        error=mod( OrigDims(1),Partition )
      case (second(1))
        error=mod( OrigDims(2),Partition )
      case (third(1))
        error=mod( OrigDims(3),Partition )
      case default
        call ThrowException('SplitIndexOfTensor3','Index is inappropriate',WhichIndex(1),CriticalError)
        return
    end select
    if(error.ne.0.or.Partition.le.0.or.Partition.gt.maxval(OrigDims)) then
        call ThrowException('SplitIndexOfTensor3','Requested partition seems incorrect',Partition,CriticalError)
        return
     endif

    !Now for real, I need to repeat this code or I have to repeat the error checking
    select case (WhichIndex(1))
      case (first(1))
        NewDims=[ Partition, OrigDims(1)/Partition, OrigDims(2), OrigDims(3) ]
      case (second(1))
        NewDims=[ OrigDims(1), Partition, OrigDims(2)/Partition, OrigDims(3)]
      case (third(1))
        NewDims=[ OrigDims(1), OrigDims(2), Partition, OrigDims(3)/Partition]
    end select

    aTensor=new_Tensor(reshape(this%data,NewDims))

end function SplitIndexOfTensor3

!##################################################################

function Pad_Tensor2(this,newDims) result(reshapedTensor)
    class(Tensor2),intent(IN) :: this
    integer,intent(IN) :: newDims(2)
    type(Tensor2) :: reshapedTensor
    integer :: oldDims(2),minDims(2)

    if(this%Initialized) then
        reshapedTensor=new_Tensor(newDims(1),newDims(2),ZERO)
        oldDims=this%GetDimensions()
        minDims(1)=min(oldDims(1),newDims(1))
        minDims(2)=min(oldDims(2),newDims(2))
        reshapedTensor%data(1:minDims(1),1:minDims(2))=this%data(1:minDims(1),1:minDims(2))
    else
        call ThrowException('Pad_Tensor2','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return

end function Pad_Tensor2
!##################################################################

function ConjugateTensor1(this) result(thisdagger)
   class(Tensor1),intent(IN) :: this
   type(Tensor1) :: thisdagger

   if(this%Initialized) then
       thisdagger=new_Tensor1_fromData(dconjg(this%data))
   else
       call ThrowException('Conjugate1','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTensor1

function ConjugateTensor2(this) result(thisdagger)
   class(Tensor2),intent(IN) :: this
   type(Tensor2) :: thisdagger

   if(this%Initialized) then
       thisdagger=new_Tensor2_fromData(dconjg(this%data))
   else
       call ThrowException('Conjugate2','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTensor2

function ConjugateTensor3(this) result(thisdagger)
   class(Tensor3),intent(IN) :: this
   type(Tensor3) :: thisdagger

   if(this%Initialized) then
       thisdagger=new_Tensor3_fromData(dconjg(this%data))
   else
       call ThrowException('Conjugate3','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTensor3

function ConjugateTensor4(this) result(thisdagger)
   class(Tensor4),intent(IN) :: this
   type(Tensor4) :: thisdagger

   if(this%Initialized) then
       thisdagger=new_Tensor4_fromData(dconjg(this%data))
   else
       call ThrowException('Conjugate4','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTensor4


!##################################################################

function ConjugateTranspose2(this) result(thisdagger)
   class(Tensor2),intent(IN) :: this
   type(Tensor2) :: thisdagger

   if(this%Initialized) then
       thisdagger=new_Tensor2_fromData(dconjg(transpose(this%data)))
   else
       call ThrowException('ConjugateTranspose','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTranspose2

function ConjugateTranspose3(this,permutation) result(thisdagger)
   class(Tensor3),intent(IN) :: this
   type(Tensor3) :: thisdagger
   integer,intent(IN) :: permutation(3)

   if(this%Initialized) then
        thisdagger=TensorTranspose(Conjugate(this),permutation)
   else
       call ThrowException('ConjugateTranspose3','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTranspose3

function ConjugateTranspose4(this,permutation) result(thisdagger)
   class(Tensor4),intent(IN) :: this
   type(Tensor4) :: thisdagger
   integer,intent(IN) :: permutation(4)

   if(this%Initialized) then
        thisdagger=TensorTranspose(Conjugate(this),permutation)
   else
       call ThrowException('ConjugateTranspose4','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function ConjugateTranspose4



function TensorTranspose2(this) result(thisdagger)
   class(Tensor2),intent(IN) :: this
   type(Tensor2) :: thisdagger

   if(this%Initialized) then
       thisdagger=new_Tensor2_fromData(transpose(this%data))
   else
       call ThrowException('TensorTranspose2','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function TensorTranspose2

! HELP :: each permutation index says to where the dimension is moved
! example :: (3,1,2) means that the first dimension is now the third, the second the first,
!            and the last the second.
function TensorTranspose3(this,permutation) result(thisdagger)
   class(Tensor3),intent(IN) :: this
   type(Tensor3) :: thisdagger
   integer,intent(IN) :: permutation(3)
   integer :: dims(3),newdims(3),n

   if(this%Initialized) then
        if(6.eq.permutation(1)*permutation(2)*permutation(3).and.6.eq.(permutation(1)+permutation(2)+permutation(3))) then
            dims=shape(this%data)
            do n=1,3
                newdims(permutation(n))=dims(n)
            enddo
            thisdagger=new_Tensor3_fromData(reshape(this%data,newdims,ORDER=permutation))
        else
            call ThrowException('TensorTranspose3','Order must be permutation of 1,2,3',NoErrorCode,CriticalError)
        endif
   else
       call ThrowException('TensorTranspose3','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function TensorTranspose3

function TensorTranspose4(this,permutation) result(thisTransposed)
   class(Tensor4),intent(IN) :: this
   type(Tensor4) :: thisTransposed
   integer,intent(IN) :: permutation(4)
   integer :: dims(4),newdims(4),n

   if(this%Initialized) then
        if(24.eq.permutation(1)*permutation(2)*permutation(3)*permutation(4).and.10.eq.sum(permutation)) then
            dims=shape(this%data)
            do n=1,size(dims)
                newdims(permutation(n))=dims(n)
            enddo
            thisTransposed=new_Tensor4_fromData(reshape(this%data,newdims,ORDER=permutation))
        else
            call ThrowException('TensorTranspose4','Order must be permutation of 1,2,3,4',NoErrorCode,CriticalError)
        endif
   else
       call ThrowException('TensorTranspose4','Tensor not initialized',NoErrorCode,CriticalError)
   endif
   return
end function TensorTranspose4


!##################################################################

complex(8) function Tensor2Trace(this) result(theTrace)
    class(Tensor2),intent(IN) :: this
    integer :: sumIndex
    integer :: dims(2)

    if(this%Initialized) then
        dims=shape(this%data)
        theTrace=ZERO
        do sumIndex=1,min(dims(1),dims(2))
            theTrace=theTrace+this%data(sumIndex,sumIndex)
        enddo
     else
        call ThrowException('Tensor2Trace','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Tensor2Trace

!##################################################################

 function Tensor3Trace(this,survivingDimension) result(TracedTensor)
    class(Tensor3),intent(IN) :: this
    type(Tensor1) :: TracedTensor
    integer, intent(IN) :: survivingDimension(1)
    integer :: sumIndex,survivingIndex
    integer :: dims(3)

    if(this%Initialized) then
        dims=shape(this%data)
        select case (survivingDimension(1))
          case (FIRST(1))
             TracedTensor=new_Tensor(dims(1),ZERO)
             do sumIndex=1,min(dims(2),dims(3))
              do survivingIndex=1,dims(1)
                TracedTensor%data(survivingIndex)=TracedTensor%data(survivingIndex)+this%data(survivingIndex,sumIndex,sumIndex)
              enddo
             enddo
          case (SECOND(1))
             TracedTensor=new_Tensor(dims(2),ZERO)
             do sumIndex=1,min(dims(1),dims(3))
              do survivingIndex=1,dims(2)
                TracedTensor%data(survivingIndex)=TracedTensor%data(survivingIndex)+this%data(sumIndex,survivingIndex,sumIndex)
              enddo
             enddo
          case (THIRD(1))
             TracedTensor=new_Tensor(dims(3),ZERO)
             do survivingIndex=1,dims(3)
              do sumIndex=1,min(dims(2),dims(1))
                TracedTensor%data(survivingIndex)=TracedTensor%data(survivingIndex)+this%data(sumIndex,sumIndex,survivingIndex)
              enddo
             enddo
          case default
           call ThrowException('Tensor3Trace','Index is inappropriate',survivingDimension(1),CriticalError)
           return
         end select
      else
        call ThrowException('Tensor3Trace','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Tensor3Trace

!##################################################################

 function Tensor4Trace(this,survivingDimensions) result(TracedTensor)
    class(Tensor4),intent(IN) :: this
    type(Tensor2) :: TracedTensor
    integer, intent(IN) :: survivingDimensions(2)
    integer :: sumIndex,survivingIndex1,survivingIndex2
    integer :: dims(4)

    if(this%Initialized) then
        dims=shape(this%data)
        select case (survivingDimensions(1))
          case (FIRST(1))
             select case (survivingDimensions(2))
               case(SECOND(1))
	             TracedTensor=new_Tensor(dims(1),dims(2),ZERO)
	             do sumIndex=1,min(dims(3),dims(4))
	              do survivingIndex2=1,dims(2)
	              do survivingIndex1=1,dims(1)
	                TracedTensor%data(survivingIndex1,survivingIndex2)=TracedTensor%data(survivingIndex1,survivingIndex2)+ &
	                  & this%data(survivingIndex1,survivingIndex2,sumIndex,sumIndex)
	              enddo
	              enddo
	             enddo
               case(THIRD(1))
                 TracedTensor=new_Tensor(dims(1),dims(3),ZERO)
                 do sumIndex=1,min(dims(2),dims(4))
                  do survivingIndex2=1,dims(3)
                  do survivingIndex1=1,dims(1)
                    TracedTensor%data(survivingIndex1,survivingIndex2)=TracedTensor%data(survivingIndex1,survivingIndex2)+ &
                      & this%data(survivingIndex1,sumIndex,survivingIndex2,sumIndex)
                  enddo
                  enddo
                 enddo
               case(FOURTH(1))
                 TracedTensor=new_Tensor(dims(1),dims(4),ZERO)
                 do sumIndex=1,min(dims(2),dims(3))
                  do survivingIndex2=1,dims(4)
                  do survivingIndex1=1,dims(1)
                    TracedTensor%data(survivingIndex1,survivingIndex2)=TracedTensor%data(survivingIndex1,survivingIndex2)+ &
                      & this%data(survivingIndex1,sumIndex,sumIndex,survivingIndex2)
                  enddo
                  enddo
                 enddo
	         case default
	           call ThrowException('Tensor4Trace','2nd Index is inappropriate',survivingDimensions(2),CriticalError)
	           return
             end select
          case (SECOND(1))
             select case (survivingDimensions(2))
               case(THIRD(1))
                 TracedTensor=new_Tensor(dims(2),dims(3),ZERO)
                 do sumIndex=1,min(dims(1),dims(4))
                  do survivingIndex2=1,dims(3)
                  do survivingIndex1=1,dims(2)
                    TracedTensor%data(survivingIndex1,survivingIndex2)=TracedTensor%data(survivingIndex1,survivingIndex2)+ &
                      & this%data(sumIndex,survivingIndex1,survivingIndex2,sumIndex)
                  enddo
                  enddo
                 enddo
               case(FOURTH(1))
                 TracedTensor=new_Tensor(dims(2),dims(4),ZERO)
                 do sumIndex=1,min(dims(1),dims(3))
                  do survivingIndex2=1,dims(4)
                  do survivingIndex1=1,dims(2)
                    TracedTensor%data(survivingIndex1,survivingIndex2)=TracedTensor%data(survivingIndex1,survivingIndex2)+ &
                      & this%data(sumIndex,survivingIndex1,sumIndex,survivingIndex2)
                  enddo
                  enddo
                 enddo
             case default
                call ThrowException('Tensor4Trace','2nd Index is inappropriate',survivingDimensions(2),CriticalError)
                return
             end select
          case (THIRD(1))
             select case (survivingDimensions(2))
               case(FOURTH(1))
                 TracedTensor=new_Tensor(dims(3),dims(4),ZERO)
                 do sumIndex=1,min(dims(1),dims(2))
                  do survivingIndex2=1,dims(4)
                  do survivingIndex1=1,dims(3)
                    TracedTensor%data(survivingIndex1,survivingIndex2)=TracedTensor%data(survivingIndex1,survivingIndex2)+ &
                      & this%data(sumIndex,sumIndex,survivingIndex1,survivingIndex2)
                  enddo
                  enddo
                 enddo
             case default
               call ThrowException('Tensor4Trace','2nd Index is inappropriate',survivingDimensions(2),CriticalError)
               return
             end select
          case default
           call ThrowException('Tensor4Trace','1st Index is inappropriate',survivingDimensions(1),CriticalError)
           return
         end select
      else
        call ThrowException('Tensor4Trace','Tensor not initialized',NoErrorCode,CriticalError)
     endif
end function Tensor4Trace

!##################################################################
!##################################################################
!##################################################################

    function UnfoldTensor4(this,survivingIndex) result(aMatrix)
        class(Tensor4),intent(IN) :: this
        integer,intent(IN) :: survivingIndex
        type(Tensor2) :: aMatrix
        integer :: dims(4),tempdims(4),reordering(4),newdims(2)

        if(this%Initialized) then
            reordering=[1,2,3,4]
            reordering(1)=survivingIndex
            reordering(survivingIndex)=1
            dims=shape(this%data)
            tempdims=dims
            tempdims(1)=dims(survivingIndex)
            tempdims(survivingIndex)=dims(1)
            newdims(1)=tempdims(1)
            newdims(2)=product(tempdims(2:4))
            aMatrix=new_Tensor( newdims(1), newdims(2) , ZERO )
            aMatrix%data=reshape( reshape(this%data, tempdims, ORDER=reordering), newdims)
        else
            call ThrowException('UnfoldTensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function UnfoldTensor4
!##################################################################

    function UnfoldTensor3(this,survivingIndex) result(aMatrix)
        class(Tensor3),intent(IN) :: this
        integer,intent(IN) :: survivingIndex
        type(Tensor2) :: aMatrix
        integer :: dims(3),tempdims(3),reordering(3),newdims(2)

        if(this%Initialized) then
            !This will be the new order of the matrix
            reordering=[1,2,3]
            reordering(1)=survivingIndex
            reordering(survivingIndex)=1
            dims=shape(this%data)
            tempdims=dims
            tempdims(1)=dims(survivingIndex)
            tempdims(survivingIndex)=dims(1)
            newdims(1)=tempdims(1)
            newdims(2)=product(tempdims(2:3))
            aMatrix=new_Tensor( newdims(1), newdims(2) , ZERO )
            aMatrix%data=reshape( reshape(this%data, tempdims, ORDER=reordering), newdims)
        else
            call ThrowException('UnfoldTensor3','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function UnfoldTensor3

!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################

  subroutine SingularValueDecompositionTensor2(this,U,Sigma,vTransposed)
     class(Tensor2),intent(IN) :: this
     type(Tensor2),intent(OUT) :: U,Sigma,vTransposed
     complex(8),allocatable :: CopyOfInput(:,:) !This extra copy is here because ZGESDD destroys the input matrix
     real(8),allocatable :: DiagonalPart(:)
     integer :: LeftDimension,RightDimension
     integer :: Error=Normal
     !Lapack ugly variables
     integer :: Lwork,LRWork,LIWork,info,idx
     complex(8),allocatable :: Work(:)
     real(8),allocatable :: RWork(:)
     integer(8),allocatable :: IWork(:)
     character,parameter :: Jobz='S' !Always get the minimum only, hopefully the rest of the matrix is zeroed out

     !Prepare matrices according to input dimensions
     LeftDimension=size(this%data,1); RightDimension=size(this%data,2)
     U=new_Tensor(LeftDimension,LeftDimension,ZERO)
     Sigma=new_Tensor(LeftDimension,RightDimension,ZERO)
     vTransposed=new_Tensor(RightDimension,RightDimension,ZERO)
     allocate(DiagonalPart(min(LeftDimension,RightDimension)))
     !This doubling of memory allocation is because ZGESDD destroys the input matrix
     allocate (CopyOfInput(LeftDimension,RightDimension))
     CopyOfInput=this%data

     !Recommended values of memory allocation from LAPACK documentation
     LWork=(Min(LeftDimension,RightDimension)*(Min(LeftDimension,RightDimension)+2)+Max(LeftDimension,RightDimension))
     LRWork=5*Min(LeftDimension,RightDimension)*(Min(LeftDimension,RightDimension)+1)
     LIWork=8*Min(LeftDimension,RightDimension)

     allocate(Work(LWork),RWork(LRWork),IWork(LIWork),STAT=Error)
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Could not allocate memory',Error,CriticalError)
        return
     endif
     !For some reason I need to call LAPACK with LWork=-1 first
     !And find out the optimum work storage, otherwise it returns an error
     LWork=-1
     call ZGESDD(JOBZ, LeftDimension, RightDimension, CopyOfInput, LeftDimension, DiagonalPart, U%data, &
          & LeftDimension,vTransposed%data,RightDimension,WORK,LWORK,RWORK,IWORK,Error )
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Lapack search call returned error in ZGESDD',Error,CriticalError)
        return
     endif

     !And now call with right value of LWork
     LWork=Int(Work(1))
     deallocate(Work)
     Allocate(Work(LWork))
     call ZGESDD(JOBZ, LeftDimension, RightDimension, CopyOfInput, LeftDimension, DiagonalPart, U%data, &
          & LeftDimension,vTransposed%data,RightDimension,WORK,LWORK,RWORK,IWORK,Error )
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Lapack returned error in ZGESDD',Error,CriticalError)
        return
     endif

    !Manually insert the diagonal, move this to a routine
     do idx=1,min(LeftDimension,RightDimension)
         Sigma%data(idx,idx)=DiagonalPart(idx)
     enddo

     !Clean up: manually remove memory otherwise compiler gives warnings
     deallocate(Work,RWork,IWork,DiagonalPart,STAT=Error)
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Problems in deallocation',Error,CriticalError)
        return
     endif

   end subroutine SingularValueDecompositionTensor2

!##################################################################
!##################################################################

!##################################################################
!##################################################################
!##################################################################
!##################################################################

   function SolveLinearProblem_LAPACK_vector(theMatrix,theVector, LS_Tolerance) result(theSolution)
    class(Tensor2),intent(IN) :: theMatrix
    class(Tensor1),intent(IN) :: theVector
    real(8),intent(IN) :: LS_Tolerance
    type(Tensor1) :: theSolution

    theSolution = theVector

    call ZGELSD( theMatrix%data, theSolution%data, LS_Tolerance )

   end function SolveLinearProblem_LAPACK_vector


   function SolveLinearProblem_LAPACK_matrix(theMatrix,theVectorAsMatrix, LS_Tolerance) result(theSolution)
    class(Tensor2),intent(IN) :: theMatrix
    class(Tensor2),intent(IN) :: theVectorAsMatrix
    real(8),intent(IN) :: LS_Tolerance
    type(Tensor2) :: theSolution
    type(Tensor2) :: CopyOfMatrix
    integer :: matrixSize(2),vectorSize(2)
    !These variables are for LAPACK F77 Interface
    real(8),allocatable :: GELOutputS(:)
    integer :: rank,info,tempDim
    complex(8),allocatable :: work(:)
    real(8),allocatable :: rwork(:)
    integer,allocatable :: iwork(:)
    integer :: lwork,lrwork,liwork

    matrixSize=shape(theMatrix%data)
    vectorSize=shape(theVectorAsMatrix%data)
    theSolution=theVectorAsMatrix
    CopyOfMatrix=theMatrix

    allocate(GELoutputS(minval(matrixSize)))
    !Call MKL Lapack solver
    tempDim=1
    allocate (work(tempDim),rwork(tempDim),iwork(tempDim))
    !Find out optimum lwork
    lwork=-1
    call ZGELSD ( matrixSize(1), matrixSize(2), vectorsize(2), CopyOfMatrix%data, matrixSize(1), theSolution%data, vectorsize(1), &
        & GELOutputS, LS_Tolerance, rank, work, lwork, rwork, iwork, info)
    lwork=int(work(1))
    lrwork=int(rwork(1))
    liwork=iwork(1)
    deallocate(work,rwork,iwork)
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    !Now make final call
    call ZGELSD ( matrixSize(1), matrixSize(2), vectorsize(2), CopyOfMatrix%data, matrixSize(1), theSolution%data, vectorsize(1), &
        & GELOutputS, LS_Tolerance, rank, work, lwork, rwork, iwork, info)

   end function SolveLinearProblem_LAPACK_matrix


end module Tensor_Class