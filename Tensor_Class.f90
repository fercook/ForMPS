!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010
module Tensor_Class

  use ErrorHandling
  use Constants

  implicit none
! Need to learn how to make operators public.
  !private
  public :: new_Tensor
  public :: operator(*),assignment(=)
  public :: operator(.diff.),operator(.absdiff.)
  public :: operator(.equaldims.),operator(.equaltype.)
  public :: Conjugate,TensorTranspose,ConjugateTranspose
  public :: JoinIndicesOf,SplitIndexOf

  integer,parameter :: Max_Combined_Dimension = 100000

  type,private :: Tensor
  	private
  	integer :: Initialized=.false.
  contains
!  	procedure :: IsInitialized => Is_Tensor_init !Commented out because of Ifort bug
    procedure :: delete => delete_Tensor
  	procedure :: print => print_Tensor
	procedure :: PrintDimensions => Print_Tensor_Dimensions
    procedure :: getDimensions => getDimensions_Of_Tensor
    procedure :: Norm => Norm_Of_Tensor
  end type Tensor

  type,public,extends(Tensor) :: Tensor1
  	private
  	complex(8),allocatable :: data(:)
  end type Tensor1

  type,public,extends(Tensor) :: Tensor2
  	private
  	complex(8),allocatable :: data(:,:)
  contains
    procedure :: SVD => SingularValueDecomposition
    procedure :: SplitIndex => SplitIndexOfTensor2
!    procedure :: dagger => ConjugateTranspose
  end type Tensor2

  type,public,extends(Tensor) :: Tensor3
  	private
  	complex(8),allocatable :: data(:,:,:)
  contains
    procedure :: JoinIndices => JoinIndicesOfTensor3
  end type Tensor3

  type,public,extends(Tensor) :: Tensor4
    private
    complex(8),allocatable :: data(:,:,:,:)
  contains
    procedure :: JoinIndices => JoinIndicesOfTensor4
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
     	  & Tensor1_dotProduct_Tensor1
  end interface

  interface assignment (=)
     module procedure new_Tensor1_fromAssignment,  &
          & new_Tensor2_fromAssignment, new_Tensor3_fromAssignment, new_Tensor4_fromAssignment
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

  interface JoinIndicesOf
  	module procedure JoinIndicesOfTensor3,JoinIndicesOfTensor4
  end interface

  interface SplitIndexOf
    module procedure SplitIndexOfTensor2
  end interface

  interface Conjugate
    module procedure ConjugateTensor1,ConjugateTensor2,ConjugateTensor3,ConjugateTensor4
  end interface

  interface TensorTranspose
    module procedure TensorTranspose2,TensorTranspose3 !,TensorTranspose4
  end interface

  interface ConjugateTranspose
    module procedure ConjugateTranspose2,ConjugateTranspose3
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
        call ThrowException('new_Tensor2_withConstant','Dimensions are larger than maximum',NoErrorCode,CriticalError)
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
        call ThrowException('new_Tensor4_withConstant','Dimensions are larger than maximum',NoErrorCode,CriticalError)
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
     type(Tensor1),intent(in) :: tensor
     type(Tensor1) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor1_fromTensor1

   function new_Tensor2_fromTensor2 (tensor) result (this)
     type(Tensor2),intent(in) :: tensor
     type(Tensor2) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor2_fromTensor2

   function new_Tensor3_fromTensor3 (tensor) result (this)
     type(Tensor3),intent(in) :: tensor
     type(Tensor3) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2),size(tensor%data,3)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor3_fromTensor3

   function new_Tensor4_fromTensor4 (tensor) result (this)
     type(Tensor4),intent(in) :: tensor
     type(Tensor4) this
     integer error

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2),size(tensor%data,3),size(tensor%data,4)))
     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor4_fromTensor4

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
   integer function delete_Tensor (this) result(error)
     class(Tensor),intent(INOUT) :: this   !!<<TYPE>>!!

     error=Warning

     if(.not.this%Initialized) then
        call ThrowException('delete_Tensor','Trying to delete an uninitialized tensor',NoErrorCode,error)
        return
     endif

     ! Need to check the type to deallocate memory
	 select type (Typed_this => this)
	 	class is (Tensor1)
	 		deallocate(Typed_this%data)
	 	class is (Tensor2)
	 		deallocate(Typed_this%data)
	 	class is (Tensor3)
	 		deallocate(Typed_this%data)
        class is (Tensor4)
            deallocate(Typed_this%data)
	 	class is (Tensor)
	 end select

     !Flip flag
     this%Initialized=.false.

     error=Normal

   end function delete_Tensor
!##################################################################


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
   subroutine Print_Tensor(this,error)
     class(Tensor),intent(IN) :: this
     integer i,j,k
     integer,optional :: error

     If(present(error)) error = Warning

     if(.not.(this%Initialized)) then
        call ThrowException('PrintTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     endif

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


   subroutine Print_Tensor_Dimensions(this)
     class(Tensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     if(.not.(this%Initialized)) then
        call ThrowException('Print Tensor','Tensor not initialized',NoErrorCode,Warning)
        return
     endif

	 select type (Typed_this => this)
	 	class is (Tensor1)
        	print *,'Vector dimension:',size(Typed_this%data,1)
	 	class is (Tensor2)
        	print *,'Matrix dimensions:',size(Typed_this%data,1),'x',size(Typed_this%data,2)
	 	class is (Tensor3)
        	print *,'3-Tensor data:',size(Typed_this%data,1),'x',size(Typed_this%data,2), &
        		& 'x',size(Typed_this%data,3)
        class is (Tensor4)
            print *,'4-Tensor data:',size(Typed_this%data,1),'x',size(Typed_this%data,2), &
                & 'x',size(Typed_this%data,3),'x',size(Typed_this%data,4)
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
!	 			Dims=new_Tensor1_fromData(one*[size(Typed_this%data,1)])
	 		class is (Tensor2)
                allocate(Dims(2))
                Dims=shape(Typed_this%data)
!	 			Dims=new_Tensor1_fromData(one*[size(Typed_this%data,1),size(Typed_this%data,2)] )
	 		class is (Tensor3)
                allocate(Dims(3))
                Dims=shape(Typed_this%data)
!	 			Dims=new_Tensor1_fromData(one*[size(Typed_this%data,1),size(Typed_this%data,2),size(Typed_this%data,3)] )
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

!!************************
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
!!************************

   function number_times_Tensor1(constant, aTensor) result(this)
     complex(8),intent(IN) :: constant
     type(Tensor1),intent(IN) :: aTensor
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
     type(Tensor2),intent(IN) :: aTensor
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
     type(Tensor3),intent(IN) :: aTensor
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
     type(Tensor4),intent(IN) :: aTensor
     type(Tensor4) :: this

     if(.not.aTensor%Initialized) then
        call ThrowException('Number_times_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
     else
         this=new_Tensor(constant*aTensor%data)
     endif
     return

   end function Number_times_Tensor4

   function Tensor2_matmul_Tensor2(tensorA,tensorB) result(this)
      class(Tensor2),intent(IN) :: tensorA,tensorB
      type(Tensor2) :: this

      if(TensorA%Initialized.and.TensorB%Initialized) then
         this=new_Tensor(matmul(tensorA%data,tensorB%data))
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

!!##################################################################
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
    equals=.true.  !!TODO

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

!
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


function SplitIndexOfTensor2(this,WhichIndex,Partition) result (aTensor)
    integer,intent(IN) :: WhichIndex(:),Partition
    class(tensor2),intent(IN) :: this
    type(tensor3) :: aTensor
    integer :: i,NewDims(3),OrigDims(2),error

     if(.not.(this%Initialized)) then
        call ThrowException('SplitIndexOfTensor2','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif
     if( (.not.(WhichIndex.equalvector.FIRST)).and.(.not.(WhichIndex.equalvector.SECOND)) ) then
        call ThrowException('SplitIndexOfTensor2','Index is inappropriate',WhichIndex(1),CriticalError)
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

  subroutine SingularValueDecomposition(this,U,Sigma,vTransposed,ErrorCode)
     class(Tensor2),intent(IN) :: this
     type(Tensor2),intent(OUT) :: U,Sigma,vTransposed
     integer,intent(OUT),optional :: ErrorCode
     complex(8),allocatable :: CopyOfInput(:,:) !This is here because ZGESDD destroys the input matrix
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
     !This is here because ZGESDD destroys the input matrix
     allocate (CopyOfInput(LeftDimension,RightDimension))
     CopyOfInput=this%data

     !Recommended values of memory allocation from LAPACK documentation
     LWork=(Min(LeftDimension,RightDimension)*(Min(LeftDimension,RightDimension)+2)+Max(LeftDimension,RightDimension))
     LRWork=5*Min(LeftDimension,RightDimension)*(Min(LeftDimension,RightDimension)+1)
     LIWork=8*Min(LeftDimension,RightDimension)

     allocate(Work(LWork),RWork(LRWork),IWork(LIWork),STAT=Error)
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Could not allocate memory',Error,CriticalError)
        if(present(ErrorCode)) ErrorCode=Error
        return
     endif
     !For some reason I need to call LAPACK with LWork=-1 first
     !And find out the optimum work storage, otherwise it returns an error
     LWork=-1
     call ZGESDD(JOBZ, LeftDimension, RightDimension, CopyOfInput, LeftDimension, DiagonalPart, U%data, &
          & LeftDimension,vTransposed%data,RightDimension,WORK,LWORK,RWORK,IWORK,Error )
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Lapack search call returned error in ZGESDD',Error,CriticalError)
        if(present(ErrorCode)) ErrorCode=Error
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
        if(present(ErrorCode)) ErrorCode=Error
        return
     endif

    !Manually insert the diagonal, move this to a routine
     do idx=1,min(LeftDimension,RightDimension)
         Sigma%data(idx,idx)=DiagonalPart(idx)
     enddo

     !Clean up
     deallocate(Work,RWork,IWork,DiagonalPart,STAT=Error)
     If (Error.ne.Normal) then
        call ThrowException('SingularValueDecomposition','Problems in deallocation',Error,CriticalError)
        if(present(ErrorCode)) ErrorCode=Error
        return
     endif

     if(present(ErrorCode)) ErrorCode=Normal

   end subroutine SingularValueDecomposition




end module Tensor_Class
