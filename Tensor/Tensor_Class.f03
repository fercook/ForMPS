!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010
module Tensor_Class

  use ErrorHandling
  use Constants
  use Matrix_Helper

  implicit none
  private

  integer,parameter :: Max_Combined_Dimension = 100000

  type,private :: Tensor
  	private
  	integer :: Initialized=.false.
  contains
!  	procedure :: IsInitialized => Is_Tensor_init
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
  end type Tensor2

  type,public,extends(Tensor) :: Tensor3
  	private
  	complex(8),allocatable :: data(:,:,:)
  end type Tensor3

!###############################
!#####  Operators and methods
!###############################

  interface new_Tensor
     module procedure new_Tensor1_Random,new_Tensor1_fromTensor1,new_Tensor1_fromData,new_Tensor1_withConstant, &
		& new_Tensor2_Random,new_Tensor2_fromTensor2,new_Tensor2_fromData,new_Tensor2_withConstant, &
		& new_Tensor3_Random,new_Tensor3_fromTensor3,new_Tensor3_fromData,new_Tensor3_withConstant
  end interface

  interface operator (*)
     module procedure &
     	  & number_times_Tensor1,number_times_Tensor2,number_times_Tensor3
  end interface

  interface assignment (=)
     module procedure new_Tensor1_fromAssignment,  &
          & new_Tensor2_fromAssignment, new_Tensor3_fromAssignment !,&
          !& new_Tensor4_fromAssignment, new_Tensor5_fromAssignment
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

  interface JoinIndices
  	module procedure JoinIndicesOfTensor3
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

     !initialize internal variables
     !initialize data
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

     if(dim1*dim2.lt.Max_Combined_Dimension) then
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

     if(dim1*dim2*dim3.lt.Max_Combined_Dimension) then
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

!##################################################################
   function new_Tensor1_fromData (originalData) result (this)
     complex(8),intent(in) :: originalData(:)
     integer :: dim1
     type(Tensor1) this

     dim1=size(originalData,1)

     if(dim1.lt.Max_Combined_Dimension) then
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

     if(dim1*dim2.lt.Max_Combined_Dimension) then
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
     if(dim1*dim2*dim3.lt.Max_Combined_Dimension) then
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

!##################################################################
   function new_Tensor1_withConstant (dim1,constant) result (this)
     integer,intent(in) :: dim1
     complex(8),intent(in) :: constant
     type(Tensor1) this

     if(dim1.lt.Max_Combined_Dimension) then
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

     if(dim1*dim2.lt.Max_Combined_Dimension) then
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

     if(dim1*dim2*dim3.lt.Max_Combined_Dimension) then
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
   integer function Print_Tensor(this) result(error)
     class(Tensor),intent(IN) :: this
     integer i,j,k

     error = Warning

     if(.not.(this%Initialized)) then
        call ThrowException('PrintTensor','Tensor not initialized',NoErrorCode,error)
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
	 	class is (Tensor)
	 		print *,'No data in raw tensor'
	 end select

     error=Normal

   end function Print_Tensor


   integer function Print_Tensor_Dimensions(this) result(error)
     class(Tensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     error = Warning

     if(.not.(this%Initialized)) then
        call ThrowException('Print Tensor','Tensor not initialized',NoErrorCode,error)
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
	 	class is (Tensor)
	 		print *,'No data in raw tensor'
	 end select

     error=Normal

   end function Print_Tensor_Dimensions
!##################################################################


	function getDimensions_Of_Tensor(this) result(Dims)
		class(Tensor),intent(IN) :: this
		type(Tensor1) :: Dims

	    if(.not.(this%Initialized)) then
    	    call ThrowException('getDimensions','Tensor not initialized',NoErrorCode,CriticalError)
        	return
	    endif

	 	select type (Typed_this => this)
	 		class is (Tensor1)
	 			Dims=new_Tensor1_fromData(one*[size(Typed_this%data,1)])
	 		class is (Tensor2)
	 			Dims=new_Tensor1_fromData(one*[size(Typed_this%data,1),size(Typed_this%data,2)] )
	 		class is (Tensor3)
	 			Dims=new_Tensor1_fromData(one*[size(Typed_this%data,1),size(Typed_this%data,2),size(Typed_this%data,3)] )
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
    equals=.true.

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
!!#######################################################################################
!!#######################################################################################
!! This are very important functions as most of the algorithm time is spent updating the
!! matrices using this Left and Right products with tensors.
!! They are also used heavily for computing expectation values, so optimization here might be key.
!! Convention:
!!              L(betaR,alphaR) = \sum_i B_i^\dagger . L_in A_i
!!                              = \sum_i \sum_betaL \sum_alphaL  B^*_{i,betaR,betaL} Lin_{betaL,alphaL} A_{i,alphaL,alphaR)
!!
!!              R(alphaL,betaL) = \sum_i A_i . RL_in . B_i^\dagger
!!                              = \sum_i \sum_alphaR \sum_betaR A_{i,alphaL,alphaR) Rin_{alphaR,betaR} B^*_{i,betaR,betaL}
!!
!!#######################################################################################
!!#######################################################################################
!
!   function MPS_Left_Product(TensorA,TensorB,matrixin) result(matrixout)
!     type(MPSTensor),intent(IN) :: TensorA,TensorB
!     type(MPSTensor) :: matrixout
!     type(MPSTensor),intent(IN),optional :: matrixin
!     type(MPSTensor) :: TempMatrix,L_in_matrix
!     integer :: s,i,k,j,l
!     complex(8) :: TEMP
!
!     if((.not.tensorA%Initialized).and.(.not.tensorB%Initialized)) then
!        call ThrowException('MPSLeftProduct','Tensors not initialized',NoErrorCode,CriticalError)
!        return
!     endif
!     if (TensorA%spin.ne.TensorB%spin) then
!        call ThrowException('MPSLeftProduct','Tensors have different spin',NoErrorCode,CriticalError)
!        return
!     endif
!     if (present(matrixin)) then
!        if(matrixin%Initialized) then
!           L_in_matrix=new_MPSTensor(matrixin)
!        else
!           call ThrowException('MPSLeftProduct','Matrix is not initialized',NoErrorCode,CriticalError)
!           return
!        endif
!     else
!        L_in_matrix=new_MPSTensor(MatrixSpin,TensorB%DLeft,TensorA%DLeft,one)
!     endif
!
!     matrixout=new_MPSTensor(MatrixSpin, TensorB%DRight,TensorA%DRight, zero)
!     TempMatrix=new_MPSTensor(MatrixSpin,TensorB%DLeft ,TensorA%DRight, zero)
!
!     !The multiplications are done by hand because I could not get ZGEMM to work properly
!     do s=1,TensorA%spin
!        Tempmatrix%data(:,:,MatrixSpin)=matmul(L_in_matrix%data(:,:,MatrixSpin),TensorA%data(:,:,s))
!        MatrixOut%data(:,:,MatrixSpin)=MatrixOut%data(:,:,MatrixSpin)+matmul(dconjg(transpose(TensorB%data(:,:,s))),Tempmatrix%data(:,:,MatrixSpin))
!    enddo
!    return
!  end function MPS_Left_Product
!
!
!!##################################################################
!!##################################################################
!
!  function MPS_Right_Product(TensorA,TensorB,matrixin) result(matrixout)
!    type(MPSTensor),intent(IN) :: TensorA,TensorB
!    type(MPSTensor) :: matrixout
!    type(MPSTensor),intent(IN),optional :: matrixin
!    type(MPSTensor) :: TempMatrix,R_in_matrix
!    integer :: s,i,k,j,l
!    complex(8) :: TEMP
!
!    if((.not.tensorA%Initialized).and.(.not.tensorB%Initialized)) then
!       call ThrowException('MPSRightProduct','Tensors not initialized',NoErrorCode,CriticalError)
!       return
!    endif
!    if (TensorA%spin.ne.TensorB%spin) then
!       call ThrowException('MPSRightProduct','Tensors have different spin',NoErrorCode,CriticalError)
!       return
!    endif
!
!    matrixout=new_MPSTensor(MatrixSpin, TensorA%DLeft,TensorB%DLeft, zero)
!    TempMatrix=new_MPSTensor(MatrixSpin,TensorA%DLeft ,TensorB%DRight, zero)
!
!    if (present(matrixin)) then
!       if(matrixin%Initialized) then
!          R_in_matrix=new_MPSTensor(matrixin)
!       else
!          call ThrowException('MPSRightProduct','Matrix is not initialized',NoErrorCode,CriticalError)
!          return
!       endif
!    else
!       R_in_matrix=new_MPSTensor(MatrixSpin,TensorA%DRight,TensorB%DRight,one)
!    endif
!
!    !The multiplications are done by hand because I could not get ZGEMM to work properly
!    do s=1,TensorA%spin
!       Tempmatrix%data(:,:,MatrixSpin)=matmul(TensorA%data(:,:,s),R_in_matrix%data(:,:,MatrixSpin))
!       MatrixOut%data(:,:,MatrixSpin)=MatrixOut%data(:,:,MatrixSpin)+matmul(Tempmatrix%data(:,:,MatrixSpin),transpose(dconjg(TensorB%data(:,:,s))))
!    enddo
!    return
!  end function MPS_Right_Product
!
!
!aTensor%JoinIndices(FirstAndSecond,Third)
!aTensor%JoinIndices(Second,FirstAndThird)
!aTensor%JoinIndices(ThirdAndSecond,First)
!aTensor%JoinIndices(ThirdAndFirst,Second)
!aTensor%JoinIndices(SecondAndFirst,Third)

function JoinIndicesOfTensor3(this,firstindex,secondindex) result (aTensor)
	integer,intent(IN) :: firstindex,secondindex
	type(tensor3),intent(IN) :: this
	type(tensor2) :: aTensor
	integer :: i,dim(3)
	logical :: isNormalOrder
	complex(8),allocatable :: anArray(:,:)

     if(.not.(this%Initialized)) then
        call ThrowException('JoinIndicesOfTensor3','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

	dim=size(this%data)

	if (firstindex .eq. FIRST .or. secondindex .eq. FIRST) then
		allocate(anArray(dim(1),dim(2)*dim(3)))
		if (secondindex .eq. SECONDANDTHIRD .or. firstindex .eq. SECONDANDTHIRD ) then
			do i=1,dim(1)
				anArray(i,:)=reshape( this%data(i,:,:), [dim(2)*dim(3)])
			enddo
			isNormalOrder=.true.
		else if (secondindex .eq. THIRDANDSECOND .or. firstindex .eq. THIRDANDSECOND) then
			do i=1,dim(1)
				anArray(i,:)=reshape( transpose(this%data(i,:,:)), [dim(2)*dim(3)])
			enddo
			isNormalOrder=.false.
		else
	        call ThrowException('JoinIndicesOfTensor3','Incorrect joint index',secondindex-firstindex,CriticalError)
    	    return
		endif

	else if(firstindex .eq. SECOND .or. secondindex .eq. SECOND) then
		allocate(anArray(dim(2),dim(1)*dim(3)))
		if (secondindex .eq. FIRSTANDTHIRD .or. firstindex .eq. FIRSTANDTHIRD ) then
			do i=1,dim(2)
				anArray(i,:)=reshape( this%data(:,i,:), [dim(1)*dim(3)])
			enddo
			isNormalOrder=.true.
		else if (secondindex .eq. THIRDANDFIRST .or. firstindex .eq. THIRDANDFIRST) then
			do i=1,dim(2)
				anArray(i,:)=reshape( transpose(this%data(:,i,:)), [dim(1)*dim(3)])
			enddo
			isNormalOrder=.false.
		else
	        call ThrowException('JoinIndicesOfTensor3','Incorrect joint index',secondindex-firstindex,CriticalError)
    	    return
		endif

	else if(firstindex .eq. THIRD .or. secondindex .eq. THIRD) then
		allocate(anArray(dim(3),dim(1)*dim(2)))
		if (secondindex .eq. FIRSTANDSECOND .or. firstindex .eq. FIRSTANDSECOND ) then
			do i=1,dim(3)
				anArray(i,:)=reshape( this%data(:,:,i), [dim(1)*dim(2)])
			enddo
			isNormalOrder=.true.
		else if (secondindex .eq. SECONDANDFIRST .or. firstindex .eq. SECONDANDFIRST) then
			do i=1,dim(3)
				anArray(i,:)=reshape( transpose(this%data(:,:,i)), [dim(1)*dim(2)])
			enddo
			isNormalOrder=.false.
		else
	        call ThrowException('JoinIndicesOfTensor3','Incorrect joint index',secondindex-firstindex,CriticalError)
    	    return
		endif
	else
        call ThrowException('JoinIndicesOfTensor3','Incorrect single index',secondindex-firstindex,CriticalError)
    	return
	endif

	!Finally allocate the tensor checking for the reverse order or normal order
	if (isNormalOrder) then
		aTensor=new_Tensor(anArray)
	else
		aTensor=new_Tensor(Transpose(anArray))
	endif

end function JoinIndicesOfTensor3
!
!  subroutine CollapseSpinWithBond(this,collapsed,whichDimension)
!    class(MPSTensor),intent(IN) :: this
!    complex(8),intent(OUT) :: collapsed(:,:)
!    integer,intent(IN) :: whichDimension
!    integer :: s,alpha,beta,leftIndex,rightIndex,leftStep,rightStep,leftDimension,rightDimension
!
!    if(.not.this%Initialized) then
!       call ThrowException('CollapseSpinWithBond','Tensor not initialized',NoErrorCode,CriticalError)
!       return
!    endif
!
!    if (whichDimension.eq.FirstDimension) then
!       leftStep=this%DLeft
!       rightStep=0
!       leftDimension=(this%spin*this%DLeft)
!       rightDimension=(this%DRight)
!    else if (whichDimension.eq.SecondDimension) then
!       leftStep=0
!       rightStep=this%DRight
!       leftDimension=(this%DLeft)
!       rightDimension=(this%spin*this%DRight)
!    else
!       call ThrowException('CollapseSpinWithBond','Wrong Dimension parameter',whichDimension,CriticalError)
!       return
!    endif
!    if ((size(collapsed,1).ne.leftDimension).and.(size(collapsed,2).ne.rightDimension)) then
!       call ThrowException('CollapseSpinWithBond','Matrix for collapsed tensor does not have right dimensions' &
!            & ,whichDimension,CriticalError)
!       return
!    endif
!
!    !This always puts the spin before the bond dimension,
!    !      [(s,alpha),(beta)]   or  [(alpha),(s,beta)]
!    do s=1,this%spin
!       do beta=1,this%DRight
!          rightIndex=beta+(s-1)*rightStep
!          do alpha=1,this%DLeft
!             leftIndex=alpha+(s-1)*leftStep
!             collapsed(leftIndex,rightIndex)=this%data(alpha,beta,s)
!          enddo
!       enddo
!    enddo
!
!  end subroutine CollapseSpinWithBond
!


!  function new_Tensor3_fromTensor2(dim1,dim2,dim3,matrix) result(this)
!    type(Tensor3) :: this
!    type(Tensor2),intent(IN) :: matrix
!    integer,intent(IN) :: dim1,dim2,dim3
!    integer :: alpha,beta,s
!    integer :: leftIndex,rightIndex,leftStep,rightStep,leftDimension,rightDimension
!    character(100) :: Message
!    character(3) :: ScratchMessage
!
!    if (size(matrix%data,1).eq.dim1*dim2) then
!       leftStep=dim2
!       rightStep=0
!    else if (size(matrix%data,2).eq.dim1*dim3) then
!       leftStep=0
!       rightStep=dim3
!    else
!       !TODO: Awful code follows, perhaps implement a message class or something that joins chars and nums
!       write (ScratchMessage,'(I3)') spin
!       Message='s='//trim(adjustl(ScratchMessage))
!       write (ScratchMessage,'(I3)') LeftBond
!       Message=trim(adjustl(Message))//', DL:'//trim(adjustl(ScratchMessage))
!       write (ScratchMessage,'(I3)') RightBond
!       Message=trim(adjustl(Message))//', DR:'//trim(adjustl(ScratchMessage))
!       write (ScratchMessage,'(I3)') size(matrix,1)
!       Message=trim(adjustl(Message))//'; received:'//trim(adjustl(ScratchMessage))
!       write (ScratchMessage,'(I3)') size(matrix,2)
!       Message=trim(adjustl(Message))//' x '//trim(adjustl(ScratchMessage))
!       call ThrowException('new_Tensor3_fromTensor2','Wrong dimensions='//Message, &
!            & NoErrorCode,CriticalError)
!       return
!    endif
!
!    if(this%Initialized) deallocate(this%data)
!    allocate(this%data(dim1,dim2,dim3))
!
!    this%data=zero
!    do s=1,Spin
!       do beta=1,RightBond
!          rightIndex=beta+(s-1)*rightStep
!          do alpha=1,LeftBond
!             leftIndex=alpha+(s-1)*leftStep
!             this%data(alpha,beta,s)=matrix(leftIndex,rightIndex)
!          enddo
!       enddo
!    enddo
!
!    this%Initialized=.true.
!
!  end function new_Tensor3_fromMatrix

!!##################################################################
!!###  IMPORTANT NOTE ON THE INTERFACE OF Matrix_times_MPSTensor:
!!###  One of the arguments must be a "Tmatrix", i.e. a MPSTensor with spin=MatrixSpin=1
!!###  The order of the arguments is VERY important
!!###  If the first argument is a matrix, then the result is matrix.tensor
!!###  otherwise the result is tensor.matrix
!!###
!   function Matrix_times_MPSTensor(tensorA, tensorB) result(this)
!     type(MPSTensor),intent(IN) :: tensorA,tensorB
!     type(MPSTensor) this
!     integer :: s
!
!     if(tensorA%Initialized.and.tensorB%Initialized) then
!        if(tensorA%DRight.eq.tensorB%DLeft) then
!           !The trick of using a tensor as a matrix is used here:
!           if (tensorA%spin.eq.MatrixSpin) then
!              this = new_MPSTensor(tensorB%spin,tensorA%DLeft,tensorB%DRight,zero)
!              do s=1,tensorB%spin
!                 this%data(:,:,s)=matmul(tensorA%data(:,:,MatrixSpin),tensorB%data(:,:,s)) !+this%data(:,:,s)
!              enddo
!           else if (tensorB%spin.eq.MatrixSpin) then
!              this = new_MPSTensor(tensorA%spin,tensorA%DLeft,tensorB%DRight,zero)
!              do s=1,tensorA%spin
!                 this%data(:,:,s)=matmul(tensorA%data(:,:,s),tensorB%data(:,:,MatrixSpin)) !+this%data(:,:,s)
!              enddo
!           else
!              call ThrowException('Matrix_times_MPSTensor','One of the arguments must be a *matrix* (MPSTensor with spin 1)' &
!                   & ,NoErrorCode,CriticalError)
!           endif
!        else
!           call ThrowException('Matrix_times_MPSTensor','Dimensions of the tensors do not match',NoErrorCode,CriticalError)
!        endif
!     else
!        call ThrowException('Matrix_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
!     endif
!     return
!   end function Matrix_times_MPSTensor
!


end module Tensor_Class
