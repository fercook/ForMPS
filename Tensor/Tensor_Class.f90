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
  	procedure :: IsInitialized => Is_Tensor_init
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
     module procedure new_Tensor1_Random,new_Tensor1_fromTensor1,new_Tensor1_fromData,new_Tensor1_Constant, &
		& new_Tensor2_Random,new_Tensor2_fromTensor2,new_Tensor2_fromData,new_Tensor2_Constant, &
		& new_Tensor3_Random,new_Tensor3_fromTensor3,new_Tensor3_fromData,new_Tensor3_Constant
  end interface

  interface operator (*)
     module procedure &
     	  & Integer_times_Tensor1,Integer_times_Tensor2,Integer_times_Tensor3, &
          & Real_times_Tensor1,Real_times_Tensor2, Real_times_Tensor3, &
          & Complex_times_Tensor1, Complex_times_Tensor2, Complex_times_Tensor3, &
          & Real8_times_Tensor1, Real8_times_Tensor2, Real8_times_Tensor3, &
          & Real8_times_Tensor4, Real8_times_Tensor5,  &
          & Complex8_times_Tensor1, Complex8_times_Tensor2, Complex8_times_Tensor3
!          & Tensor2_times_Tensor3 !&
!         & Integer_times_Tensor4,Integer_times_Tensor5, &
!         & Real_times_Tensor4, Real_times_Tensor5, &
!         & Complex_times_Tensor4, Complex_times_Tensor5, &
!         & Complex8_times_Tensor4, Complex8_times_Tensor5
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

     !initialize internal variables
     !initialize data
     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data=randomtensorR+II*randomtensorC

     this%Initialized=.true.

   end function new_Tensor3_Random

!##################################################################
   function new_Tensor3_withData (originalData) result (this)
     complex(8),intent(in) :: originalData(:,:,:)
     integer :: dim1,dim2,dim3
     type(Tensor3) this

	dim1=size(originalData,1)
	dim2=size(originalData,2)
	dim3=size(originalData,3)

     if(dim1*dim2*dim3.lt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor3_withData','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1) then
        call ThrowException('new_Tensor3_withData','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))

     this%data=originalData
     this%Initialized=.true.

   end function new_Tensor3_withData

!##################################################################
   function new_Tensor3_withConstant (dim1,dim2,dim3,constant) result (this)
     integer,intent(in) :: dim1,dim2,dim3
     complex(8),intent(in) :: constant
     type(Tensor3) this

     if(dim1*dim2*dim3.lt.Max_Combined_Dimension) then
        call ThrowException('new_Tensor3_withData','Dimensions are larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(dim1.lt.1.or.dim2.lt.1.or.dim3.lt.1) then
        call ThrowException('new_Tensor3_withData','One dimension is smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(dim1,dim2,dim3))

     this%data=constant
     this%Initialized=.true.

   end function new_Tensor3_withConstant

!##################################################################
   function new_MPSTensor3_fromMPSTensor3 (tensor) result (this)
     type(Tensor3),intent(in) :: tensor
     type(Tensor3) this
     integer error

     error=tensor%isInitialized()
     if (WasThereError()) call ProcessException('new_Tensor3_fromTensor3')

     if(this%Initialized) deallocate(this%data)
     allocate(this%data(size(tensor%data,1),size(tensor%data,2),size(tensor%data,3)))

     this%data=tensor%data
     this%Initialized=.true.

   end function new_Tensor3_fromTensor3

   subroutine new_Tensor3_fromAssignment(lhs,rhs)
     class(Tensor3),intent(out) :: lhs
     type(Tensor3),intent(in) :: rhs

     if(.not.rhs%Initialized) then
        call ThrowException('new_Tensor3_fromAssignment','Original tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

!    The following check should be
!          if(lhs%Initialized) deallocate(lhs%data)
!    But instead I check for allocated because of a bug in GFortran,
!    http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43969
     if(allocated(lhs%data)) deallocate(lhs%data)

     allocate(lhs%data(size(rhs%data,1),size(rhs%data,1),size(rhs%data,1)))

     lhs%data=rhs%data
     lhs%Initialized=.true.

   end subroutine new_Tensor3_fromAssignment





  function new_Tensor3_fromArray(dim1,dim2,dim3,matrix) result(this)
    type(Tensor3) :: this
    complex(8),intent(IN) :: matrix(:,:)
    integer,intent(IN) :: dim1,dim2,dim3
    integer :: alpha,beta,s,spin
    integer :: leftIndex,rightIndex,leftStep,rightStep,leftDimension,rightDimension
    character(100) :: Message
    character(3) :: ScratchMessage

    if (size(matrix,1).eq.dim1*dim2) then
       leftStep=dim2
       rightStep=0
    else if (size(matrix,2).eq.dim1*dim3) then
       leftStep=0
       rightStep=dim3
    else
       !TODO: Awful code follows, perhaps implement a message class or something that joins chars and nums
       write (ScratchMessage,'(I3)') spin
       Message='s='//trim(adjustl(ScratchMessage))
       write (ScratchMessage,'(I3)') LeftBond
       Message=trim(adjustl(Message))//', DL:'//trim(adjustl(ScratchMessage))
       write (ScratchMessage,'(I3)') RightBond
       Message=trim(adjustl(Message))//', DR:'//trim(adjustl(ScratchMessage))
       write (ScratchMessage,'(I3)') size(matrix,1)
       Message=trim(adjustl(Message))//'; received:'//trim(adjustl(ScratchMessage))
       write (ScratchMessage,'(I3)') size(matrix,2)
       Message=trim(adjustl(Message))//' x '//trim(adjustl(ScratchMessage))
       call ThrowException('new_MPSTensor_fromMatrix','Wrong dimensions='//Message, &
            & NoErrorCode,CriticalError)
       return
    endif

    if(this%Initialized) deallocate(this%data)
    allocate(this%data(dim1,dim2,dim3))

    this%data=zero
!    do s=1,Spin
!       do beta=1,RightBond
!          rightIndex=beta+(s-1)*rightStep
!          do alpha=1,LeftBond
!             leftIndex=alpha+(s-1)*leftStep
!             this%data(alpha,beta,s)=matrix(leftIndex,rightIndex)
!          enddo
!       enddo
!    enddo

    this%Initialized=.true.

  end function new_Tensor3_fromMatrix

!######################################    delete
   integer function delete_MPSTensor (this) result(error)
     !!class(MPSTensor),intent(INOUT) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(INOUT) :: this   !!<<TYPE>>!!

     error=Warning

     if(.not.this%Initialized) then
        call ThrowException('delete_MPSTensor','Trying to delete an uninitialized tensor',NoErrorCode,error)
        return
     endif

     !Erase info
     this%spin=0
     this%DLeft=0
     this%DRight=0
     !Erase data
     deallocate(this%data)
     !Flip flag
     this%Initialized=.false.

     error=Normal

   end function delete_MPSTensor
!##################################################################


integer function InitializationCheck(this) result(error)
    !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
    class(MPSTensor),intent(IN) :: this !!<<TYPE>>!!

    if (.not.this%Initialized) then
       error=CriticalError
       call ThrowException('Internal Routine ','Uninitialized tensor',NoErrorCode,error)
    else
       error=Normal
    endif

  end function InitializationCheck


!######################################     print
   integer function Print_MPSTensor(this) result(error)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     error = Warning

     if(.not.(this%Initialized)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,error)
        return
     endif

     do i=1,this%spin
        print *,'State :',i
        do j=1,this%DLeft
           print *,(this%data(j,k,i),k=1,this%DRight)
        enddo
     enddo

     error=Normal

   end function Print_MPSTensor



   integer function Print_MPSTensor_Dimensions(this) result(error)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     error = Warning

     if(.not.(this%Initialized)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,error)
        return
     endif

     print *,'Spin = ',this%spin
     print *,'DL = ',this%DLeft
     print *,'DR = ',this%DRight

     error=Normal

   end function Print_MPSTensor_Dimensions
!##################################################################

!##################################################################
!###########       Accessor methods
!##################################################################
   integer function spin_MPSTensor(this) result(s)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!

    if(.not.(this%Initialized)) then
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

     if(.not.(this%Initialized)) then
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

     if(.not.(this%Initialized)) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DR=this%DRight
     endif

   end function DRight_MPSTensor

!##################################################################
!#######################        Products by things
!##################################################################

   real(8) function Norm_Of_MPSTensor(this)
     class(MPSTensor),intent(IN) :: this   !!<<TYPE>>!!
     integer :: s,alpha,beta

     Norm_Of_MPSTensor=0.0d0
     do s=1,this%spin
        do beta=1,this%DRight
           do alpha=1,this%DLeft
              Norm_Of_MPSTensor=Norm_Of_MPSTensor+abs(this%data(alpha,beta,s))
           enddo
        enddo
     enddo

     return
   end function Norm_Of_MPSTensor


   function Integer_times_MPSTensor(constant, tensor) result(this)
     integer,intent(IN) :: constant
     type(MPSTensor),intent(IN) :: tensor
     type(MPSTensor) this

     if(tensor%Initialized) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data)
        return
     else
        call ThrowException('Integer_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Integer_times_MPSTensor

!##################################################################
   function Real_times_MPSTensor(constant, tensor) result(this)
     real,intent(IN) :: constant
     type(MPSTensor),intent(IN) :: tensor
     type(MPSTensor) this

     if(tensor%Initialized) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data)
        return
     else
        call ThrowException('Real_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Real_times_MPSTensor

!##################################################################
   function Real8_times_MPSTensor(constant, tensor) result(this)
     real(8),intent(IN) :: constant
     type(MPSTensor),intent(IN) :: tensor
     type(MPSTensor) this

     if(tensor%Initialized) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data)
        return
     else
        call ThrowException('Real8_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Real8_times_MPSTensor

!##################################################################
   function Complex_times_MPSTensor(constant, tensor) result(this)
     complex,intent(IN) :: constant
     type(MPSTensor),intent(IN) :: tensor
     type(MPSTensor) this

     if(tensor%Initialized) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data)
        return
     else
        call ThrowException('Complex_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Complex_times_MPSTensor

!##################################################################
   function Complex8_times_MPSTensor(constant, tensor) result(this)
     complex(8),intent(IN) :: constant
     type(MPSTensor),intent(IN) :: tensor
     type(MPSTensor) this

     if(tensor%Initialized) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data)
        return
     else
        call ThrowException('Complex8_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Complex8_times_MPSTensor


!##################################################################
!###  IMPORTANT NOTE ON THE INTERFACE OF Matrix_times_MPSTensor:
!###  One of the arguments must be a "Tmatrix", i.e. a MPSTensor with spin=MatrixSpin=1
!###  The order of the arguments is VERY important
!###  If the first argument is a matrix, then the result is matrix.tensor
!###  otherwise the result is tensor.matrix
!###
   function Matrix_times_MPSTensor(tensorA, tensorB) result(this)
     type(MPSTensor),intent(IN) :: tensorA,tensorB
     type(MPSTensor) this
     integer :: s

     if(tensorA%Initialized.and.tensorB%Initialized) then
        if(tensorA%DRight.eq.tensorB%DLeft) then
           !The trick of using a tensor as a matrix is used here:
           if (tensorA%spin.eq.MatrixSpin) then
              this = new_MPSTensor(tensorB%spin,tensorA%DLeft,tensorB%DRight,zero)
              do s=1,tensorB%spin
                 this%data(:,:,s)=matmul(tensorA%data(:,:,MatrixSpin),tensorB%data(:,:,s)) !+this%data(:,:,s)
              enddo
           else if (tensorB%spin.eq.MatrixSpin) then
              this = new_MPSTensor(tensorA%spin,tensorA%DLeft,tensorB%DRight,zero)
              do s=1,tensorA%spin
                 this%data(:,:,s)=matmul(tensorA%data(:,:,s),tensorB%data(:,:,MatrixSpin)) !+this%data(:,:,s)
              enddo
           else
              call ThrowException('Matrix_times_MPSTensor','One of the arguments must be a *matrix* (MPSTensor with spin 1)' &
                   & ,NoErrorCode,CriticalError)
           endif
        else
           call ThrowException('Matrix_times_MPSTensor','Dimensions of the tensors do not match',NoErrorCode,CriticalError)
        endif
     else
        call ThrowException('Matrix_times_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
     endif
     return
   end function Matrix_times_MPSTensor


!##################################################################
!##################################################################
!##################################################################
!##################################################################

   real function Difference_btw_MPSTensors(tensor1, tensor2) result(diff)
     type(MPSTensor),intent(IN) :: tensor1,tensor2
     integer :: n,alpha,beta

     diff=0.0d0
     if(tensor1%Initialized.and.tensor2%Initialized) then
        if(tensor1.equaldims.tensor2) then
           do n=1,tensor1%spin
              do beta=1,tensor1%DRight
                 do alpha=1,tensor1%DLeft
                    diff=diff+abs(tensor1%data(alpha,beta,n)-tensor2%data(alpha,beta,n))
                 enddo
              enddo
           enddo
        else
           call ThrowException('Difference_btw_MPSTensors','Tensors of different shape',NoErrorCode,CriticalError)
        endif
        return
     else
        call ThrowException('Difference_btw_MPSTensors','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Difference_btw_MPSTensors


   real function Difference_btw_MPSTensors_WithAbsoluteValue(tensor1, tensor2) result(diff)
     type(MPSTensor),intent(IN) :: tensor1,tensor2
     integer :: n,alpha,beta

     diff=0.0d0
     if(tensor1%Initialized.and.tensor2%Initialized) then
        if(tensor1.equaldims.tensor2) then
           do n=1,tensor1%spin
              do beta=1,tensor1%DRight
                 do alpha=1,tensor1%DLeft
                    diff=diff+abs(abs(tensor1%data(alpha,beta,n))-abs(tensor2%data(alpha,beta,n)))
                 enddo
              enddo
           enddo
        else
           call ThrowException('Difference_btw_MPSTensors','Tensors of different shape',NoErrorCode,CriticalError)
        endif
        return
     else
        call ThrowException('Difference_btw_MPSTensors','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function Difference_btw_MPSTensors_WithAbsoluteValue


!##################################################################

   logical function MPSTensors_are_of_equal_Shape(tensor1,tensor2) result(equals)
     type(MPSTensor),intent(IN) :: tensor1,tensor2

     if(tensor1%Initialized.and.tensor2%Initialized) then
        equals=(tensor1%spin.eq.tensor2%spin).and.(tensor1%DLeft.eq.tensor2%DLeft).and.(tensor1%DRight.eq.tensor2%DRight)
        return
     else
        call ThrowException('MPSTensors_are_of_equal_Shape','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif

   end function MPSTensors_are_of_equal_Shape

!#######################################################################################
!#######################################################################################
! This are very important functions as most of the algorithm time is spent updating the
! matrices using this Left and Right products with tensors.
! They are also used heavily for computing expectation values, so optimization here might be key.
! Convention:
!              L(betaR,alphaR) = \sum_i B_i^\dagger . L_in A_i
!                              = \sum_i \sum_betaL \sum_alphaL  B^*_{i,betaR,betaL} Lin_{betaL,alphaL} A_{i,alphaL,alphaR)
!
!              R(alphaL,betaL) = \sum_i A_i . RL_in . B_i^\dagger
!                              = \sum_i \sum_alphaR \sum_betaR A_{i,alphaL,alphaR) Rin_{alphaR,betaR} B^*_{i,betaR,betaL}
!
!#######################################################################################
!#######################################################################################

   function MPS_Left_Product(TensorA,TensorB,matrixin) result(matrixout)
     type(MPSTensor),intent(IN) :: TensorA,TensorB
     type(MPSTensor) :: matrixout
     type(MPSTensor),intent(IN),optional :: matrixin
     type(MPSTensor) :: TempMatrix,L_in_matrix
     integer :: s,i,k,j,l
     complex(8) :: TEMP

     if((.not.tensorA%Initialized).and.(.not.tensorB%Initialized)) then
        call ThrowException('MPSLeftProduct','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif
     if (TensorA%spin.ne.TensorB%spin) then
        call ThrowException('MPSLeftProduct','Tensors have different spin',NoErrorCode,CriticalError)
        return
     endif
     if (present(matrixin)) then
        if(matrixin%Initialized) then
           L_in_matrix=new_MPSTensor(matrixin)
        else
           call ThrowException('MPSLeftProduct','Matrix is not initialized',NoErrorCode,CriticalError)
           return
        endif
     else
        L_in_matrix=new_MPSTensor(MatrixSpin,TensorB%DLeft,TensorA%DLeft,one)
     endif

     matrixout=new_MPSTensor(MatrixSpin, TensorB%DRight,TensorA%DRight, zero)
     TempMatrix=new_MPSTensor(MatrixSpin,TensorB%DLeft ,TensorA%DRight, zero)

     !The multiplications are done by hand because I could not get ZGEMM to work properly
     do s=1,TensorA%spin
        Tempmatrix%data(:,:,MatrixSpin)=matmul(L_in_matrix%data(:,:,MatrixSpin),TensorA%data(:,:,s))
        MatrixOut%data(:,:,MatrixSpin)=MatrixOut%data(:,:,MatrixSpin)+matmul(dconjg(transpose(TensorB%data(:,:,s))),Tempmatrix%data(:,:,MatrixSpin))
    enddo
    return
  end function MPS_Left_Product


!##################################################################
!##################################################################

  function MPS_Right_Product(TensorA,TensorB,matrixin) result(matrixout)
    type(MPSTensor),intent(IN) :: TensorA,TensorB
    type(MPSTensor) :: matrixout
    type(MPSTensor),intent(IN),optional :: matrixin
    type(MPSTensor) :: TempMatrix,R_in_matrix
    integer :: s,i,k,j,l
    complex(8) :: TEMP

    if((.not.tensorA%Initialized).and.(.not.tensorB%Initialized)) then
       call ThrowException('MPSRightProduct','Tensors not initialized',NoErrorCode,CriticalError)
       return
    endif
    if (TensorA%spin.ne.TensorB%spin) then
       call ThrowException('MPSRightProduct','Tensors have different spin',NoErrorCode,CriticalError)
       return
    endif

    matrixout=new_MPSTensor(MatrixSpin, TensorA%DLeft,TensorB%DLeft, zero)
    TempMatrix=new_MPSTensor(MatrixSpin,TensorA%DLeft ,TensorB%DRight, zero)

    if (present(matrixin)) then
       if(matrixin%Initialized) then
          R_in_matrix=new_MPSTensor(matrixin)
       else
          call ThrowException('MPSRightProduct','Matrix is not initialized',NoErrorCode,CriticalError)
          return
       endif
    else
       R_in_matrix=new_MPSTensor(MatrixSpin,TensorA%DRight,TensorB%DRight,one)
    endif

    !The multiplications are done by hand because I could not get ZGEMM to work properly
    do s=1,TensorA%spin
       Tempmatrix%data(:,:,MatrixSpin)=matmul(TensorA%data(:,:,s),R_in_matrix%data(:,:,MatrixSpin))
       MatrixOut%data(:,:,MatrixSpin)=MatrixOut%data(:,:,MatrixSpin)+matmul(Tempmatrix%data(:,:,MatrixSpin),transpose(dconjg(TensorB%data(:,:,s))))
    enddo
    return
  end function MPS_Right_Product



!#######################################################################################
!#######################################################################################
!
!                                    HELPER CODE
!
!#######################################################################################
!#######################################################################################

  subroutine CollapseSpinWithBond(this,collapsed,whichDimension)
    class(MPSTensor),intent(IN) :: this
    complex(8),intent(OUT) :: collapsed(:,:)
    integer,intent(IN) :: whichDimension
    integer :: s,alpha,beta,leftIndex,rightIndex,leftStep,rightStep,leftDimension,rightDimension

    if(.not.this%Initialized) then
       call ThrowException('CollapseSpinWithBond','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    if (whichDimension.eq.FirstDimension) then
       leftStep=this%DLeft
       rightStep=0
       leftDimension=(this%spin*this%DLeft)
       rightDimension=(this%DRight)
    else if (whichDimension.eq.SecondDimension) then
       leftStep=0
       rightStep=this%DRight
       leftDimension=(this%DLeft)
       rightDimension=(this%spin*this%DRight)
    else
       call ThrowException('CollapseSpinWithBond','Wrong Dimension parameter',whichDimension,CriticalError)
       return
    endif
    if ((size(collapsed,1).ne.leftDimension).and.(size(collapsed,2).ne.rightDimension)) then
       call ThrowException('CollapseSpinWithBond','Matrix for collapsed tensor does not have right dimensions' &
            & ,whichDimension,CriticalError)
       return
    endif

    !This always puts the spin before the bond dimension,
    !      [(s,alpha),(beta)]   or  [(alpha),(s,beta)]
    do s=1,this%spin
       do beta=1,this%DRight
          rightIndex=beta+(s-1)*rightStep
          do alpha=1,this%DLeft
             leftIndex=alpha+(s-1)*leftStep
             collapsed(leftIndex,rightIndex)=this%data(alpha,beta,s)
          enddo
       enddo
    enddo

  end subroutine CollapseSpinWithBond


end module Tensor_Class
