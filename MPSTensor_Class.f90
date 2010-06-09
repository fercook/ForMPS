!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010

!!  This module contains a class of objects with three legs as used by MPS algorithms.
!!  One leg is "special", it is the physical leg of the tensor.
!!  In notation it is nice to have the spin as the first index, but in the implementation 
!!  it is the third, as it makes it better for performance

!!TODO: check performance of allocatable data 
!!TODO: Reimplement extendind the type from a Matrix class that performs some low level functions

module MPSTensor_Class

  use ErrorHandling
  use Constants
  use Matrix_Helper

  implicit none

  integer,parameter :: MAX_spin = 2
  integer,parameter :: MAX_D = 100

!###############################
!#####  The class main object
!###############################
  type MPSTensor
     private
     integer spin,DLeft,DRight
     logical :: initialized_=.false.
     complex(8),allocatable :: data_(:,:,:) !!$TODO: Change to extend a matrix (maybe abstract) type
   contains
     procedure :: delete => delete_MPSTensor
     procedure :: print => print_MPSTensor
     procedure :: PrintDimensions => Print_MPSTensor_Dimensions
     procedure :: getDRight => DRight_MPSTensor
     procedure :: getDLeft => DLeft_MPSTensor
     procedure :: getSpin => Spin_MPSTensor
     procedure :: LCanonize => Left_Canonize_MPSTensor
     procedure :: RCanonize => Right_Canonize_MPSTensor 
     procedure :: isInitialized => InitializationCheck
     procedure :: CopyFrom => new_MPSTensor_fromAssignment
     procedure :: Norm => Norm_Of_MPSTensor
     procedure :: ApplyOperator => Apply_Operator_From_Matrix
!     generic :: assignment(=) => CopyFrom
  end type MPSTensor

!###############################
!#####  Operators and methods
!###############################
  interface new_MPSTensor
     module procedure new_MPSTensor_Random,new_MPSTensor_fromMPSTensor,new_MPSTensor_withData, &
          & new_MPSTensor_withConstant,new_MPSTensor_fromMatrix
  end interface

  interface operator (*)
     module procedure Integer_times_MPSTensor,Real_times_MPSTensor,Complex_times_MPSTensor,Real8_times_MPSTensor, &
          & Complex8_times_MPSTensor,Matrix_times_MPSTensor
  end interface

  interface assignment (=)
     module procedure new_MPSTensor_fromAssignment
  end interface

  interface operator (.diff.)
     module procedure Difference_btw_MPSTensors
  end interface

  interface operator (.absdiff.)
     module procedure Difference_btw_MPSTensors_WithAbsoluteValue
  end interface

  interface operator (.equaldims.)
     module procedure  MPSTensors_are_of_equal_Shape
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
     integer :: n,alpha,beta
     real(8) :: randomtensorR(DLeft,DRight,spin),randomtensorC(DLeft,DRight,spin)

     if(spin.gt.MAX_spin.or.DLeft.gt.MAX_D.or.DRight.gt.MAX_D) then
        call ThrowException('new_MPSTensor_Random','spin or bond dimension larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(spin.lt.1.or.DLeft.lt.1.or.DRight.lt.1) then
        call ThrowException('new_MPSTensor_Random','spin or bond dimension smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     !initialize internal variables
     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight
     !initialize data
     if(this%initialized_) deallocate(this%data_)
     allocate(this%data_(DLeft,DRight,spin))

     Call random_number(randomtensorR)
     call random_number(randomtensorC)

     This%data_=randomtensorR+II*randomtensorC

     this%initialized_=.true.

   end function new_MPSTensor_Random

!##################################################################
   function new_MPSTensor_withData (spin,DLeft,DRight,originalData) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     complex(8),intent(in) :: originalData(:,:,:)
     type(MPSTensor) this
     integer :: n,alpha,beta

     if(spin.gt.MAX_spin.or.DLeft.gt.MAX_D.or.DRight.gt.MAX_D) then
        call ThrowException('new_MPSTensor_withData','spin or bond dimension larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(spin.lt.1.or.DLeft.lt.1.or.DRight.lt.1) then
        call ThrowException('new_MPSTensor_withData','spin or bond dimension smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight

     if(this%initialized_) deallocate(this%data_)
     allocate(this%data_(DLeft,DRight,spin))
     this%data_=zero
     do n=1,spin
        do beta=1,DRight
           do alpha=1,DLeft
              this%data_(alpha,beta,n)=originalData(alpha,beta,n)
           enddo
        enddo
     enddo
     this%initialized_=.true.

   end function new_MPSTensor_withData

!##################################################################
   function new_MPSTensor_withConstant (spin,DLeft,DRight,constant) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     complex(8),intent(in) :: constant
     type(MPSTensor) this
     integer :: n,alpha,beta

     if(spin.gt.MAX_spin.or.DLeft.gt.MAX_D.or.DRight.gt.MAX_D) then
        call ThrowException('new_MPSTensor_withData','spin or bond dimension larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(spin.lt.1.or.DLeft.lt.1.or.DRight.lt.1) then
        call ThrowException('new_MPSTensor_withData','spin or bond dimension smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight

     if(this%initialized_) deallocate(this%data_)
     allocate(this%data_(DLeft,DRight,spin))
     this%data_=zero
     do n=1,spin
        do beta=1,DRight
           do alpha=1,DLeft
              this%data_(alpha,beta,n)=constant
           enddo
        enddo
     enddo
     this%initialized_=.true.

   end function new_MPSTensor_withConstant

!##################################################################
   function new_MPSTensor_fromMPSTensor (tensor) result (this)
     type(MPSTensor),intent(in) :: tensor
     type(MPSTensor) this
     integer error

     error=tensor%isInitialized()
     if (WasThereError()) call ProcessException('new_MPSTensor_fromMPSTensor')

     this%spin=tensor%spin
     this%DLeft=tensor%DLeft
     this%DRight=tensor%DRight
     if(this%initialized_) deallocate(this%data_)
     allocate(this%data_(this%DLeft,this%DRight,this%spin))
     this%data_=zero
     this%data_=tensor%data_
     this%initialized_=.true.

   end function new_MPSTensor_fromMPSTensor

   subroutine new_MPSTensor_fromAssignment(lhs,rhs)
     class(MPSTensor),intent(out) :: lhs
     type(MPSTensor),intent(in) :: rhs

     if(.not.rhs%initialized_) then
        call ThrowException('new_MPSTensor_fromAssignment','Original tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

     lhs%spin=rhs%spin
     lhs%DLeft=rhs%DLeft
     lhs%DRight=rhs%DRight

!    The following check should be
!          if(lhs%initialized_) deallocate(lhs%data_)
!    But instead I check for allocated because of a bug in GFortran,
!    http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43969
     if(allocated(lhs%data_)) deallocate(lhs%data_)

     allocate(lhs%data_(lhs%DLeft,lhs%DRight,lhs%spin))
     lhs%data_=zero
     lhs%data_=rhs%data_
     lhs%initialized_=.true.

   end subroutine new_MPSTensor_fromAssignment

  function new_MPSTensor_fromMatrix(spin,LeftBond,RightBond,matrix) result(this)
    complex(8),intent(IN) :: matrix(:,:)
    type(MPSTensor) :: this
    integer :: LeftBond,RightBond
    integer :: alpha,beta,s,spin
    integer :: leftIndex,rightIndex,leftStep,rightStep,leftDimension,rightDimension
    character(100) :: Message
    character(3) :: ScratchMessage

    if (size(matrix,1).eq.Spin*LeftBond) then
       leftStep=LeftBond
       rightStep=0
    else if (size(matrix,2).eq.Spin*RightBond) then
       leftStep=0
       rightStep=RightBond
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

    this%spin=spin
    this%DLeft=LeftBond
    this%DRight=RightBond
    if(this%initialized_) deallocate(this%data_)
    allocate(this%data_(this%DLeft,this%DRight,this%spin))
    this%data_=zero
    do s=1,Spin
       do beta=1,RightBond
          rightIndex=beta+(s-1)*rightStep
          do alpha=1,LeftBond
             leftIndex=alpha+(s-1)*leftStep
             this%data_(alpha,beta,s)=matrix(leftIndex,rightIndex)
          enddo
       enddo
    enddo
    
    this%initialized_=.true.
            
  end function new_MPSTensor_fromMatrix

!######################################    delete
   integer function delete_MPSTensor (this) result(error)
     !!class(MPSTensor),intent(INOUT) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(INOUT) :: this   !!<<TYPE>>!!

     error=Warning

     if(.not.this%initialized_) then
        call ThrowException('delete_MPSTensor','Trying to delete an uninitialized tensor',NoErrorCode,error)
        return
     endif
     
     !Erase info
     this%spin=0
     this%DLeft=0
     this%DRight=0
     !Erase data
     deallocate(this%data_)
     !Flip flag
     this%initialized_=.false.     

     error=Normal

   end function delete_MPSTensor
!##################################################################


integer function InitializationCheck(this) result(error)
    !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
    class(MPSTensor),intent(IN) :: this !!<<TYPE>>!!

    if (.not.this%initialized_) then    
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

     if(.not.(this%initialized_)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,error)
        return
     endif

     do i=1,this%spin
        print *,'State :',i
        do j=1,this%DLeft
           print *,(this%data_(j,k,i),k=1,this%DRight)
        enddo
     enddo

     error=Normal

   end function Print_MPSTensor



   integer function Print_MPSTensor_Dimensions(this) result(error)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     error = Warning

     if(.not.(this%initialized_)) then
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
 
    if(.not.(this%initialized_)) then
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

     if(.not.(this%initialized_)) then
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

     if(.not.(this%initialized_)) then
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
              Norm_Of_MPSTensor=Norm_Of_MPSTensor+abs(this%data_(alpha,beta,s))
           enddo
        enddo
     enddo

     return
   end function Norm_Of_MPSTensor


   function Integer_times_MPSTensor(constant, tensor) result(this)
     integer,intent(IN) :: constant
     type(MPSTensor),intent(IN) :: tensor
     type(MPSTensor) this

     if(tensor%initialized_) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data_)
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

     if(tensor%initialized_) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data_)
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

     if(tensor%initialized_) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data_)
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

     if(tensor%initialized_) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data_)
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

     if(tensor%initialized_) then
        this = new_MPSTensor(tensor%spin,tensor%DRight,tensor%DLeft,constant*tensor%data_)
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

     if(tensorA%initialized_.and.tensorB%initialized_) then
        if(tensorA%DRight.eq.tensorB%DLeft) then
           !The trick of using a tensor as a matrix is used here:
           if (tensorA%spin.eq.MatrixSpin) then
              this = new_MPSTensor(tensorB%spin,tensorA%DLeft,tensorB%DRight,zero)
              do s=1,tensorB%spin
                 this%data_(:,:,s)=matmul(tensorA%data_(:,:,MatrixSpin),tensorB%data_(:,:,s)) !+this%data_(:,:,s)
              enddo
           else if (tensorB%spin.eq.MatrixSpin) then
              this = new_MPSTensor(tensorA%spin,tensorA%DLeft,tensorB%DRight,zero)
              do s=1,tensorA%spin
                 this%data_(:,:,s)=matmul(tensorA%data_(:,:,s),tensorB%data_(:,:,MatrixSpin)) !+this%data_(:,:,s)
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



   function Apply_Operator_From_Matrix(this,matrix) result(aTensor)
     class(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     type(MPSTensor) :: aTensor
     complex(8),intent(IN) :: matrix(:,:)
     integer alpha,beta

     if(this%initialized_) then
        if(size(matrix,1).eq.this%spin.and.size(matrix,2).eq.this%spin) then
           aTensor = new_MPSTensor(this%spin,this%DLeft,this%DRight,zero)
           do beta=1,this%DRight
              do alpha=1,this%DLeft
                 aTensor%data_(alpha,beta,:)=matmul(matrix,this%data_(alpha,beta,:))
              enddo
           enddo
        else
           call ThrowException('Apply_Operator_From_Matrix','Operator is not of the rigt size',size(matrix,1)-this%spin,CriticalError)
        endif
     else
        call ThrowException('Apply_Operator_From_Matrix','Tensor not initialized',NoErrorCode,CriticalError)
     endif

   end function Apply_Operator_From_Matrix





!##################################################################
!##################################################################
!##################################################################
!##################################################################

   real function Difference_btw_MPSTensors(tensor1, tensor2) result(diff)
     type(MPSTensor),intent(IN) :: tensor1,tensor2
     integer :: n,alpha,beta
     
     diff=0.0d0
     if(tensor1%initialized_.and.tensor2%initialized_) then
        if(tensor1.equaldims.tensor2) then
           do n=1,tensor1%spin
              do beta=1,tensor1%DRight
                 do alpha=1,tensor1%DLeft
                    diff=diff+abs(tensor1%data_(alpha,beta,n)-tensor2%data_(alpha,beta,n))
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
     if(tensor1%initialized_.and.tensor2%initialized_) then
        if(tensor1.equaldims.tensor2) then
           do n=1,tensor1%spin
              do beta=1,tensor1%DRight
                 do alpha=1,tensor1%DLeft
                    diff=diff+abs(abs(tensor1%data_(alpha,beta,n))-abs(tensor2%data_(alpha,beta,n)))
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

     if(tensor1%initialized_.and.tensor2%initialized_) then
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
     
     if((.not.tensorA%initialized_).and.(.not.tensorB%initialized_)) then
        call ThrowException('MPSLeftProduct','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif     
     if (TensorA%spin.ne.TensorB%spin) then
        call ThrowException('MPSLeftProduct','Tensors have different spin',NoErrorCode,CriticalError)
        return
     endif
     if (present(matrixin)) then
        if(matrixin%initialized_) then
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
        Tempmatrix%data_(:,:,MatrixSpin)=matmul(L_in_matrix%data_(:,:,MatrixSpin),TensorA%data_(:,:,s))
        MatrixOut%data_(:,:,MatrixSpin)=MatrixOut%data_(:,:,MatrixSpin)+matmul(dconjg(transpose(TensorB%data_(:,:,s))),Tempmatrix%data_(:,:,MatrixSpin))
    enddo
    return 
  end function MPS_Left_Product

   function MPS_Left_ProductAlloc(TensorA,TensorB,matrixin) result(matrixout)
     type(MPSTensor),intent(IN) :: TensorA,TensorB
     complex(8),allocatable :: matrixout(:,:)
     complex(8),intent(IN),optional :: matrixin(:,:)
     complex(8),allocatable :: TempMatrix(:,:),L_in_matrix(:,:)
     integer :: s,i,k,j,l
     complex(8) :: TEMP
     
     if((.not.tensorA%initialized_).and.(.not.tensorB%initialized_)) then
        call ThrowException('MPSLeftProduct','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif     
     if (TensorA%spin.ne.TensorB%spin) then
        call ThrowException('MPSLeftProduct','Tensors have different spin',NoErrorCode,CriticalError)
        return
     endif
     if (present(matrixin)) then
        allocate(L_in_matrix(size(matrixin,1),size(matrixin,2)))
        L_in_matrix=matrixin
     else
        allocate(L_in_matrix(TensorB%DLeft,TensorA%DLeft))
        L_in_matrix=one
     endif

     allocate(matrixout(TensorB%DRight,TensorA%DRight))
     allocate(TempMatrix(TensorB%DLeft ,TensorA%DRight))
     matrixout=zero
     TempMatrix=zero

     !The multiplications are done by hand because I could not get ZGEMM to work properly
     do s=1,TensorA%spin
        Tempmatrix=matmul(L_in_matrix,TensorA%data_(:,:,s))
        MatrixOut=MatrixOut+matmul(dconjg(transpose(TensorB%data_(:,:,s))),Tempmatrix)
    enddo
    return 
  end function MPS_Left_ProductAlloc

!##################################################################
!##################################################################

  function MPS_Right_Product(TensorA,TensorB,matrixin) result(matrixout)
    type(MPSTensor),intent(IN) :: TensorA,TensorB
    type(MPSTensor) :: matrixout
    type(MPSTensor),intent(IN),optional :: matrixin
    type(MPSTensor) :: TempMatrix,R_in_matrix
    integer :: s,i,k,j,l
    complex(8) :: TEMP
    
    if((.not.tensorA%initialized_).and.(.not.tensorB%initialized_)) then
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
       if(matrixin%initialized_) then
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
       Tempmatrix%data_(:,:,MatrixSpin)=matmul(TensorA%data_(:,:,s),R_in_matrix%data_(:,:,MatrixSpin))
       MatrixOut%data_(:,:,MatrixSpin)=MatrixOut%data_(:,:,MatrixSpin)+matmul(Tempmatrix%data_(:,:,MatrixSpin),transpose(dconjg(TensorB%data_(:,:,s))))
    enddo
    return 
  end function MPS_Right_Product
  

!##################################################################
!##################################################################
! Left Site Canonization -- Returns the matrix that needs to be multiplied
! to the adjacent site on the RIGHT
!##################################################################
!##################################################################

  function Left_Canonize_MPSTensor(this) result(matrix)
    !!class(MPSTensor),intent(INOUT) :: this !!<<CLASS>>!!
    class(MPSTensor),intent(INOUT) :: this  !!<<TYPE>>!!
    type(MPSTensor) :: matrix
    complex(8), allocatable :: U(:,:),vTransposed(:,:),collapsedTensor(:,:)
    real(8),allocatable :: Sigma(:)
    integer :: Spin,LeftBond,RightBond
    integer :: newLeftBond,newRightBond
    integer :: jj,kk

    if(.not.this%initialized_) then
       call ThrowException('Left_Canonize_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    Spin=this%spin
    LeftBond=this%DLeft
    RightBond=this%DRight

    allocate(collapsedTensor(Spin*LeftBond,RightBond))
    allocate(U(Spin*LeftBond,Spin*LeftBond))
    allocate(Sigma(Min(Spin*LeftBond,RightBond)))
    allocate(vTransposed(RightBond,RightBond))

    call CollapseSpinWithBond(this,collapsedTensor,FirstDimension)
    if (WasThereError()) then
       call ThrowException('Left_Canonize_MPSTensor','Could not collapse the tensor',NoErrorCode,CriticalError)
       return
    endif

    if(Spin*LeftBond.gt.MAX_D) then
       call ThrowException('Left_Canonize_MPSTensor','Working dimension larger than Maximum',NoErrorCode,CriticalError)
       return
    endif

    kk= SingularValueDecomposition(CollapsedTensor,U,Sigma,vTransposed)

    newLeftBond=LeftBond
    newRightBond=Min(Spin*LeftBond,RightBond)

    this=new_MPSTensor(Spin,newLeftBond,newRightBond,U)
    if (WasThereError()) then
       call ThrowException('Left_Canonize_MPSTensor','Could not split the matrix',NoErrorCode,CriticalError)
       return
    endif

    !matrix is reshaped to fit the product with the tensor on the right
    matrix=new_MPSTensor(MatrixSpin,newRightBond,RightBond, & 
         & reshape(vecmat(Sigma,vTransposed) , [newRightBond,RightBond,MatrixSpin], &
         & Pad= [ (zero, kk=1,newRightBond*RightBond*MatrixSpin) ]   ) )  !! Pad with zeros at the end

  end function Left_Canonize_MPSTensor
  


! Right Site Canonization -- Returns the matrix that needs to be multiplied
! to the adjacent site on the LEFT
!##################################################################
!##################################################################

  function Right_Canonize_MPSTensor(this) result(matrix)
    !!class(MPSTensor),intent(INOUT) :: this !!<<CLASS>>!!
    class(MPSTensor),intent(INOUT) :: this  !!<<TYPE>>!!
    type(MPSTensor) :: matrix
    complex(8), allocatable :: U(:,:),vTransposed(:,:),collapsedTensor(:,:)
    real(8),allocatable :: Sigma(:)
    integer :: Spin,LeftBond,RightBond
    integer :: newLeftBond,newRightBond
    integer :: jj,kk

    if(.not.this%initialized_) then
       call ThrowException('Right_Canonize_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    Spin=this%spin
    LeftBond=this%DLeft
    RightBond=this%DRight

    allocate(collapsedTensor(LeftBond,Spin*RightBond))
    allocate(U(LeftBond,LeftBond))
    allocate(Sigma(Min(LeftBond,Spin*RightBond)))
    allocate(vTransposed(Spin*RightBond,Spin*RightBond))

    call CollapseSpinWithBond(this,collapsedTensor,SecondDimension)
    if (WasThereError()) then
       call ThrowException('Right_Canonize_MPSTensor','Could not collapse the tensor',NoErrorCode,CriticalError)
       return
    endif

    if(Spin*RightBond.gt.MAX_D) then
       call ThrowException('Left_Canonize_MPSTensor','Working dimension larger than Maximum',NoErrorCode,CriticalError)
       return
    endif

    kk= SingularValueDecomposition(CollapsedTensor,U,Sigma,vTransposed)

    newLeftBond=Min(Spin*RightBond,LeftBond)
    newRightBond=RightBond

    this=new_MPSTensor(Spin,newLeftBond,newRightBond,conjg(vTransposed))
    if (WasThereError()) then
       call ThrowException('Right_Canonize_MPSTensor','Could not split the matrix',NoErrorCode,CriticalError)
       return
    endif

    matrix=new_MPSTensor(MatrixSpin,LeftBond,newLeftBond, & 
         & reshape(matvec(U,Sigma) , [LeftBond,newLeftBond,MatrixSpin], &
         & Pad= [ (zero, kk=1,LeftBond*newLeftBond*MatrixSpin) ]   ) )  

  end function Right_Canonize_MPSTensor

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

    if(.not.this%initialized_) then
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
             collapsed(leftIndex,rightIndex)=this%data_(alpha,beta,s)
          enddo
       enddo
    enddo
             
  end subroutine CollapseSpinWithBond


 end module MPSTensor_Class
