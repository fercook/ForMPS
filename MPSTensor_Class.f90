!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010

!!  This module contains a class of objects with three legs as used by MPS algorithms.
!!  One leg is "special", it is the physical leg of the tensor.
!!  In notation it is nice to have the spin as the first index, but in the implementation 
!!  it is the third, as it makes it better for performance

!!TODO: Convert the fixed dimension of the matrices in the object to allocatable, check performance
!!TODO: Reimplement extendind the type from a Matrix class that performs some low level functions

module MPSTensor_Class

  use ErrorHandling
  use Constants
  use Matrix_Helper

  implicit none

  integer,parameter :: MAX_spin = 2, MAX_D = 100

!###############################
!#####  The class main object
!###############################
  type MPSTensor
     private
     integer spin_,DLeft_,DRight_ 
     logical :: initialized_=.false.
     complex(8),allocatable :: data_(:,:,:) !!$TODO: Change to extend a matrix (maybe abstract) type
   contains
     procedure :: delete => delete_MPSTensor
     procedure :: print => print_MPSTensor
     procedure :: PrintDimensions => Print_MPSTensor_Dimensions
     procedure :: DRight => DRight_MPSTensor
     procedure :: DLeft => DLeft_MPSTensor
     procedure :: Spin => Spin_MPSTensor
     procedure :: LCanonize => Left_Canonize_MPSTensor
     procedure :: RCanonize => Right_Canonize_MPSTensor 
     procedure :: isInitialized => InitializationCheck
     procedure :: CopyFrom => new_MPSTensor_fromAssignment
     procedure :: Norm => Norm_Of_MPSTensor
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
     this%spin_=spin
     this%DLeft_=DLeft
     this%DRight_=DRight
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

     this%spin_=spin
     this%DLeft_=DLeft
     this%DRight_=DRight

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

     this%spin_=spin
     this%DLeft_=DLeft
     this%DRight_=DRight

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

     this%spin_=tensor%spin_
     this%DLeft_=tensor%DLeft_
     this%DRight_=tensor%DRight_
     if(this%initialized_) deallocate(this%data_)
     allocate(this%data_(this%DLeft_,this%DRight_,this%spin_))
     this%data_=zero
     this%data_=tensor%data_
     this%initialized_=.true.

   end function new_MPSTensor_fromMPSTensor

   subroutine new_MPSTensor_fromAssignment(lhs,rhs)
     TYPEORCLASS(MPSTensor),intent(out) :: lhs
     type(MPSTensor),intent(in) :: rhs

     if(.not.rhs%initialized_) then
        call ThrowException('new_MPSTensor_fromAssignment','Original tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

     lhs%spin_=rhs%spin_
     lhs%DLeft_=rhs%DLeft_
     lhs%DRight_=rhs%DRight_

!    The following check should be
!          if(lhs%initialized_) deallocate(lhs%data_)
!    But instead I check for allocated because of a bug in GFortran,
!    http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43969
     if(allocated(lhs%data_)) deallocate(lhs%data_)

     allocate(lhs%data_(lhs%DLeft_,lhs%DRight_,lhs%spin_))
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

    this%spin_=spin
    this%DLeft_=LeftBond
    this%DRight_=RightBond
    if(this%initialized_) deallocate(this%data_)
    allocate(this%data_(this%DLeft_,this%DRight_,this%spin_))
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
     TYPEORCLASS(MPSTensor),intent(INOUT) :: this   !!<<TYPE>>!!

     error=Warning

     if(.not.this%initialized_) then
        call ThrowException('delete_MPSTensor','Trying to delete an uninitialized tensor',NoErrorCode,error)
        return
     endif
     
     !Erase info
     this%spin_=0
     this%DLeft_=0
     this%DRight_=0
     !Erase data
     deallocate(this%data_)
     !Flip flag
     this%initialized_=.false.     

     error=Normal

   end function delete_MPSTensor
!##################################################################


integer function InitializationCheck(this) result(error)
    !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
    TYPEORCLASS(MPSTensor),intent(IN) :: this !!<<TYPE>>!!

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
     TYPEORCLASS(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     error = Warning

     if(.not.(this%initialized_)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,error)
        return
     endif

     do i=1,this%spin_
        print *,'State :',i
        do j=1,this%DLeft_
           print *,(this%data_(j,k,i),k=1,this%DRight_)
        enddo
     enddo

     error=Normal

   end function Print_MPSTensor



   integer function Print_MPSTensor_Dimensions(this) result(error)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     TYPEORCLASS(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
     integer i,j,k

     error = Warning

     if(.not.(this%initialized_)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,error)
        return
     endif

     print *,'Spin = ',this%spin_
     print *,'DL = ',this%DLeft_
     print *,'DR = ',this%DRight_

     error=Normal

   end function Print_MPSTensor_Dimensions
!##################################################################   

!##################################################################
!###########       Accessor methods
!##################################################################
   integer function Spin_MPSTensor(this) result(s)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     TYPEORCLASS(MPSTensor),intent(IN) :: this  !!<<TYPE>>!!
 
    if(.not.(this%initialized_)) then
        call ThrowException('Spin','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        s=this%spin_
     endif

   end function Spin_MPSTensor
!##################################################################

   integer function DLeft_MPSTensor(this) result(DL)
     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     TYPEORCLASS(MPSTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%initialized_)) then
        call ThrowException('DLeft','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DL=this%DLeft_
     endif

   end function DLeft_MPSTensor
!##################################################################
   integer function DRight_MPSTensor(this) result(DR)

     !!class(MPSTensor),intent(IN) :: this !!<<CLASS>>!!
     TYPEORCLASS(MPSTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%initialized_)) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DR=this%DRight_
     endif

   end function DRight_MPSTensor

!##################################################################
!#######################        Products by things
!##################################################################

   real(8) function Norm_Of_MPSTensor(this)
     TYPEORCLASS(MPSTensor),intent(IN) :: this   !!<<TYPE>>!!
     integer :: s,alpha,beta

     Norm_Of_MPSTensor=0.0d0
     do s=1,this%spin_
        do beta=1,this%DRight_
           do alpha=1,this%DLeft_
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
        this = new_MPSTensor(tensor%spin_,tensor%DRight_,tensor%DLeft_,constant*tensor%data_)
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
        this = new_MPSTensor(tensor%spin_,tensor%DRight_,tensor%DLeft_,constant*tensor%data_)
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
        this = new_MPSTensor(tensor%spin_,tensor%DRight_,tensor%DLeft_,constant*tensor%data_)
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
        this = new_MPSTensor(tensor%spin_,tensor%DRight_,tensor%DLeft_,constant*tensor%data_)
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
        this = new_MPSTensor(tensor%spin_,tensor%DRight_,tensor%DLeft_,constant*tensor%data_)
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
        if(tensorA%DRight_.eq.tensorB%DLeft_) then
           !The trick of using a tensor as a matrix is used here:
           if (tensorA%spin_.eq.MatrixSpin) then
              this = new_MPSTensor(tensorB%spin_,tensorA%DLeft_,tensorB%DRight_,zero)
              do s=1,tensorB%spin_
                 call mymatmul(tensorA%data_(:,:,MatrixSpin),tensorB%data_(:,:,s),this%data_(:,:,s), &
                      & tensorA%DLeft_,tensorB%DLeft_,tensorB%DRight_,'N')
              enddo
           else if (tensorB%spin_.eq.MatrixSpin) then
              this = new_MPSTensor(tensorA%spin_,tensorA%DLeft_,tensorB%DRight_,zero)
              do s=1,tensorA%spin_
                 call mymatmul(tensorA%data_(:,:,s),tensorB%data_(:,:,MatrixSpin),this%data_(:,:,s), &
                      & tensorA%DLeft_,tensorB%DLeft_,tensorB%DRight_,'N')
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
     if(tensor1%initialized_.and.tensor2%initialized_) then
        if(tensor1.equaldims.tensor2) then
           do n=1,tensor1%spin_
              do beta=1,tensor1%DRight_
                 do alpha=1,tensor1%DLeft_
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
           do n=1,tensor1%spin_
              do beta=1,tensor1%DRight_
                 do alpha=1,tensor1%DLeft_
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
        equals=(tensor1%spin_.eq.tensor2%spin_).and.(tensor1%DLeft_.eq.tensor2%DLeft_).and.(tensor1%DRight_.eq.tensor2%DRight_)
        return 
     else
        call ThrowException('MPSTensors_are_of_equal_Shape','Tensors not initialized',NoErrorCode,CriticalError)
        return
     endif     

   end function MPSTensors_are_of_equal_Shape

!#######################################################################################
!#######################################################################################
! This are very important functions as most of the algorithm time is spent updating the
! matrices using this Left and Right products with tensors
! They are also used heavily for computing expectation values, so optimization here might be key
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
     if (TensorA%Spin_.ne.TensorB%spin_) then
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
        L_in_matrix=new_MPSTensor(MatrixSpin,TensorB%DLeft_,TensorA%DLeft_,one)
     endif

     matrixout=new_MPSTensor(MatrixSpin, TensorB%DRight_,TensorA%DRight_, zero)
     TempMatrix=new_MPSTensor(MatrixSpin,TensorB%DLeft_ ,TensorA%DRight_, zero)

     !The multiplications are done by hand because I could not get ZGEMM to work properly
     do s=1,TensorA%Spin_
        Tempmatrix%data_=0.0d0        
        call mymatmul(L_in_matrix%data_(:,:,MatrixSpin),TensorA%data_(:,:,s),Tempmatrix%data_(:,:,MatrixSpin), &
             & TensorB%DLeft_,TensorA%DLeft_,TensorA%DRight_,'N')
        call mymatmul(TensorB%data_(:,:,s),Tempmatrix%data_(:,:,MatrixSpin), MatrixOut%data_(:,:,MatrixSpin), &
             & TensorB%DRight_,TensorB%DLeft_,TensorA%DRight_,'A')
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
    
    if((.not.tensorA%initialized_).and.(.not.tensorB%initialized_)) then
       call ThrowException('MPSRightProduct','Tensors not initialized',NoErrorCode,CriticalError)
       return
    endif
    if (TensorA%Spin_.ne.TensorB%spin_) then
       call ThrowException('MPSRightProduct','Tensors have different spin',NoErrorCode,CriticalError)
       return
    endif
    
    matrixout=new_MPSTensor(MatrixSpin, TensorA%DLeft_,TensorB%DLeft_, zero)
    TempMatrix=new_MPSTensor(MatrixSpin,TensorA%DLeft_ ,TensorB%DRight_, zero)
    
    if (present(matrixin)) then
       if(matrixin%initialized_) then
          R_in_matrix=new_MPSTensor(matrixin)
       else
          call ThrowException('MPSRightProduct','Matrix is not initialized',NoErrorCode,CriticalError)
          return           
       endif
    else
       R_in_matrix=new_MPSTensor(MatrixSpin,TensorA%DRight_,TensorB%DRight_,one)
    endif
    
    !The multiplications are done by hand because I could not get ZGEMM to work properly
    do s=1,TensorA%Spin_
       Tempmatrix%data_=0.0d0  
       call mymatmul(TensorA%data_(:,:,s),R_in_matrix%data_(:,:,MatrixSpin),Tempmatrix%data_(:,:,MatrixSpin), &
            & TensorA%DLeft_,TensorA%DRight_,TensorB%DRight_,'N')
       call mymatmul(Tempmatrix%data_(:,:,MatrixSpin),TensorB%data_(:,:,s),MatrixOut%data_(:,:,MatrixSpin), &
            & TensorA%DLeft_,TensorB%DRight_,TensorB%DLeft_,'B')
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
    TYPEORCLASS(MPSTensor),intent(INOUT) :: this  !!<<TYPE>>!!
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

    Spin=this%spin_
    LeftBond=this%DLeft_
    RightBond=this%DRight_

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
    TYPEORCLASS(MPSTensor),intent(INOUT) :: this  !!<<TYPE>>!!
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

    Spin=this%spin_
    LeftBond=this%DLeft_
    RightBond=this%DRight_

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




!!$ Mathematica code for canonization
!!$
!!$ Options[MPSCanonizeSite] = {Direction -> "Right", UseMatrix -> True};
!!$ SetAttributes[MPSCanonizeSite, HoldAll];
!!$ MPSCanonizeSite[tensor_, matrix_, OptionsPattern[]] := 
!!$  Module[{sense = OptionValue[Direction], 
!!$    usematrix = OptionValue[UseMatrix], 
!!$    numTensors, \[Chi]L, \[Chi]R, \[Chi], u, v, t, newTensor},(* 
!!$   Start by multiplying the tensor with the matrix from the previous site *)
!!$   If[sense == "Right",
!!$    If[usematrix, newTensor = tensor.matrix, newTensor = tensor];
!!$    {\[Chi]L, \[Chi]R} = Dimensions[newTensor[[1]]];
!!$    \[Chi] = Max[\[Chi]L, \[Chi]R];
!!$    (* SVD of the new tensor, putting [chiL, spin*
!!$    chiR] *)

!!$    {u, v, t} = 
!!$     SingularValueDecomposition[Flatten[newTensor, {{2}, {1, 3}}]];
!!$    (* Prepare new right matrix *)
!!$    
!!$    matrix = 
!!$     PadRight[
!!$      u.v, {Min[\[Chi], \[Chi]L], Min[\[Chi], Length[t], \[Chi]L]}];
!!$    (* Form the new tensor with the first row of t^
!!$    dagger *)
!!$    (Partition[
!!$       ConjugateTranspose[
!!$        t], {Min[\[Chi], Length[t], \[Chi]L], \[Chi]R}][[1, All]])
!!$    , (* LEFT CANONIZATION *)
!!$    If[usematrix, newTensor = matrix.# & /@ tensor, 
!!$     newTensor = tensor];
!!$    {\[Chi]L, \[Chi]R} = Dimensions[newTensor[[1]]];
!!$    \[Chi] = Max[\[Chi]L, \[Chi]R];
!!$    (* SVD of the new tensor, putting [chiL*spin, 
!!$    chiR] *)
!!$    {u, v, t} = 
!!$     SingularValueDecomposition[Flatten[newTensor, {{1, 2}, {3}}]];
!!$    (* Prepare new right matrix *)
!!$    
!!$    matrix = 
!!$     PadRight[
!!$      v.ConjugateTranspose[t], {Min[\[Chi], Length[u], \[Chi]R], 
!!$       Min[\[Chi], \[Chi]R]}];
!!$    (* Form the new tensor with the first column of u *)
!!$    \
!!$    (Partition[u, {\[Chi]L, Min[\[Chi], Length[u], \[Chi]R]}][[All, 1]])
!!$    ]
!!$   ];



!#######################################################################################
!#######################################################################################
! 
!                                    HELPER CODE
!
!#######################################################################################
!#######################################################################################

  subroutine CollapseSpinWithBond(this,collapsed,whichDimension)
    TYPEORCLASS(MPSTensor),intent(IN) :: this
    complex(8),intent(OUT) :: collapsed(:,:)
    integer,intent(IN) :: whichDimension
    integer :: s,alpha,beta,leftIndex,rightIndex,leftStep,rightStep,leftDimension,rightDimension

    if(.not.this%initialized_) then
       call ThrowException('CollapseSpinWithBond','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    if (whichDimension.eq.FirstDimension) then
       leftStep=this%DLeft_
       rightStep=0
       leftDimension=(this%spin_*this%DLeft_)
       rightDimension=(this%DRight_)
    else if (whichDimension.eq.SecondDimension) then
       leftStep=0
       rightStep=this%DRight_
       leftDimension=(this%DLeft_)
       rightDimension=(this%spin_*this%DRight_)
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
    do s=1,this%Spin_
       do beta=1,this%DRight_
          rightIndex=beta+(s-1)*rightStep
          do alpha=1,this%DLeft_
             leftIndex=alpha+(s-1)*leftStep
             collapsed(leftIndex,rightIndex)=this%data_(alpha,beta,s)
          enddo
       enddo
    enddo
             
  end subroutine CollapseSpinWithBond


 end module MPSTensor_Class
