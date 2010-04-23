!!  This module describes an object with three legs as used by MPS algorithms.
!!  One leg is "special", it is the physical leg of the tensor.

module MPSTensor_Class

  use ErrorHandling
  use Constants

  implicit none

  !Set maximum bond dimensions and spin
  integer,parameter :: MAX_spin = 2, MAX_D = 100

  type MPSTensor
     private
     integer spin_,DLeft_,DRight_ !D is bond dimension
     logical :: initialized_=.false.
     complex(8) data_(MAX_D,MAX_D,MAX_spin) !Now uses fixed max dimensions and internal variables to adjust
   contains
     procedure delete => delete_MPSTensor
     procedure print => print_MPSTensor
     procedure DRight => DRight_MPSTensor
     procedure DLeft => DLeft_MPSTensor
     procedure Spin => Spin_MPSTensor
  end type MPSTensor

  interface new_MPSTensor
     module procedure new_MPSTensor_Random,new_MPSTensor_fromMPSTensor,new_MPSTensor_withData
  end interface

  interface operator (*)
     module procedure Integer_times_MPSTensor,Real_times_MPSTensor,Complex_times_MPSTensor,Real8_times_MPSTensor,Complex8_times_MPSTensor
  end interface

  interface operator (.diff.)
     module procedure Difference_btw_MPSTensors
  end interface

  interface operator (.equaldims.)
     module procedure  MPSTensors_are_of_equal_Shape
  end interface


 contains

   
!###########  new
   function new_MPSTensor_Random (spin,DLeft,DRight) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     type(MPSTensor) :: this
     integer :: n,alpha,beta
     integer,save :: iseed = 101

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
     do n=1,spin
        do beta=1,DRight
           do alpha=1,DLeft
              this%data_(alpha,beta,n)=ran(iseed)+II*ran(iseed)
           enddo
        enddo
     enddo
     !Flip flag
     this%initialized_=.true.

   end function new_MPSTensor_Random
!###########

   function new_MPSTensor_withData (spin,DLeft,DRight,originalData) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     complex(8),intent(in) :: originalData(:,:,:)
     type(MPSTensor) this

     if(spin.gt.MAX_spin.or.DLeft.gt.MAX_D.or.DRight.gt.MAX_D) then
        call ThrowException('new_MPSTensor_withData','spin or bond dimension larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(spin.lt.1.or.DLeft.lt.1.or.DRight.lt.1) then
        call ThrowException('new_MPSTensor_withData','spin or bond dimension smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     !initialize internal variables
     this%spin_=spin
     this%DLeft_=DLeft
     this%DRight_=DRight
     !initialize data
     this%data_=originalData
     !flip flag
     this%initialized_=.true.

   end function new_MPSTensor_withData


   function new_MPSTensor_fromMPSTensor (tensor) result (this)
     type(MPSTensor),intent(in) :: tensor
     type(MPSTensor) this

     if(.not.tensor%initialized_) then
        call ThrowException('new_MPSTensor_fromMPSTensor','Original tensor not initialized',NoErrorCode,CriticalError)
        return
     endif

     !initialize internal variables
     this%spin_=tensor%spin_
     this%DLeft_=tensor%DLeft_
     this%DRight_=tensor%DRight_
     !initialize data
     this%data_=tensor%data_
     !flip flag
     this%initialized_=.true.

   end function new_MPSTensor_fromMPSTensor


!########### delete
   integer function delete_MPSTensor (this) result(error)
     class(MPSTensor),intent(INOUT) :: this

     error=Warning

     if(.not.this%initialized_) then
        call ThrowException('delete_MPSTensor','Trying to delete an uninitialized tensor',NoErrorCode,Warning)
        error = Warning
        return
     endif
     
     !Erase info
     this%spin_=0
     this%DLeft_=0
     this%DRight_=0
     !Erase data
     this%data_=0.0d0
     !Flip flag
     this%initialized_=.false.     

     error=Normal

   end function delete_MPSTensor
!###########

!###########  print
   integer function Print_MPSTensor(this) result(error)
     class(MPSTensor),intent(IN) :: this
     integer i,j,k

     error = Warning

     if(.not.(this%initialized_)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,Warning)
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
!###########       


!###########       Accessor methods
   integer function Spin_MPSTensor(this) result(s)
     class(MPSTensor),intent(IN) :: this
 
    if(.not.(this%initialized_)) then
        call ThrowException('Spin','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        s=this%spin_
     endif

   end function Spin_MPSTensor

   integer function DLeft_MPSTensor(this) result(DL)
     class(MPSTensor),intent(IN) :: this

     if(.not.(this%initialized_)) then
        call ThrowException('DLeft','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DL=this%DLeft_
     endif

   end function DLeft_MPSTensor

   integer function DRight_MPSTensor(this) result(DR)

     class(MPSTensor),intent(IN) :: this

     if(.not.(this%initialized_)) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DR=this%DRight_
     endif

   end function DRight_MPSTensor

!############        Products by a constant
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

!############        Products by a constant
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

!############        Products by a constant
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

!############        Products by a constant
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

!############        Products by a constant
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

   real function Difference_btw_MPSTensors(tensor1, tensor2) result(diff)
     type(MPSTensor),intent(IN) :: tensor1,tensor2

     if(tensor1%initialized_.and.tensor2%initialized_) then
        if(tensor1.equaldims.tensor2) then
           diff = sum(abs(tensor1%data_-tensor2%data_))
        else
           call ThrowException('Difference_btw_MPSTensors','Tensors of different shape',NoErrorCode,CriticalError)
        endif
        return 
     else
        call ThrowException('Difference_btw_MPSTensors','Tensor not initialized',NoErrorCode,CriticalError)
        return
     endif     

   end function Difference_btw_MPSTensors

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

 end module MPSTensor_Class
