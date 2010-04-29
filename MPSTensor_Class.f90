!!  This module describes an object with three legs as used by MPS algorithms.
!!  One leg is "special", it is the physical leg of the tensor.

module MPSTensor_Class

  use ErrorHandling
  use Constants

  implicit none

  !Set maximum bond dimensions and spin
  integer,parameter :: MAX_spin = 2, MAX_D = 100

!###############################
!#####  The class main object
!###############################
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
     procedure LCanonize => Left_Canonize_MPSTensor
     procedure RCanonize => Right_Canonize_MPSTensor
  end type MPSTensor

!###############################
!#####  Operators and methods
!###############################
  interface new_MPSTensor
     module procedure new_MPSTensor_Random,new_MPSTensor_fromMPSTensor,new_MPSTensor_withData,new_MPSTensor_withConstant
  end interface

  interface operator (*)
     module procedure Integer_times_MPSTensor,Real_times_MPSTensor,Complex_times_MPSTensor,Real8_times_MPSTensor, &
          & Complex8_times_MPSTensor,Matrix_times_MPSTensor
  end interface

  interface operator (.diff.)
     module procedure Difference_btw_MPSTensors
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
     this%data_=zero
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

     !initialize internal variables
     this%spin_=spin
     this%DLeft_=DLeft
     this%DRight_=DRight
     !initialize data
     this%data_=zero
     do n=1,spin
        do beta=1,DRight
           do alpha=1,DLeft
              this%data_(alpha,beta,n)=originalData(alpha,beta,n)
           enddo
        enddo
     enddo
     !flip flag
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

     !initialize internal variables
     this%spin_=spin
     this%DLeft_=DLeft
     this%DRight_=DRight
     !initialize data
     this%data_=zero
     do n=1,spin
        do beta=1,DRight
           do alpha=1,DLeft
              this%data_(alpha,beta,n)=constant
           enddo
        enddo
     enddo
     !flip flag
     this%initialized_=.true.

   end function new_MPSTensor_withConstant

!##################################################################
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
     this%data_=zero
     this%data_=tensor%data_
     !flip flag
     this%initialized_=.true.

   end function new_MPSTensor_fromMPSTensor


!######################################    delete
   integer function delete_MPSTensor (this) result(error)
     class(MPSTensor),intent(INOUT) :: this

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
     this%data_=zero
     !Flip flag
     this%initialized_=.false.     

     error=Normal

   end function delete_MPSTensor
!##################################################################

!######################################     print
   integer function Print_MPSTensor(this) result(error)
     class(MPSTensor),intent(IN) :: this
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
!##################################################################   

!##################################################################
!###########       Accessor methods
!##################################################################
   integer function Spin_MPSTensor(this) result(s)
     class(MPSTensor),intent(IN) :: this
 
    if(.not.(this%initialized_)) then
        call ThrowException('Spin','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        s=this%spin_
     endif

   end function Spin_MPSTensor
!##################################################################

   integer function DLeft_MPSTensor(this) result(DL)
     class(MPSTensor),intent(IN) :: this

     if(.not.(this%initialized_)) then
        call ThrowException('DLeft','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DL=this%DLeft_
     endif

   end function DLeft_MPSTensor
!##################################################################
   integer function DRight_MPSTensor(this) result(DR)

     class(MPSTensor),intent(IN) :: this

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
!###  One of the arguments must be a "matrix", i.e. a MPSTensor with spin=1
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
           !Recognize the matrix because it has spin_=1
           if (tensorA%spin_.eq.1) then
              this = new_MPSTensor(tensorB%spin_,tensorA%DLeft_,tensorB%DRight_,zero)
              do s=1,tensorB%spin_
                 call mymatmul(tensorA%data_(:,:,1),tensorB%data_(:,:,s),this%data_(:,:,s), &
                      & tensorA%DLeft_,tensorB%DLeft_,tensorB%DRight_,'N')
              enddo
           else if (tensorB%spin_.eq.1) then
              this = new_MPSTensor(tensorA%spin_,tensorA%DLeft_,tensorB%DRight_,zero)
              do s=1,tensorA%spin_
                 call mymatmul(tensorA%data_(:,:,s),tensorB%data_(:,:,1),this%data_(:,:,s), &
                      & tensorA%DLeft_,tensorB%DLeft_,tensorB%DRight_,'N')
              enddo              
           else
              call ThrowException('Matrix_times_MPSTensor','One of the arguments must be a *matrix* (MPSTensor with spin 1)',NoErrorCode,CriticalError)
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
        L_in_matrix=new_MPSTensor(1,TensorB%DLeft_,TensorA%DLeft_,one)
     endif

     matrixout=new_MPSTensor(1, TensorB%DRight_,TensorA%DRight_, zero)
     TempMatrix=new_MPSTensor(1,TensorB%DLeft_ ,TensorA%DRight_, zero)

     !The multiplications are done by hand because I could not get ZGEMM to work properly
     do s=1,TensorA%Spin_
        !TempMatrix = InternalMatrix * TensorA
        Tempmatrix%data_=0.0d0        
        call mymatmul(L_in_matrix%data_(:,:,1),TensorA%data_(:,:,s),Tempmatrix%data_(:,:,1), &
             & TensorB%DLeft_,TensorA%DLeft_,TensorA%DRight_,'N')
        ! matrixout = TensorB^\dagger * TempMatrix
        call mymatmul(TensorB%data_(:,:,s),Tempmatrix%data_(:,:,1), MatrixOut%data_(:,:,1), &
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
    
    matrixout=new_MPSTensor(1, TensorA%DLeft_,TensorB%DLeft_, zero)
    TempMatrix=new_MPSTensor(1,TensorA%DLeft_ ,TensorB%DRight_, zero)
    
    if (present(matrixin)) then
       if(matrixin%initialized_) then
          R_in_matrix=new_MPSTensor(matrixin)
       else
          call ThrowException('MPSRightProduct','Matrix is not initialized',NoErrorCode,CriticalError)
          return           
       endif
    else
       R_in_matrix=new_MPSTensor(1,TensorA%DRight_,TensorB%DRight_,one)
    endif
    
    !The multiplications are done by hand because I could not get ZGEMM to work properly
    do s=1,TensorA%Spin_
       Tempmatrix%data_=0.0d0  
       !TempMatrix =  TensorA * InternalMatrix
       call mymatmul(TensorA%data_(:,:,s),R_in_matrix%data_(:,:,1),Tempmatrix%data_(:,:,1), &
            & TensorA%DLeft_,TensorA%DRight_,TensorB%DRight_,'N')
       ! matrixout =  TempMatrix * TensorB^\dagger
       call mymatmul(Tempmatrix%data_(:,:,1),TensorB%data_(:,:,s),MatrixOut%data_(:,:,1), &
            & TensorA%DLeft_,TensorB%DRight_,TensorB%DLeft_,'B')
    enddo
    return 
  end function MPS_Right_Product
  

!##################################################################
!##################################################################
! Site Canonization -- Returns the matrix that needs to be multiplied
! to the adjacent site
!##################################################################
!##################################################################

  function Left_Canonize_MPSTensor(this) result(matrix)
    class(MPSTensor),intent(INOUT) :: this
    type(MPSTensor) :: matrix
    complex(8), allocatable :: u(:,:),v(:,:),w(:,:),collapsedTensor(:,:)
    integer :: dims

    if(.not.this%initialized_) then
       call ThrowException('Left_Canonize_MPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
       return
    endif

    allocate(collapsedTensor((this%spin_)*(this%DLeft_),this%DRight_))
    allocate(u((this%spin_)*(this%DLeft_),(this%spin_)*(this%DLeft_)))
    allocate(v((this%spin_)*(this%DLeft_),this%DRight_))
    allocate(w(this%DRight_,this%DRight_))

    

  end function Left_Canonize_MPSTensor
  

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


   subroutine mymatmul(A,B,C,indexL,indexC,indexR,mode)
     complex(8) :: A(:,:),B(:,:),C(:,:)
     integer :: indexL,indexC,indexR
     character*1 :: mode
     integer :: I,J,K,L
     complex(8) TEMP
     ! mode = 'N' is normal multiplication C = A * B + C
     ! mode = 'A' is with A dagged, C = A^+ * B + C
     ! mode = 'B' is with B dagged, C = A * B^+ + C
     if (mode.eq.'N'.or.mode.eq.'n') then
        !C = A * B + C
        DO J = 1,indexR
           DO L = 1,indexC
              IF (B(L,J).NE.ZERO) THEN
                 TEMP = B(L,J)
                 DO I = 1,indexL
                    C(I,J) = C(I,J) + A(I,L)*TEMP
                 enddo
              END IF
           enddo
        enddo
     else if (mode.eq.'A'.or.mode.eq.'a') then
        ! C = A^+ * B + C
        DO J = 1, indexR
           DO I = 1,indexL
              TEMP = ZERO
              DO L = 1,indexC
                 TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
              enddo
             C(I,J) = TEMP + C(I,J)
           enddo
        enddo
     else if (mode.eq.'B'.or.mode.eq.'b') then
        ! C = A * B^+ + C
        DO J = 1,indexR
           DO L = 1,indexC
              IF (B(J,L).NE.ZERO) THEN
                 TEMP = DCONJG(B(J,L))
                 DO I = 1,indexL
                    C(I,J) = C(I,J) + A(I,L)*TEMP
                 enddo
              END IF
           enddo
        enddo

     endif
   end subroutine mymatmul




 end module MPSTensor_Class
