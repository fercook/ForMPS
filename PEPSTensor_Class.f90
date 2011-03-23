!!   Copyright 2010 Fernando M. Cucchietti
!
!    This file is part of ForMPS
!
!    ForMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ForMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ForMPS.  If not, see <http://www.gnu.org/licenses/>.

module PEPSTensor_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use MPSTensor_Class

  implicit none

!  private

!###############################
!#####  The class main object
!###############################
  type,public,extends(Tensor5) :: PEPSTensor
     private
   contains
     procedure,public :: getDRight => DRight_PEPSTensor
     procedure,public :: getDLeft => DLeft_PEPSTensor
     procedure,public :: getDUp => DUp_PEPSTensor
     procedure,public :: getDDown => DDown_PEPSTensor
     procedure,public :: getMaxBondDimension => Get_MaxBondDimensionPEPSTensor
     procedure,public :: getSpin => Spin_PEPSTensor
     procedure,public :: getBonds => Get_PEPSBondList
     procedure,public :: PrintDimensions => Print_PEPSTensor_Dimensions
     procedure,public :: ApplyOperator => Apply_Operator_To_PEPS_Spin_Dimension
     procedure,public :: ApplyMatrixToBond => Apply_Matrix_To_PEPS
     procedure,public :: ApplyMPOToBond => Apply_Tensor4_To_PEPS
     procedure,public :: Collapse => Collapse_PEPS_Into_Tensor4
     procedure,public :: CompactBonds => CompactPEPSBondDimensions
     procedure,public :: CollapseAllIndicesBut => CollapsePEPSandSpinIndicesBut
     procedure,public :: asMPSTensor => DropTwoBondsAndReturnMPSTensor
     procedure,public :: JoinSpinWith => JoinPEPSSpinWithBond
     procedure,public :: HOSVD => HighOrderSVDofPEPS
  end type PEPSTensor

!###############################
!#####  Operators and methods
!###############################
  interface new_PEPSTensor
     module procedure new_PEPSTensor_Random,new_PEPSTensor_fromPEPSTensor, &
          & new_PEPSTensor_withConstant, new_PEPSTensor_with_SplitData, &
          & new_PEPSTensor_fromTensor5, new_PEPSTensor_fromTensor5_Transposed,new_PEPSTensor_fromTensor3
  end interface

  interface assignment (=)
     module procedure new_PEPSTensor_fromAssignment
  end interface

  interface operator (.apply.)
     module procedure Apply_Operator_To_PEPS_Spin_Dimension
  end interface

  interface ApplyMatrixToBond
     module procedure Apply_Matrix_To_PEPS
  end interface

  interface ApplyMPOToBond
     module procedure Apply_Tensor4_To_PEPS
  end interface

  interface CollapsePEPS
     module procedure Collapse_PEPS_Into_Tensor4, Collapse_Two_PEPS_into_Tensor4
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
   function new_PEPSTensor_Random (spin,DLeft,DRight,Dup,DDown) result (this)
     integer,intent(in) :: spin,DLeft,DRight,Dup,DDown
     type(PEPSTensor) :: this

     this=new_Tensor(DLeft,DRight,Dup,DDown,spin)

   end function new_PEPSTensor_Random

!##################################################################
   function new_PEPSTensor_with_SplitData (originalDataUP,originalDataDOWN) result (this)
     complex(8),intent(in) :: originalDataUP(:,:,:,:),originalDataDOWN(:,:,:,:)
     complex(8),allocatable :: JoinedData(:,:,:,:,:)
     type(PEPSTensor) this
     integer :: upDims(4),downDims(4)

     upDims=shape(originalDataUP)
     downDims=shape(originalDataDOWN)

     if (upDims.equalvector.downDims) then
        allocate (JoinedData(upDims(1),upDims(2),upDims(3),upDims(4),2))
        JoinedData(:,:,:,:,1)=originalDataUP
        JoinedData(:,:,:,:,2)=originalDataDOWN

        this=new_Tensor(JoinedData)

      else
        call ThrowException('new_PEPSTensor_with_SplitData','Up and down data have different size',upDims(1)-downDims(2),CriticalError)
      endif
   end function new_PEPSTensor_with_SplitData

!##################################################################
   function new_PEPSTensor_withConstant (spin,DLeft,DRight,DUp,DDown,constant) result (this)
     integer,intent(in) :: spin,DLeft,DRight,DUp,DDown
     complex(8),intent(in) :: constant
     type(PEPSTensor) :: this

     this=new_Tensor(DLeft,DRight,DUp,DDown,spin,constant)

   end function new_PEPSTensor_withConstant

!##################################################################
   function new_PEPSTensor_fromPEPSTensor (tensor) result (this)
     class(PEPSTensor),intent(in) :: tensor
     type(PEPSTensor) this

     this=new_Tensor(tensor)

   end function new_PEPSTensor_fromPEPSTensor

!##################################################################
   subroutine new_PEPSTensor_fromAssignment(lhs,rhs)
     class(PEPSTensor),intent(out) :: lhs
     type(PEPSTensor),intent(in) :: rhs

     lhs=new_Tensor(rhs)

   end subroutine new_PEPSTensor_fromAssignment

!##################################################################
   function new_PEPSTensor_fromTensor5 (tensor) result (this)
     type(Tensor5),intent(in) :: tensor
     type(PEPSTensor) this

     this=new_Tensor(tensor)

   end function new_PEPSTensor_fromTensor5

   function new_PEPSTensor_fromTensor5_Transposed(tensor, &
       &  whichDimIsSpin,whichDimIsLeft,whichDimIsRight,whichDimIsUp,whichDimIsDown ) result (this)
     type(Tensor5),intent(in) :: tensor
     integer, intent(IN) :: whichDimIsSpin,whichDimIsLeft,whichDimIsRight,whichDimIsUp,whichDimIsDown
     type(PEPSTensor) this
     integer :: newDims(5)

     newDims(whichDimIsLeft)=1
     newDims(whichDimIsRight)=2
     newDims(whichDimIsUp)=3
     newDims(whichDimIsDown)=4
     newDims(whichDimIsSpin)=5
     this=TensorTranspose(tensor,newDims)

   end function new_PEPSTensor_fromTensor5_Transposed

!##################################################################
   function new_PEPSTensor_fromTensor3(tensor,whichDimIsSpin,DirectionOfBonds,transverseBondsSize) result (this)
      class(Tensor3),intent(in) :: tensor
      integer, intent(IN) :: whichDimIsSpin,DirectionOfBonds
      integer, intent(IN),optional :: transverseBondsSize(2)
      type(PEPSTensor) :: this
      integer :: newDims(5),oldDims(3),reorderedDims(3),reorderWithTransverse(5)

      if (present(transverseBondsSize).and.whichDimIsSpin.ne.3) then
         !I have not checked what happens when spin is not third dimension and requesting decomposition
         call ThrowException('new_PEPSTensor_fromTensor3','When decomposing spin into transverse bonds spin should be last dimension',NoErrorCode,Warning)
      endif

      oldDims=tensor%GetDimensions()

      if (present(transverseBondsSize)) then
        newDims(1)=oldDims(1)
        newDims(2)=oldDims(2)
        newDims(3)=transverseBondsSize(1)
        newDims(4)=transverseBondsSize(2)
        newDims(5)=oldDims(3) / product(transverseBondsSize) !spin
        if (DirectionOfBonds.eq.HORIZONTAL) then
            reorderWithTransverse=[1,2,3,4,5]
        else if (DirectionOfBonds.eq.VERTICAL) then
            reorderWithTransverse=[4,3,1,2,5]
        endif
        this=TensorTranspose(TensorReshape( tensor ,newDims) ,reorderWithTransverse)
      else
         !Find in which order the tensor has to enter so that the spin is in the last dimension
         reorderedDims=[1,2,3]
         reorderedDims(whichDimIsSpin)=3
         reorderedDims(3)=whichDimIsSpin
         if (DirectionOfBonds.eq.HORIZONTAL) then
            newDims(1)=oldDims(reorderedDims(1))
            newDims(2)=oldDims(reorderedDims(2))
            newDims(3)=integerONE
            newDims(4)=integerONE
            newDims(5)=oldDims(reorderedDims(3))
         else if (DirectionOfBonds.eq.VERTICAL) then
            newDims(1)=integerONE
            newDims(2)=integerONE
            newDims(3)=oldDims(reorderedDims(1))
            newDims(4)=oldDims(reorderedDims(2))
         endif
         this=TensorReshape( TensorTranspose(tensor,reorderedDims) , newDims)
      endif

   end function new_PEPSTensor_fromTensor3
!##################################################################
!###########       Accessor methods
!##################################################################
   integer function spin_PEPSTensor(this) result(s)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

    if(.not.(this%IsInitialized())) then
        call ThrowException('spin_PEPSTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        s=dims(5)
     endif

   end function spin_PEPSTensor
!##################################################################
   function Get_PEPSBondList(this) result(bonds)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)
     integer :: bonds(4)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DLeft_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        bonds=dims(1:4)
     endif

   end function Get_PEPSBondList
!##################################################################
   integer function DLeft_PEPSTensor(this) result(DL)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DLeft_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DL=dims(1)
     endif

   end function DLeft_PEPSTensor
!##################################################################
   integer function DRight_PEPSTensor(this) result(DR)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DR=dims(2)
     endif

   end function DRight_PEPSTensor

!##################################################################

   integer function DUp_PEPSTensor(this) result(DU)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DUp_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DU=dims(3)
     endif

   end function DUp_PEPSTensor

!##################################################################
   integer function DDown_PEPSTensor(this) result(DD)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DDown_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DD=dims(4)
     endif

   end function DDown_PEPSTensor
!##################################################################

   integer function Get_MaxBondDimensionPEPSTensor(this) result(maxDimension)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('Get_MaxBondDimensionPEPSTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        maxDimension=maxval(dims)
     endif

   end function Get_MaxBondDimensionPEPSTensor

       subroutine Print_PEPSTensor_Dimensions(this,message)
        class(PEPSTensor) :: this
        character*(*),optional :: message
        integer :: dims(5)

        if((this%IsInitialized())) then
            if (present(message) ) write(*,'(A)') message
            dims=this%GetDimensions()
            write(*,'("Spin: ",I3)') ,dims(5)
            write(*,'("DLeft: ",I3,", DRight: ",I3)') dims(1),dims(2)
            write(*,'("DUp: ",I3,", DDown: ",I3)') dims(3),dims(4)
        else
            call ThrowException('Print_PEPSDims','Tensor not initialized',NoErrorCode,Warning)
            return
        endif

   end subroutine Print_PEPSTensor_Dimensions


!##################################################################
!#######################        Products by things
!##################################################################

   function Apply_Operator_To_PEPS_Spin_Dimension(this,anOperator) result(aTensor)
     class(PEPSTensor),intent(IN) :: this
     type(PEPSTensor) :: aTensor
     class(SpinOperator),intent(IN) :: anOperator
     integer :: opDims(2),tensorDims(5)

     opDims=anOperator%GetDimensions()
     if(this%IsInitialized()) then
        if(opDims(1).eq.this%getSpin().and.opDims(2).eq.this%getSpin()) then
           aTensor = new_Tensor(this*anOperator)
        else
           call ThrowException('Apply_Operator_To_Spin_Dimension','Operator is not of the right size',opDims(1)-opDims(2),CriticalError)
        endif
     else
        call ThrowException('Apply_Operator_To_Spin_Dimension','Tensor not initialized',NoErrorCode,CriticalError)
     endif

   end function Apply_Operator_To_PEPS_Spin_Dimension



   function Apply_Matrix_To_PEPS(this,aMatrix,whichDimension) result(aTensor)
     class(PEPSTensor),intent(IN) :: this
     type(PEPSTensor) :: aTensor
     class(Tensor2),intent(IN) :: aMatrix
     integer :: whichDimension

     if(this%IsInitialized()) then
        select case(whichDimension)
            case(LEFT)
                aTensor=nModeProduct(aMatrix,this,FIRST)
            case(RIGHT)
                aTensor=nModeProduct(aMatrix,this,SECOND)
            case(UP)
                aTensor=nModeProduct(aMatrix,this,THIRD)
            case(DOWN)
                aTensor=nModeProduct(aMatrix,this,FOURTH)
            case default
                call ThrowException('Apply_Matrix_To_PEPS','dim must be LEFT/RIGHT/UP/DOWN',whichDimension,CriticalError)
        end select
     else
        call ThrowException('Apply_Matrix_To_PEPS','Tensor not initialized',NoErrorCode,CriticalError)
     endif

   end function Apply_Matrix_To_PEPS


   function Apply_Tensor4_To_PEPS(this,aTensor4,whichBondInPEPS) result(aNewPEPS)
      class(PEPSTensor),intent(IN) :: this
      class(Tensor4),intent(IN) :: aTensor4 !Tensor is assumed in convention L/R/U/D
      integer,intent(IN) :: whichBondInPEPS
      type(PEPSTensor) :: aNewPEPS

      if(this%IsInitialized()) then
        select case(whichBondInPEPS)
            case(LEFT)
               aNewPEPS=new_PEPSTensor( TensorTranspose ( MultAndCollapse ( TensorTranspose(aTensor4,[4,1,2,3]), TensorTranspose(this,[1,4,2,3,5]) ), [3,4,1,2,5]) )
            case(RIGHT)
                aNewPEPS=new_PEPSTensor( TensorTranspose( MultAndCollapse ( TensorTranspose(aTensor4,[1,4,2,3]), TensorTranspose(this,[4,1,2,3,5]) ), [3,4,2,1,5]) )
            case(UP)
                aNewPEPS=new_PEPSTensor(                  MultAndCollapse ( TensorTranspose(aTensor4,[2,3,4,1]), TensorTranspose(this,[2,3,1,4,5]) ) )
            case(DOWN)
                aNewPEPS=new_PEPSTensor( TensorTranspose( MultAndCollapse ( TensorTranspose(aTensor4,[2,3,1,4]), TensorTranspose(this,[2,3,4,1,5]) ), [1,2,4,3,5]) )
            case default
                call ThrowException('Apply_Matrix_To_PEPS','dim must be LEFT/RIGHT/UP/DOWN',whichBondInPEPS,CriticalError)
        end select
     else
        call ThrowException('ApplyMPOToBond','Tensor not initialized',NoErrorCode,CriticalError)
     endif

   end function Apply_Tensor4_To_PEPS


!##################################################################
!##################################################################
!##################################################################

    function Collapse_PEPS_Into_Tensor4(this,anOperator) result (aTensor)
        class(PEPSTensor),intent(IN) :: this
        class(SpinOperator),intent(IN),optional :: anOperator
        type(Tensor4) :: aTensor

        if(this%IsInitialized()) then
            if(present(anOperator)) then
                aTensor= MirrorCompact(this, this.apply.anOperator, FIFTH)
            else
                aTensor= MirrorCompact(this, this, FIFTH)
            endif
        else
            call ThrowException('Collapse_PEPS_Into_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function Collapse_PEPS_Into_Tensor4

!##################################################################

    function Collapse_Two_PEPS_into_Tensor4(upPEPS,downPEPS,anOperator) result (aTensor)
        class(PEPSTensor),intent(IN) :: upPEPS,downPEPS
        class(SpinOperator),intent(IN),optional :: anOperator
        type(Tensor4) :: aTensor

        if(upPEPS%IsInitialized().and.downPEPS%IsInitialized()) then
            if(present(anOperator)) then
                aTensor= MirrorCompact(upPEPS.apply.anOperator, downPEPS, FIFTH)
            else
                aTensor= MirrorCompact(upPEPS, downPEPS, FIFTH)
            endif
        else
            call ThrowException('Collapse_Two_PEPS_into_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function Collapse_Two_PEPS_into_Tensor4

!##################################################################

   function CollapsePEPSandSpinIndicesBut(aPEPS,survivingIndex) result(aMAtrix)
      class(PEPSTensor),intent(IN) :: aPEPS
      integer,intent(IN) :: survivingIndex
      type(Tensor2) :: aMatrix
      type(Tensor2) :: tensorAsMatrix

      if(aPEPS%Isinitialized()) then
         select case (survivingIndex)
            case(LEFT)
               tensorAsMatrix=JoinIndicesOf(TensorTranspose(aPEPS,[5,1,2,3,4]))
            case(RIGHT)
               tensorAsMatrix=JoinIndicesOf(TensorTranspose(aPEPS,[1,5,2,3,4]))
            case(UP)
               tensorAsMatrix=JoinIndicesOf(TensorTranspose(aPEPS,[1,2,5,3,4]))
            case(DOWN)
               tensorAsMatrix=JoinIndicesOf(TensorTranspose(aPEPS,[1,2,3,5,4]))
            case default
               call ThrowException('CollapsePEPSandSpinIndicesBut','Index must be L/R/U/D',NoErrorCode,CriticalError)
         end select
         aMatrix= TensorTranspose(tensorAsMatrix) * Conjugate(tensorAsMatrix) !Result matrix is [ index,conjg(index) ]

      else
         call ThrowException('CollapsePEPSandSpinIndicesBut','Tensor not initialized',NoErrorCode,CriticalError)
      endif

   end function CollapsePEPSandSpinIndicesBut

!##################################################################

    function CompactPEPSBondDimensions(aPEPS) result(aMatrix)
        class(PEPSTensor),intent(IN) :: aPEPS
        type(Tensor2) :: aMatrix

        if(aPEPS%Isinitialized()) then
            aMatrix=aPEPS%JoinIndices()
        else
            call ThrowException('CompactPEPSBondDimensions','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function CompactPEPSBondDimensions

!##################################################################

    function JoinPEPSSpinWithBond(aPEPS,whichDirection) result(aTensor)
        class(PEPSTensor),intent(IN) :: aPEPS
        integer,intent(IN) :: whichDirection
        type(Tensor4) :: aTensor

        if(aPEPS%Isinitialized()) then
          select case (whichDirection)
            case (LEFT)
               aTensor=Flatten(aPEPS,[1,5],[2],[3],[4])
            case (RIGHT)
               aTensor=Flatten(aPEPS,[1],[2,5],[3],[4])
            case (UP)
               aTensor=Flatten(aPEPS,[1],[2],[3,5],[4])
            case (DOWN)
               aTensor=Flatten(aPEPS,[1],[2],[3],[4,5])
          end select
        else
            call ThrowException('CompactPEPSBondDimensions','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function JoinPEPSSpinWithBond

!##################################################################
   function DropTwoBondsAndReturnMPSTensor(aPEPS,whichIsLeft,whichIsRight) result(aMPSTensor)
      class(PEPSTensor),intent(IN) :: aPEPS
      integer,intent(IN) :: whichIsLeft,whichIsRight
      type(MPSTensor) :: aMPSTensor
      integer :: n,joinedDims(3),tempInt

        if(aPEPS%Isinitialized()) then
          tempInt=1
          do n=1,4 !loop to skip over the used dimensions and join all the others
            if (n.ne.whichIsLeft .and. n.ne.whichIsRight ) then !5 is the spin dimension...
               joinedDims(tempInt)=n
               tempInt=tempInt+1
            endif
          enddo
          joinedDims(3)=5
          aMPSTensor=new_MPSTensor(Flatten(aPEPS,[whichIsLeft],[whichIsRight],joinedDims),3,1,2) ! 3,1,2 == Spin, Left, Right
        else
            call ThrowException('DropTwoBondsAndReturnMPSTensor','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function DropTwoBondsAndReturnMPSTensor


!##################################################################

    subroutine HighOrderSVDofPEPS(aPEPS,CoreTensor,UMatrices,maxBondDim,whichDirections,SigmaMatrices)
        class(PEPSTensor),intent(IN) :: aPEPS
        integer,intent(IN),optional :: maxBondDim(:),whichDirections(:)
        type(Tensor2) ,intent(OUT):: UMatrices(:)
        type(Tensor2) ,intent(OUT),optional :: SigmaMatrices(:)
        type(PEPSTensor), intent(OUT) :: CoreTensor
        integer :: oldBondDims(4)
        integer :: numOfDimsToSVD

        if(aPEPS%Isinitialized()) then
            if (present(maxBondDim) .and. present(whichDirections) ) then
               if (present(SigmaMatrices)) then
                  call SingularValueDecomposition(aPEPS,CoreTensor,UMatrices, maxBondDim, whichDirections ,SigmaMatrices)
               else
                  call SingularValueDecomposition(aPEPS,CoreTensor,UMatrices, maxBondDim, whichDirections )
               endif
            else if (present(maxBondDim) .and. .not.present(whichDirections)) then
               call SingularValueDecomposition(aPEPS,CoreTensor,UMatrices, maxBondDim, [LEFT, RIGHT, UP, DOWN])
            else
               !Return order is Left-Right-Up-Down
               oldBondDims=aPEPS%GetBonds()
               call SingularValueDecomposition(aPEPS,CoreTensor,UMatrices, oldBondDims, [LEFT, RIGHT, UP, DOWN])
            endif
        else
            call ThrowException('HighOrderSVDofPEPS','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end subroutine HighOrderSVDofPEPS


end module PEPSTensor_Class
