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

module Multiplicator2D_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use MPSTensor_Class
  use MPOTensor_Class
  use MPS_Class
  use MPO_Class
  use Multiplicator_Class
  use MPSAlgorithms_Class
  use PEPSTensor_Class
  use PEPOTensor_Class
  use PEPS_Class
  use PEPO_Class

  implicit none
  !private

  integer,parameter :: DefaultApproximationBond=100

    type, public :: Multiplicator2D
        !private
        integer :: XLength,YLength
        integer :: MaximumApproximationBond = DefaultApproximationBond
        logical :: Initialized=.false.
        complex(8) :: UpperProductOfNorms=ONE,LowerProductOfNorms=ONE
        type(MPS),allocatable :: MPS_Above(:)
        type(MPS),allocatable :: MPS_Below(:)
        type(MPO),allocatable :: RowsAsMPO(:)
        type(Multiplicator_With_MPO) :: RowMultiplicator
        type(PEPS),pointer :: PEPS_Above => null()
        type(PEPS),pointer :: PEPS_Below => null()
        type(PEPO),pointer :: PEPO_Center => null()
        logical :: IsPEPOUsed=.false.
        logical,pointer :: int_HasPEPSChangedAt(:,:) => null()
    contains
        procedure,public :: LeftAt => Multiplicator2D_Left
        procedure,public :: RightAt => Multiplicator2D_Right
        procedure,public :: AboveAt => Multiplicator2D_Above
        procedure,public :: BelowAt => Multiplicator2D_Below
        procedure,public :: FullContraction => Overlap_PEPSAboveBelow
        procedure,public :: SetMaxApproxBond => Set_MaximumApproximationBond
        procedure,public :: GetMaxApproxBond => Get_MaximumApproximationBond
        procedure,private :: LowerMPSAtRow => Multiplicator2D_RowAsLowerMPS
        procedure,private :: UpperMPSAtRow => Multiplicator2D_RowAsUpperMPS
        procedure :: PrepareRowAsMPO => PrepareRowAsMPOLazyEvaluation
        procedure,private :: SplitSpinIntoUpperandLowerBonds => Split_MPSTensor_IntoTensor4
        procedure,private :: HasRowChanged => HasARowOfANY_PEPSChanged
        procedure,private :: HasTensorChanged => HasATensorOfANY_PEPSChanged
        procedure,private :: CheckPointPEPS => CheckpointALLPEPS
        procedure,public :: Reset => Reset_Multiplicator2D
        procedure,public :: Delete => Delete_Multiplicator2D
    end type Multiplicator2D

    interface new_Multiplicator2D
        module procedure new_Multiplicator2D_WithPEPSandPEPO
    end interface

    interface LeftAtSite2D
        module procedure Multiplicator2D_Left
    end interface

    interface RightAtSite2D
        module procedure Multiplicator2D_Right
    end interface

    interface BelowAtSite2D
        module procedure Multiplicator2D_Below
    end interface

    interface AboveAtSite2D
        module procedure Multiplicator2D_Above
    end interface

!##################################################################
!##################################################################

contains

!##################################################################
!##################################################################

  function new_Multiplicator2D_WithPEPSandPEPO(PEPS_A,PEPS_B,PEPO_C,MatrixToTrackChanges) result (this)
    class(PEPS),target,intent(IN) :: PEPS_A
    class(PEPS),target,intent(IN),optional :: PEPS_B
    class(PEPO),target,intent(IN),optional :: PEPO_C
    logical,target,intent(INOUT) :: MatrixToTrackChanges(:,:)
    type(Multiplicator2D) :: this
    integer :: Lengths(2),XLength,YLength


    if (PEPS_A%IsInitialized()) then
        Lengths=PEPS_A%GetSize()
        XLength=Lengths(1)
        YLength=Lengths(2)
        allocate( this%MPS_Above(0:Ylength+1),this%MPS_Below(0:YLength+1),this%RowsAsMPO(0:YLength+1) )
        !Initialize the border mps to template with 1
        this%MPS_Above(0)=new_MPS(XLength)
        this%MPS_Above(Ylength+1)=new_MPS(XLength)
        this%MPS_Below(0)=new_MPS(XLength)
        this%MPS_Below(Ylength+1)=new_MPS(XLength)
        !Set the source PEPSs
        this%PEPS_Above => PEPS_A
        if( present(PEPS_B) ) then
            if (PEPS_B%IsInitialized()) then
                this%PEPS_Below => PEPS_B
            else
                call ThrowException('new_Multiplicator2d','PEPS_B not initialized',NoErrorCode,CriticalError)
                return
            endif
        else
            this%PEPS_Below => PEPS_A
        endif
        if( present(PEPO_C) ) then
            if (PEPO_C%IsInitialized()) then
                this%PEPO_Center => PEPO_C
                this%IsPEPOUsed=.true.
            else
                call ThrowException('new_Multiplicator2d','PEPO_C not initialized',NoErrorCode,CriticalError)
                return
            endif
        else
            this%PEPO_Center => null()
            this%IsPEPOUsed=.false.
        endif
        this%Xlength = Xlength
        this%Ylength = Ylength
        this%int_HasPEPSChangedAt => MatrixToTrackChanges
        this%Initialized = .true.
    else
        call ThrowException('new_Multiplicator2d','PEPS_A not initialized',NoErrorCode,CriticalError)
    endif
  end function new_Multiplicator2D_WithPEPSandPEPO


  subroutine Delete_Multiplicator2D(this)
    class(Multiplicator2D),intent(INOUT) :: this
    integer :: n,error=0

    if(this%Initialized) then
        this%PEPS_Above => null()
        this%PEPS_Below => null()
        this%PEPO_Center => null()
        do n=0,this%Ylength+1
            if(this%MPS_Above(n)%IsInitialized()) error=this%MPS_Above(n)%Delete()
            if(this%MPS_Below(n)%IsInitialized()) error=this%MPS_Below(n)%Delete()
            if(this%RowsAsMPO(n)%IsInitialized()) error=this%RowsAsMPO(n)%Delete()
        enddo
        if(this%RowMultiplicator%IsInitialized()) call this%RowMultiplicator%Delete()
        deallocate(this%MPS_Above,this%MPS_Below, this%RowsAsMPO)
        this%int_HasPEPSChangedAt=> null()

        this%Initialized=.false.
    else
        call ThrowException('Delete_Multiplicator','Multiplicator2D is already deleted',NoErrorCode,Warning)
    endif
  end subroutine Delete_Multiplicator2D

!##################################################################

  subroutine Reset_Multiplicator2D(this,aDirection)
    class(Multiplicator2D),intent(INOUT) :: this
    integer,intent(IN) :: aDirection
    integer :: n,error=0

    if(this%Initialized) then
        select case (aDirection)
            case(UP)
                do n=1,this%Ylength
                    if(this%MPS_Above(n)%IsInitialized()) error=this%MPS_Above(n)%Delete()
                enddo
                this%UpperProductOfNorms=ONE
            case(DOWN)
                do n=1,this%Ylength
                    if(this%MPS_Below(n)%IsInitialized()) error=this%MPS_Below(n)%Delete()
                enddo
                this%LowerProductOfNorms=ONE
            case(LEFT)
                call this%RowMultiplicator%Reset(LEFT)
            case(RIGHT)
                call this%RowMultiplicator%Reset(RIGHT)
            case default
                call ThrowException('Reset_Multiplicator','Direction must be LEFT, RIGHT, UP or DOWN',aDirection,CriticalError)
        end select
    else
        call ThrowException('Reset_Multiplicator','Multiplicator is not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine Reset_Multiplicator2D

!##################################################################

    function Multiplicator2D_Left(this,siteX,siteY) result(Mult_LeftAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_LeftAtSite

        if(this%Initialized) then
            call this%PrepareRowAsMPO(siteY)
            this%RowMultiplicator = new_Multiplicator(this%UpperMPSAtRow(siteY+1),this%LowerMPSAtRow(siteY-1),this%RowsAsMPO(siteY),DONOTCONJUGATE)
            Mult_LeftAtSite=this%SplitSpinIntoUpperandLowerBonds(this%RowMultiplicator%MPSLeftAt(siteX),siteX,siteY, LEFT)
        else
            call ThrowException('Multiplicator_Left','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator2D_Left

!##################################################################
!##################################################################

    function Multiplicator2D_Right(this,siteX,siteY) result(Mult_RightAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_RightAtSite

        if(this%Initialized) then
            call this%PrepareRowAsMPO(siteY)
            this%RowMultiplicator = new_Multiplicator(this%UpperMPSAtRow(siteY+1),this%LowerMPSAtRow(siteY-1),this%RowsAsMPO(siteY),DONOTCONJUGATE)
            Mult_RightAtSite=this%SplitSpinIntoUpperandLowerBonds(this%RowMultiplicator%MPSRightAt(siteX), siteX,siteY,RIGHT)
        else
            call ThrowException('Multiplicator_Right','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator2D_Right

!##################################################################

    function Multiplicator2D_Above(this,siteX,siteY) result(Mult_AboveAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_AboveAtSite
        type(MPS),pointer :: UpperMPS

        if(this%Initialized) then
            UpperMPS => Multiplicator2D_RowAsUpperMPS(this,siteY)
            Mult_AboveAtSite=this%SplitSpinIntoUpperandLowerBonds( UpperMPS%GetTensorAt(siteX),siteX,siteY,UP )
        else
            call ThrowException('Multiplicator_Above','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator2D_Above

!##################################################################

    function Multiplicator2D_Below(this,siteX,siteY) result(Mult_BelowAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_BelowAtSite
        type(MPS),pointer :: LowerMPS

        if(this%Initialized) then
            LowerMPS => Multiplicator2D_RowAsLowerMPS(this,siteY)
            Mult_BelowAtSite=this%SplitSpinIntoUpperandLowerBonds( LowerMPS%GetTensorAt(siteX) ,siteX,siteY,DOWN )
        else
            call ThrowException('Multiplicator_Above','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator2D_Below

!##################################################################

    function Split_MPSTensor_IntoTensor4(this,anMPSTensor,siteX,siteY, aDirection) result(aTensor4)
        class(Multiplicator2D),intent(IN) :: this
        class(MPSTensor),intent(IN) :: anMPSTensor
        integer,intent(IN) :: siteX,siteY
        integer :: aDirection,oppositeDirection
        integer :: UpperBondDimension,LowerBondDimension
        type(Tensor4) :: aTensor4

        if(this%Initialized) then
            select case (aDirection)
                case (LEFT)
                    oppositeDirection=RIGHT
                case (RIGHT)
                    oppositeDirection=LEFT
                case (UP)
                    oppositeDirection=DOWN
                case (DOWN)
                    oppositeDirection=UP
            end select
            LowerBondDimension=this%PEPS_Below%GetBondAt(siteX,siteY,oppositeDirection)
            UpperBondDimension=this%PEPS_Above%GetBondAt(siteX,siteY,oppositeDirection)
            if (this%IsPEPOUsed) then
                UpperBondDimension=UpperBondDimension*this%PEPO_Center%GetBondAt(siteX,siteY,oppositeDirection)
            endif
            aTensor4=anMPSTensor%SplitSpinDimension(UpperBondDimension,LowerBondDimension)
        else
            call ThrowException('Split_MPSTensor_IntoTensor4','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end function Split_MPSTensor_IntoTensor4

!##################################################################

    recursive function Multiplicator2D_RowAsLowerMPS(this,row) result (aRowAsMPSBelow)
        class(Multiplicator2D),intent(INOUT),target  :: this
        integer,intent(IN) :: row
        type(MPS),pointer :: aRowAsMPSBelow
        complex(8) :: normOfnewMPS

        if(this%Initialized) then
            if (.not. this%MPS_Below(row)%IsInitialized() ) then
                call this%PrepareRowAsMPO(row)
                this%MPS_Below(row) = this%LowerMPSAtRow(row-1) .applyMPOTo. this%RowsAsMPO(row)
                call  this%MPS_Below(row)%Canonize(normOfnewMPS)
                this%LowerProductOfNorms=this%LowerProductOfNorms*normOfnewMPS
                this%MPS_Below(row) = Approximate(this%MPS_Below(row), this%MaximumApproximationBond)
            endif
            aRowAsMPSBelow => this%MPS_Below(row)
        else
            call ThrowException('Multiplicator_RowAsLowerMPS','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end function Multiplicator2D_RowAsLowerMPS

!##################################################################

    recursive function Multiplicator2D_RowAsUpperMPS(this,row) result (aRowAsMPSAbove)
        class(Multiplicator2D),intent(INOUT),target  :: this
        integer,intent(IN) :: row
        type(MPS),pointer :: aRowAsMPSAbove
        complex(8) :: normOfnewMPS

        if(this%Initialized) then
            if (.not. this%MPS_Above(row)%IsInitialized() ) then
                call this%PrepareRowAsMPO(row)
                this%MPS_Above(row) = this%RowsAsMPO(row) .applyMPOTo. this%UpperMPSAtRow(row+1)
                call this%MPS_Above(row)%Canonize(normOfnewMPS)
                this%UpperProductOfNorms=this%UpperProductOfNorms*normOfnewMPS
                this%MPS_Above(row) = Approximate( this%MPS_Above(row), this%MaximumApproximationBond)
            endif
            aRowAsMPSAbove => this%MPS_Above(row)
        else
            call ThrowException('Multiplicator_RowAsUpperMPS','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end function Multiplicator2D_RowAsUpperMPS

!##################################################################

    subroutine Set_MaximumApproximationBond(this,newMaxBond)
        class(Multiplicator2D),intent(INOUT) :: this
        integer :: newMaxBond
        if(this%Initialized) then
            this%MaximumApproximationBond=newMaxBond
        else
            call ThrowException('Set_MaximumApproximationBond','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end subroutine Set_MaximumApproximationBond

    integer function Get_MaximumApproximationBond(this)
        class(Multiplicator2D),intent(INOUT) :: this
        if(this%Initialized) then
            Get_MaximumApproximationBond=this%MaximumApproximationBond
        else
            call ThrowException('Get_MaximumApproximationBond','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end function Get_MaximumApproximationBond
!##################################################################

    function Overlap_PEPSAboveBelow(this) result(theOverlap)
        class(Multiplicator2D),intent(INOUT) :: this
        complex(8) :: theOverlap
        type(MPS),pointer :: UpperMPS,LowerMPS

        if(this%Initialized) then
            UpperMPS => Multiplicator2D_RowAsUpperMPS(this,1)
            LowerMPS => Multiplicator2D_RowAsLowerMPS(this,0)
            theOverlap = overlap(UpperMPS,LowerMPS) * this%UpperProductOfNorms * this%LowerProductOfNorms
        else
            call ThrowException('Overlap_PEPSAboveBelow','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end function Overlap_PEPSAboveBelow

!##################################################################    `

!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################

!##################################################################

! THIS PART MIGHT BECOME A HELPER CLASS OR SOMETHING LIKE THAT

  subroutine PrepareRowAsMPOLazyEvaluation(this,row)
    class(Multiplicator2D),intent(INOUT) :: this
    integer,intent(IN) :: row
    integer :: siteX
    type(PEPSTensor) :: TempPEPS
    type(MPOTensor) :: tempMPO

        if(this%Initialized) then
            If (.not.this%RowsAsMPO(row)%IsInitialized() .or. row.eq.0 .or. row .eq.this%YLength+1 ) this%RowsAsMPO(row)=new_MPO(this%XLength)
            if ( this%HasRowChanged(row) ) then
                do siteX=1,this%XLength
                  if (this%HasTensorChanged(siteX,row)) then
                    if (this%IsPEPOUsed) then
                        TempPEPS = this%PEPO_Center%GetTensorAt(siteX,row) .applyTo. this%PEPS_Above%GetTensorAt(siteX,row)
                        tempMPO=CollapsePEPS(tempPEPS,this%PEPS_Below%GetTensorAt(siteX,row))
                    else
                        tempMPO=CollapsePEPS(this%PEPS_Above%GetTensorAt(siteX,row),this%PEPS_Below%GetTensorAt(siteX,row))
                    endif
                    call this%RowsAsMPO(row)%SetTensorAt(siteX,tempMPO)
                    call this%CheckPointPEPS(siteX,row)
                  endif
                enddo
            endif
        else
            call ThrowException('PrepareRowAsMPOLazyEvaluation','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

  end subroutine PrepareRowAsMPOLazyEvaluation

!#################################################################

  logical function HasARowOfANY_PEPSChanged(this,row) result(hasThisRowChanged)
     class(Multiplicator2D),intent(IN) :: this
     integer,intent(IN) :: row

     hasThisRowChanged=any(this%int_HasPEPSChangedAt(:,row))

  end function HasARowOfANY_PEPSChanged

!#################################################################

  logical function HasATensorOfANY_PEPSChanged(this,col,row) result(hasThisTensorChanged)
     class(Multiplicator2D),intent(IN) :: this
     integer,intent(IN) :: col,row

     hasThisTensorChanged=this%int_HasPEPSChangedAt(col,row)

  end function HasATensorOfANY_PEPSChanged
!#################################################################

  subroutine CheckpointALLPEPS(this,col,row)
     class(Multiplicator2D),intent(INOUT) :: this
     integer,intent(IN) :: col,row

     this%int_HasPEPSChangedAt(col,row)=.false.

  end subroutine CheckpointALLPEPS

!#################################################################

end module Multiplicator2D_Class
