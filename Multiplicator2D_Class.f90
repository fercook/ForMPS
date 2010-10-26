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

    implicit none
    !private

    type, public :: Multiplicator2D
        private
        integer :: XLength,YLength
        logical :: Initialized=.false.
        type(MPS),allocatable :: MPS_Above(:)
        type(MPS),allocatable :: MPS_Below(:)
        type(Multiplicator_With_MPO) :: RowMultiplicator
        integer :: CurrentRow = UNDEFINED
        type(MPO) :: CurrentRowAsMPO
        type(PEPS),pointer :: PEPS_Above => null()
        type(PEPS),pointer :: PEPS_Below => null()
        type(PEPO),pointer :: PEPO_Center => null()
        logical :: IsPEPOUsed=.false.
    contains
        procedure,public :: LeftAt => Multiplicator_Left
        procedure,public :: RightAt => Multiplicator_Right
        procedure,public :: AboveAt => Multiplicator_Above
        procedure,public :: BelowAt => Multiplicator_Below
        procedure,private :: LowerMPSAtRow => Multiplicator_RowAsLowerMPS
        procedure,private :: UpperMPSAtRow => Multiplicator_RowAsUpperMPS
        procedure,private :: SetCurrentRow => SetCurrentRowasMPO
        procedure,private :: SplitUpperandLowerBonds => Split_MPSTensor_IntoTensor4
        procedure,public :: Reset => Reset_Multiplicator
        procedure,public :: Delete => Delete_Multiplicator
    end type Multiplicator2D

    interface new_Multiplicator2D
        module procedure new_Multiplicator2D_WithPEPSandPEPO
    end interface

    interface LeftAtSite
        module procedure Multiplicator_Left
    end interface

    interface RightAtSite
        module procedure Multiplicator_Right
    end interface

    interface BelowAtSite
        module procedure Multiplicator_Below
    end interface

    interface AboveAtSite
        module procedure Multiplicator_Above
    end interface

!##################################################################
!##################################################################

contains

!##################################################################
!##################################################################

  function new_Multiplicator2D_WithPEPSandPEPO(PEPS_A,PEPS_B,PEPO_C) result (this)
    class(PEPS),target,intent(IN) :: PEPS_A
    class(PEPS),target,intent(IN),optional :: PEPS_B
    class(PEPO),target,intent(IN),optional :: PEPO_C
    type(Multiplicator2D) :: this
    integer :: Lengths(2),XLength,YLength

    if (PEPS_A%IsInitialized()) then
        Lengths=PEPS_A%GetSize()
        XLength=Lengths(1)
        YLength=Lengths(2)
        allocate(this%MPS_Above(0:Ylength+1),this%MPS_Below(0:YLength+1))
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
                IsPEPOUsed=.true.
            else
                call ThrowException('new_Multiplicator2d','PEPS_B not initialized',NoErrorCode,CriticalError)
                return
            endif
        else
            this%PEPO_Center => null()
            IsPEPOUsed=.true.
        endif
        this%Xlength = Xlength
        this%Ylength = Ylength
        this%CurrentRow = UNDEFINED
        this%Initialized = .true.
    else
        call ThrowException('new_Multiplicator2d','PEPS not initialized',NoErrorCode,CriticalError)
    endif
  end function new_Multiplicator2D_WithPEPSandPEPO


  subroutine Delete_Multiplicator(this)
    class(Multiplicator2D),intent(INOUT) :: this
    integer :: n,error=0

    if(this%Initialized) then
        this%PEPS_Above => null()
        this%PEPS_Below => null()
        this%PEPO_Center => null()
        do n=0,this%Ylength+1
            if(this%MPS_Above(n)%IsInitialized()) error=this%MPS_Above(n)%Delete()
            if(this%MPS_Below(n)%IsInitialized()) error=this%MPS_Below(n)%Delete()
        enddo
        error=RowMultiplicator%Delete()
        this%CurrentRow = UNDEFINED
        this%Initialized=.false.
    else
        call ThrowException('Delete_Multiplicator','Multiplicator2D is already deleted',NoErrorCode,Warning)
    endif
  end subroutine Delete_Multiplicator

!##################################################################

  subroutine Reset_Multiplicator(this,aDirection)
    class(Multiplicator),intent(INOUT) :: this
    type(Direction),intent(IN) :: aDirection
    integer :: n,error=0

    if(this%Initialized) then
        select case (aDirection)
            case(UP)
                do n=1,this%Ylength
                    if(this%MPS_Above(n)%IsInitialized()) error=this%MPS_Above(n)%Delete()
                enddo
            case(DOWN)
                do n=1,this%Ylength
                    if(this%MPS_Below(n)%IsInitialized()) error=this%MPS_Below(n)%Delete()
                enddo
            case(LEFT)
                error=this%RowMultiplicator%Reset(LEFT)
            case(RIGHT)
                error=this%RowMultiplicator%Reset(RIGHT)
            case default
                call ThrowException('Reset_Multiplicator','Direction must be LEFT, RIGHT, UP or DOWN',Direction,CriticalError)
        end select
    else
        call ThrowException('Reset_Multiplicator','Multiplicator is not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine Reset_Multiplicator


!##################################################################
  subroutine SetCurrentRowasMPO(this,row)
    class(Multiplicator2D),intent(INOUT) :: this
    integer,intent(IN) :: row
    integer :: siteX
    type(PEPSTensor) :: TempPEPS
    type(MPOTensor) :: TempMPO

        if(this%Initialized) then
            CurrentRowAsMPO=new_MPO(this%XLength)
            do siteX=1,this%XLength
                if (this%IsPEPOUsed) then
                    TempPEPS = this%PEPO_Center%GetTensorAt(siteX,row) .applyTo. this%PEPS_Above%GetTensorAt(siteX,row)
                    tempMPO=CollapsePEPS(tempPEPS,this%PEPS_Below%GetTensorAt(siteX,row))
                else
                    tempMPO=CollapsePEPS(this%PEPS_Above%GetTensorAt(siteX,row),this%PEPS_Below%GetTensorAt(siteX,row))
                endif
                call this%CurrentRowAsMPO%SetTensorAt(siteX,tempMPO)
            enddo
        else
            call ThrowException('Multiplicator_Left','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

  end subroutine SetCurrentRowasMPO

!##################################################################

    function Multiplicator_Left(this,siteX,siteY) result(Mult_LeftAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_LeftAtSite

        if(this%Initialized) then
            if (this%CurrentRow.ne.siteY) then
                call this%SetCurrentRow(siteY)
                this%RowMultiplicator = new_Multiplicator(this%UpperMPSAtRow(siteY+1),this%LoweMPSAtRow(siteY-1),this%CurrentRowAsMPO,DONOTCONJUGATE)
            endif
            Mult_LeftAtSite=this%SplitUpperandLowerDims(this%RowMultiplicator%LeftAt(siteX) ,siteX,siteY, LEFT)
        else
            call ThrowException('Multiplicator_Left','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Left

!##################################################################
!##################################################################

    function Multiplicator_Right(this,siteX,siteY) result(Mult_RightAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_RightAtSite

        if(this%Initialized) then
            if (this%CurrentRow.ne.siteY) then
                call this%SetCurrentRow(siteY)
                this%RowMultiplicator = new_Multiplicator(this%UpperMPSAtRow(siteY+1),this%LowerMPSAtRow(siteY-1),this%CurrentRowAsMPO,DONOTCONJUGATE)
            endif
            Mult_RightAtSite=this%SplitUpperandLowerDims(this%RowMultiplicator%RightAt(siteX), siteX,siteY,RIGHT)
        else
            call ThrowException('Multiplicator_Right','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Right

!##################################################################

    function Multiplicator_Above(this,siteX,siteY) result(Mult_AboveAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_AboveAtSite

        if(this%Initialized) then
                Mult_AboveAtSite=this%SplitUpperandLowerDims( this%UpperMPSAtRow(siteY)%GetTensorAt(siteX),siteX,siteY,UP )
        else
            call ThrowException('Multiplicator_Above','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Above

!##################################################################

    function Multiplicator_Below(this,siteX,siteY) result(Mult_BelowAtSite)
        class(Multiplicator2D),intent(INOUT) :: this
        integer,intent(IN) :: siteX,siteY
        type(Tensor4) :: Mult_BelowAtSite

        if(this%Initialized) then
                Mult_BelowAtSite=this%SplitUpperandLowerDims( this%LowerMPSAtRow(siteY)%GetTensorAt(siteX) ,siteX,siteY,DOWN )
        else
            call ThrowException('Multiplicator_Above','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Below

!##################################################################

    function Split_MPSTensor_IntoTensor4(this,anMPSTensor,siteX,siteY, aDirection) result(aTensor4)
        class(Multiplicator2D),intent(IN) :: this
        class(MPSTensor),intent(IN) :: anMPSTensor
        integer,intent(IN) :: siteX,siteY
        type(Direction) :: aDirection
        integer :: UpperBondDimension,LowerBondDimension
        type(Tensor4) :: aTensor4

        if(this%Initialized) then
            select case (aDirection)
                case (LEFT)
                    LowerBondDimension=this%PEPS_B%GetBondAt(siteX,siteY)
                case (RIGHT)
                case (UP)
                case (DOWN)
            end select
            aTensor4=anMPSTensor%SplitSpin
        else
            call ThrowException('Split_MPSTensor_IntoTensor4','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif

    end function Split_MPSTensor_IntoTensor4

!##################################################################

    recursive function Multiplicator_RowAsLowerMPS(this,row) result (aRowAsMPSBelow)
        class(Multiplicator2D),intent(INOUT)  :: this
        integer,intent(IN) :: row
        type(MPS) :: aRowAsMPSBelow


    end function Multiplicator_RowAsLowerMPS

!##################################################################
!##################################################################gd2nv `q1    `

end module Multiplicator2D_Class
