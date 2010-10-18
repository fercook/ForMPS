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

module Multiplicator_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use MPSTensor_Class
  use MPOTensor_Class
  use MPS_Class

    implicit none
    !private

    type, public :: Multiplicator
        private
        integer :: Length
        logical :: Initialized=.false.
        type(Tensor2),allocatable :: LeftTensors(:)
        type(Tensor2),allocatable :: RightTensors(:)
        type(MPS),pointer :: MPS_Normal => null()
        type(MPS),pointer :: MPS_Conjugated => null()
    contains
        !procedure,public :: LeftAtSite => Multiplicator_Left
        !procedure,public :: RighAtSite => Multiplicator_Right
        procedure,public :: Reset => Reset_Multiplicator
        procedure,public :: Delete => Delete_Multiplicator
    end type Multiplicator

    interface new_Multiplicator
        module procedure new_Multiplicator_one_MPS,new_Multiplicator_two_MPS
    end interface

    interface LeftAtSite
        module procedure Multiplicator_Left_Clean,Multiplicator_Left_with_operator
    end interface

    interface RightAtSite
        module procedure Multiplicator_Right_Clean,Multiplicator_Right_with_operator
    end interface

!##################################################################
!##################################################################

contains

!##################################################################
!##################################################################

  function new_Multiplicator_one_MPS(MPS_A) result (this)
    type(MPS),target,intent(IN) :: MPS_A
    type(Multiplicator) :: this
    integer :: length

    if (MPS_A%IsInitialized()) then
        length=MPS_A%GetSize()
        allocate(this%LeftTensors(0:length+1),this%RightTensors(0:length+1))
        !Initialize the border matrices to 1
        this%LeftTensors(0)=new_Tensor(integerONE,integerONE,ONE)
        this%LeftTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
        this%RightTensors(0)=new_Tensor(integerONE,integerONE,ONE)
        this%RightTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
        !Set the source MPSes
        this%MPS_Normal => MPS_A
        this%MPS_Conjugated => MPS_A
        this%length = length
        this%Initialized = .true.
    else
        call ThrowException('new_Multiplicator_one_MPS','MPS not initialized',NoErrorCode,CriticalError)
    endif
  end function new_Multiplicator_one_MPS


  function new_Multiplicator_two_MPS(MPS_A,MPS_B) result (this)
    type(MPS),target,intent(IN) :: MPS_A,MPS_B
    type(Multiplicator) :: this
    integer :: length

    if (MPS_A%IsInitialized().and.MPS_B%IsInitialized()) then
        length=MPS_A%GetSize()
        if(length.eq.MPS_B%GetSize()) then
            allocate(this%LeftTensors(0:length+1),this%RightTensors(0:length+1))
            !Initialize the border matrices to 1
            this%LeftTensors(0)=new_Tensor(integerONE,integerONE,ONE)
            this%LeftTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
            this%RightTensors(0)=new_Tensor(integerONE,integerONE,ONE)
            this%RightTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
            !Set the source MPSes
            this%MPS_Normal => MPS_A
            this%MPS_Conjugated => MPS_B
            this%length = length
            this%Initialized = .true.
        else
            call ThrowException('new_Multiplicator_two_MPS','MPSs have different length',this%length-MPS_B%GetSize(),CriticalError)
        endif
    else
        call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
    endif

  end function new_Multiplicator_two_MPS


  subroutine Delete_Multiplicator(this)
    class(Multiplicator),intent(INOUT) :: this
    integer :: n,error=0

    if(this%Initialized) then
        this%MPS_Normal => null()
        this%MPS_Conjugated => null()
        do n=0,this%length+1
            if(this%LeftTensors(n)%IsInitialized()) error=this%LeftTensors(n)%Delete()
            if(this%RightTensors(n)%IsInitialized()) error=this%RightTensors(n)%Delete()
        enddo
        this%Initialized=.false.
    else
        call ThrowException('Delete_Multiplicator','Multiplicator is already deleted',NoErrorCode,Warning)
    endif
  end subroutine Delete_Multiplicator

!##################################################################

  subroutine Reset_Multiplicator(this,Direction)
    class(Multiplicator),intent(INOUT) :: this
    integer,intent(IN) :: Direction
    integer :: n,error=0

    if(this%Initialized) then
        select case (Direction)
            case(LEFT)
                do n=1,this%length
                    if(this%LeftTensors(n)%IsInitialized()) error=this%LeftTensors(n)%Delete()
                enddo
            case(RIGHT)
                do n=1,this%length
                    if(this%RightTensors(n)%IsInitialized()) error=this%RightTensors(n)%Delete()
                enddo
            case default
                call ThrowException('Reset_Multiplicator','Direction must be LEFT or RIGHT',Direction,CriticalError)
        end select
    else
        call ThrowException('Reset_Multiplicator','Multiplicator is not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine Reset_Multiplicator


!##################################################################
!##################################################################

    recursive function Multiplicator_Left_Clean(this,site) result(Mult_LeftAtSite)
        class(Multiplicator),intent(INOUT) :: this
        integer,intent(IN) :: site
        type(Tensor2) :: Mult_LeftAtSite

        if(this%Initialized) then
            if (this%LeftTensors(site)%IsInitialized()) then
                Mult_LeftAtSite=this%LeftTensors(site)
            else
                this%LeftTensors(site)=MPSLeftProduct(Multiplicator_Left_Clean(this,site-1), &
                        & this%MPS_Normal%GetTensorAt(site),this%MPS_Conjugated%GetTensorAt(site))
                Mult_LeftAtSite=this%LeftTensors(site)
            endif
        else
            call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Left_Clean


    function Multiplicator_Left_with_operator(this,site,anOperator) result(Mult_LeftAtSite)
        class(Multiplicator),intent(INOUT) :: this
        integer,intent(IN) :: site
        class(SpinOperator),intent(IN) :: anOperator
        type(Tensor2) :: Mult_LeftAtSite
        type(MPSTensor) :: OperatedTensor

        if(this%Initialized) then
            if (this%LeftTensors(site)%IsInitialized()) then
                Mult_LeftAtSite=this%LeftTensors(site)
            else
                OperatedTensor=this%MPS_Normal%GetTensorAt(site)
                this%LeftTensors(site)=MPSLeftProduct(Multiplicator_Left_Clean(this,site-1), &
                        & OperatedTensor.Apply.anOperator,this%MPS_Conjugated%GetTensorAt(site))
                Mult_LeftAtSite=this%LeftTensors(site)
            endif
        else
            call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Left_with_operator

!##################################################################

    recursive function Multiplicator_Right_Clean(this,site) result(Mult_RightAtSite)
        class(Multiplicator),intent(INOUT) :: this
        integer,intent(IN) :: site
        type(Tensor2) :: Mult_RightAtSite

        if(this%Initialized) then
            if (this%RightTensors(site)%IsInitialized()) then
                Mult_RightAtSite=this%RightTensors(site)
            else
                this%RightTensors(site)=MPSRightProduct(Multiplicator_Right_Clean(this,site+1), &
                    & this%MPS_Normal%GetTensorAt(site),this%MPS_Conjugated%GetTensorAt(site))
                Mult_RightAtSite=this%RightTensors(site)
            endif
        else
            call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Right_Clean


    function Multiplicator_Right_with_operator(this,site,anOperator) result(Mult_RightAtSite)
        class(Multiplicator),intent(INOUT) :: this
        integer,intent(IN) :: site
        class(SpinOperator),intent(IN) :: anOperator
        type(Tensor2) :: Mult_RightAtSite
        type(MPSTensor) :: OperatedTensor

        if(this%Initialized) then
            if (this%RightTensors(site)%IsInitialized()) then
                Mult_RightAtSite=this%RightTensors(site)
            else
                OperatedTensor=this%MPS_Normal%GetTensorAt(site)
                this%RightTensors(site)=MPSRightProduct(Multiplicator_Right_Clean(this,site+1), &
                    & OperatedTensor.Apply.anOperator,this%MPS_Conjugated%GetTensorAt(site))
                Mult_RightAtSite=this%RightTensors(site)
            endif
        else
            call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Right_with_operator

!##################################################################
!##################################################################

end module Multiplicator_Class
