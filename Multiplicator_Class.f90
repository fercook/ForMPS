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
  use MPO_Class

    implicit none
    !private

    type, public :: Multiplicator
        private
        integer :: Length
        logical :: Initialized=.false.
        integer :: IsMPSBelowConjugated
        type(Tensor2),allocatable :: LeftTensors(:)
        type(Tensor2),allocatable :: RightTensors(:)
        type(MPS),pointer :: MPS_Above => null()
        type(MPS),pointer :: MPS_Below => null()
    contains
        procedure,public :: LeftAt => Multiplicator_Left
        procedure,public :: RightAt => Multiplicator_Right
        procedure,public :: Reset => Reset_Multiplicator
        procedure,public :: Delete => Delete_Multiplicator
    end type Multiplicator

    type, extends(Multiplicator) :: Multiplicator_With_MPO
        private
        type(MPO),pointer :: MPO_Center => null()
    contains
        procedure,public :: MPSLeftAt => MultiplicatorMPO_Left
        procedure,public :: MPSRightAt => MultiplicatorMPO_Right
    end type Multiplicator_With_MPO

    interface new_Multiplicator
        module procedure new_Multiplicator_one_MPS,new_Multiplicator_two_MPS, &
            & new_MultiplicatorMPO_one_MPS,new_MultiplicatorMPO_two_MPS
    end interface

    interface LeftAtSite
        module procedure Multiplicator_Left,MultiplicatorMPO_Left
    end interface

    interface RightAtSite
        module procedure Multiplicator_Right,MultiplicatorMPO_Right
    end interface

!##################################################################
!##################################################################

contains

!##################################################################
!##################################################################

  function new_Multiplicator_one_MPS(MPS_A) result (this)
    class(MPS),target,intent(IN) :: MPS_A
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
        this%MPS_Above => MPS_A
        this%MPS_Below => MPS_A
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
            this%MPS_Above => MPS_A
            this%MPS_Below => MPS_B
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
        this%MPS_Above => null()
        this%MPS_Below => null()
        do n=0,this%length+1
            if(this%LeftTensors(n)%IsInitialized()) error=this%LeftTensors(n)%Delete()
            if(this%RightTensors(n)%IsInitialized()) error=this%RightTensors(n)%Delete()
        enddo
	    select type (Typed_this => this)
	        class is (Multiplicator_With_MPO)
	           typed_this%MPO_Center => null()
        end select
        this%Initialized=.false.
    else
        call ThrowException('Delete_Multiplicator','Multiplicator is already deleted',NoErrorCode,Warning)
    endif
  end subroutine Delete_Multiplicator

!##################################################################

  subroutine Reset_Multiplicator(this,aDirection)
    class(Multiplicator),intent(INOUT) :: this
    integer,intent(IN) :: aDirection
    integer :: n,error=0

    if(this%Initialized) then
        select case (aDirection)
            case(LEFT)
                do n=1,this%length
                    if(this%LeftTensors(n)%IsInitialized()) error=this%LeftTensors(n)%Delete()
                enddo
            case(RIGHT)
                do n=1,this%length
                    if(this%RightTensors(n)%IsInitialized()) error=this%RightTensors(n)%Delete()
                enddo
            case default
                call ThrowException('Reset_Multiplicator','Direction must be LEFT or RIGHT',aDirection,CriticalError)
        end select
    else
        call ThrowException('Reset_Multiplicator','Multiplicator is not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine Reset_Multiplicator


!##################################################################

    recursive function Multiplicator_Left(this,site,anOperator) result(Mult_LeftAtSite)
        class(Multiplicator),intent(INOUT) :: this
        integer,intent(IN) :: site
        class(SpinOperator),intent(IN),optional :: anOperator
        type(Tensor2) :: Mult_LeftAtSite
        type(MPSTensor) :: OperatedTensor

        if(this%Initialized) then
            if (this%LeftTensors(site)%IsInitialized()) then
                Mult_LeftAtSite=this%LeftTensors(site)
            else
                select type (Typed_this => this)
                    class is (Multiplicator_With_MPO)
		                OperatedTensor=Typed_this%MPO_Center%GetTensorAt(site).applyTo.Typed_this%MPS_Above%GetTensorAt(site)
                    class is (Multiplicator)
		                if (present(anOperator)) then
		                    OperatedTensor=Typed_this%MPS_Above%GetTensorAt(site).Apply.anOperator
		                else
		                    OperatedTensor=Typed_this%MPS_Above%GetTensorAt(site)
		                endif
                end select
                Mult_LeftAtSite=MPSLeftProduct(Multiplicator_Left(this,site-1), &
                                & OperatedTensor,this%MPS_Below%GetTensorAt(site))
                this%LeftTensors(site)=Mult_LeftAtSite
            endif
        else
            call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Left

!##################################################################
!##################################################################

    recursive function Multiplicator_Right(this,site,anOperator) result(Mult_RightAtSite)
        class(Multiplicator),intent(INOUT) :: this
        integer,intent(IN) :: site
        class(SpinOperator),intent(IN),optional :: anOperator
        type(Tensor2) :: Mult_RightAtSite
        type(MPSTensor) :: OperatedTensor

        if(this%Initialized) then
            if (this%RightTensors(site)%IsInitialized()) then
                Mult_RightAtSite=this%RightTensors(site)
            else
                select type (Typed_this => this)
                    class is (Multiplicator_With_MPO)
                        OperatedTensor=Typed_this%MPO_Center%GetTensorAt(site).applyTo.Typed_this%MPS_Above%GetTensorAt(site)
                    class is (Multiplicator)
                        if (present(anOperator)) then
                            OperatedTensor=Typed_this%MPS_Above%GetTensorAt(site).Apply.anOperator
                        else
                            OperatedTensor=Typed_this%MPS_Above%GetTensorAt(site)
                        endif
                end select
                Mult_RightAtSite=MPSRightProduct(Multiplicator_Right(this,site+1), &
                        & OperatedTensor,this%MPS_Below%GetTensorAt(site))
                this%RightTensors(site)=Mult_RightAtSite
            endif
        else
            call ThrowException('new_Multiplicator_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Multiplicator_Right

!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################
!##################################################################

  function new_MultiplicatorMPO_one_MPS(MPS_A,MPO_C,ShouldConjugate) result (this)
    class(MPS),target,intent(IN) :: MPS_A
    class(MPO),target,intent(IN) :: MPO_C
    integer,intent(IN),optional :: ShouldConjugate
    type(Multiplicator_With_MPO) :: this
    integer :: length

    if (MPS_A%IsInitialized().and.MPO_C%IsInitialized()) then
        length=MPS_A%GetSize()
        if(MPO_C%GetSize().eq.length) then
	        allocate(this%LeftTensors(0:length+1),this%RightTensors(0:length+1))
	        !Initialize the border matrices to 1
	        this%LeftTensors(0)=new_Tensor(integerONE,integerONE,ONE)
	        this%LeftTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
	        this%RightTensors(0)=new_Tensor(integerONE,integerONE,ONE)
	        this%RightTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
	        !Set the source MPSes
	        this%MPS_Above => MPS_A
	        this%MPS_Below => MPS_A
	        this%MPO_Center => MPO_C
	        this%length = length
	        if (present(ShouldConjugate)) then
                this%IsMPSBelowConjugated=ShouldConjugate
            else
                this%IsMPSBelowConjugated=NO
	        endif
	        this%Initialized = .true.
        else
            call ThrowException('new_MultiplicatorMPO_one_MPS','MPO and MPS are of different length',NoErrorCode,CriticalError)
        endif
    else
        call ThrowException('new_MultiplicatorMPO_one_MPS','MPS not initialized',NoErrorCode,CriticalError)
    endif
  end function new_MultiplicatorMPO_one_MPS


  function new_MultiplicatorMPO_two_MPS(MPS_A,MPS_B,MPO_C,ShouldConjugate) result (this)
    class(MPS),target,intent(IN) :: MPS_A,MPS_B
    class(MPO),target,intent(IN) :: MPO_C
    integer,intent(IN),optional :: ShouldConjugate
    type(Multiplicator_With_MPO) :: this
    integer :: length

    if (MPS_A%IsInitialized().and.MPS_B%IsInitialized()) then
        length=MPS_A%GetSize()
        if(length.eq.MPS_B%GetSize().and.length.eq.MPO_C%GetSize()) then
            allocate(this%LeftTensors(0:length+1),this%RightTensors(0:length+1))
            !Initialize the border matrices to 1
            this%LeftTensors(0)=new_Tensor(integerONE,integerONE,ONE)
            this%LeftTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
            this%RightTensors(0)=new_Tensor(integerONE,integerONE,ONE)
            this%RightTensors(length+1)=new_Tensor(integerONE,integerONE,ONE)
            !Set the source MPSes
            this%MPS_Above => MPS_A
            this%MPS_Below => MPS_B
            this%MPO_Center => MPO_C
            this%length = length
            if (present(ShouldConjugate)) then
                this%IsMPSBelowConjugated=ShouldConjugate
            else
                this%IsMPSBelowConjugated=NO
            endif
            this%Initialized = .true.
        else
            call ThrowException('new_MultiplicatorMPO_two_MPS','MPSs have different length',this%length-MPS_B%GetSize(),CriticalError)
        endif
    else
        call ThrowException('new_MultiplicatorMPO_two_MPS','MPS not initialized',NoErrorCode,CriticalError)
    endif

  end function new_MultiplicatorMPO_two_MPS

  function MultiplicatorMPO_Left(this,site) result(Mult_LeftAtSite)
        class(Multiplicator_With_MPO),intent(INOUT) :: this
        integer,intent(IN) :: site
        type(MPSTensor) :: Mult_LeftAtSite
        type(tensor2) ::aTensor
        type(MPOTensor) :: aMPO

        if(this%Initialized) then
            !print *,'Trying to split :',site,this%MPO_Center%GetBond(site,RIGHT)
            aMPO= this%MPO_Center%GetTensorAt(site)
            call aMPO%PrintDimensions('Inside of MPO tensor at site')
            atensor= Multiplicator_Left(this,site)
            call aTensor%PrintDimensions('   from second dim of  :')
            Mult_LeftAtSite=SplitSpinFromBond(Multiplicator_Left(this,site),SECOND,this%MPO_Center%GetBond(site,RIGHT) )
        else
            call ThrowException('MultiplicatorMPO_Left','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function MultiplicatorMPO_Left

!##################################################################
!##################################################################

  function MultiplicatorMPO_Right(this,site) result(Mult_RightAtSite)
        class(Multiplicator_with_MPO),intent(INOUT) :: this
        integer,intent(IN) :: site
        type(MPSTensor) :: Mult_RightAtSite

        if(this%Initialized) then
            Mult_RightAtSite=TensorTranspose( &
              & SplitSpinFromBond(Multiplicator_Right(this,site),FIRST,this%MPO_Center%GetBond(site,LEFT) ), [3,2,1] )
        else
            call ThrowException('MultiplicatorMPO_Right','Multiplicator not initialized',NoErrorCode,CriticalError)
        endif
    end function MultiplicatorMPO_Right

!##################################################################
!##################################################################


end module Multiplicator_Class
