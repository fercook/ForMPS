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

Module MPO_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use MPSTensor_Class
  use MPOTensor_Class
  use MPS_Class

!  private
  implicit none

!###############################
!#####  The class main object
!###############################
  type,public :: MPO
     private
     integer :: length,spin
     integer, allocatable :: BondList(:,:)
     type(MPOTensor), allocatable :: TensorCollection(:)
     logical :: Initialized=.false.
   contains
     procedure,public :: GetTensorAt => GetMPOTensorAtSite
     procedure,public :: SetTensorAt => SetMPOTensorAtSite
     procedure,public :: delete => delete_MPO
     procedure,public :: GetSize => GetMPOLength
     procedure,public :: GetSpin => GetMPOSpin
     procedure,public :: GetBond => GetMPOBond
     procedure,public :: IsInitialized => Is_MPO_Initialized
  end type MPO

  interface new_MPO
    module procedure new_MPO_Random,new_MPO_fromMPO,new_MPO_Template
  end interface

  interface assignment (=)
     module procedure new_MPO_fromAssignment
  end interface

  interface operator (.TensorAt.)
    module procedure GetMPOTensorAtSite
  end interface

  interface operator (.applyMPOTo.)
    module procedure Apply_MPO_To_MPS,Apply_MPS_To_MPO
  end interface

    contains
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function new_MPO_Random(length,spin,bond) result (this)
    integer,intent(IN) :: length,bond,spin
    type(MPO) :: this
    integer :: n

    allocate(this%TensorCollection(0:length+1))
    this%TensorCollection(0)=new_MPOTensor(spin,integerONE,integerONE)
    this%TensorCollection(length+1)=new_MPOTensor(spin,integerONE,integerONE)
    this%TensorCollection(1)=new_MPOTensor(spin,integerONE,bond)
    this%TensorCollection(length)=new_MPOTensor(spin,bond,integerONE)
    do n=2,length-1
        this%TensorCollection(n)=new_MPOTensor(spin,bond,bond)
    enddo
    this%Length=length
    this%Spin=spin
    allocate(this%BondList(0:length+1,LEFT:RIGHT))
    do n=0,length+1
        this%BondList(n,LEFT)=this%TensorCollection(n)%GetDLeft()
        this%BondList(n,RIGHT)=this%TensorCollection(n)%GetDRight()
    enddo
    this%Initialized=.true.

  end function new_MPO_Random

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function new_MPO_Template(length) result (this)
    integer,intent(IN) :: length
    type(MPO) :: this
    integer :: n,bond=integerONE,spin=integerONE

    allocate(this%TensorCollection(0:length+1))
    do n=0,length+1
        this%TensorCollection(n)=new_MPOTensor(spin,bond,bond)
    enddo
    this%Length=length
    this%Spin=spin
    allocate(this%BondList(0:length+1,LEFT:RIGHT))
    do n=0,length+1
        this%BondList(n,LEFT)=this%TensorCollection(n)%GetDLeft()
        this%BondList(n,RIGHT)=this%TensorCollection(n)%GetDRight()
    enddo
    this%Initialized=.true.

  end function new_MPO_Template
!##################################################################
   function new_MPO_fromMPO (aMPO) result (this)
     class(MPO),intent(in) :: aMPO
     type(MPO) :: this
     integer  :: n,error,length,spin,bond

     if(.not.aMPO%Initialized) then
         call ThrowException('new_MPO_fromMPO','MPO not initialized',NoErrorCode,CriticalError)
     endif
     length=aMPO%length
     spin=aMPO%spin
     if (this%Initialized) error=this%delete()
     allocate(this%TensorCollection(0:length+1))
     do n=0,length+1
         this%TensorCollection(n)=new_MPOTensor(aMPO%TensorCollection(n))
     enddo
     this%Length=length
     this%Spin=spin
     allocate(this%BondList(0:length+1,LEFT:RIGHT))
     this%BondList=aMPO%BondList
     this%initialized=.true.

   end function new_MPO_fromMPO

   subroutine new_MPO_fromAssignment(lhs,rhs)
     class(MPO),intent(out) :: lhs
     type(MPO),intent(in) :: rhs
     integer  :: n,error,length,spin,bond

     if(.not.rhs%initialized) then
         call ThrowException('new_MPO_fromAssignment','MPO not initialized',NoErrorCode,CriticalError)
     endif
     length=rhs%length
     spin=rhs%spin
     if (lhs%initialized) error=lhs%delete()
     allocate(lhs%TensorCollection(0:length+1))
     do n=0,length+1
         lhs%TensorCollection(n)=new_MPOTensor(rhs%TensorCollection(n))
     enddo
     lhs%Length=length
     lhs%Spin=spin
     allocate(lhs%BondList(0:length+1,LEFT:RIGHT))
     lhs%BondList=rhs%BondList
     lhs%initialized=.true.

   end subroutine new_MPO_fromAssignment

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   integer function delete_MPO(this) result(error)
      class(MPO),intent(INOUT) :: this
      integer :: n

     if(.not.this%initialized) then
         call ThrowException('delete_MPO','Tensor is not initalized',error,CriticalError)
     endif
     n=0
     error=Normal
     do while (error.eq.Normal.and.n.le.this%length+1)
        error=this%TensorCollection(n)%delete()
        n=n+1
     enddo
     if (error.ne.Normal) then
         call ThrowException('delete_MPO','Some error while deleting tensors !',error,CriticalError)
     endif
     deallocate(this%BondList)
     this%length=0
     this%Initialized=.false.

    end function delete_MPO


   logical function Is_MPO_initialized(this) result(AmIInitialized)
    class(MPO) :: this

    AmIInitialized=this%Initialized

   end function Is_MPO_initialized

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetMPOTensorAtSite(aMPO,site) result(aMPOTensor)
     class(MPO),intent(IN) :: aMPO
     integer,intent(IN) :: site
     type(MPOTensor) :: aMPOTensor

     if(aMPO%Initialized) then
        if(site.ge.1.or.site.le.aMPO%length) then
             aMPOTensor=aMPO%TensorCollection(site)
        else
             call ThrowException('GetMPOTensorAtSite','Site is wrong index',site,CriticalError)
        endif
     else
         call ThrowException('GetMPOTensorAtSite','MPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPOTensorAtSite
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine SetMPOTensorAtSite(thisMPO,site,aMPOTensor)
     class(MPO),intent(INOUT) :: thisMPO
     integer,intent(IN) :: site
     class(MPOTensor),intent(IN) :: aMPOTensor

     if(thisMPO%Initialized.and.aMPOTensor%IsInitialized()) then
        if(site.ge.1.or.site.le.thisMPO%length) then
             thisMPO%TensorCollection(site)=aMPOTensor
        else
             call ThrowException('SetMPOTensorAtSite','Site is wrong index',site,CriticalError)
        endif
     else
         call ThrowException('SetMPOTensorAtSite','MPO or Tensor not initialized',NoErrorCode,CriticalError)
     endif
   end Subroutine SetMPOTensorAtSite
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   integer function GetMPOLength(aMPO) result(length)
     class(MPO),intent(IN) :: aMPO
     if(aMPO%Initialized) then
         length = aMPO%length
     else
         call ThrowException('GetMPOLength','MPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPOLength


   integer function GetMPOSpin(aMPO) result(spin)
     class(MPO),intent(IN) :: aMPO
     if(aMPO%Initialized) then
         spin = aMPO%spin
     else
         call ThrowException('GetMPOSpin','MPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPOSpin

   integer function GetMPOBond(aMPO,site,aDirection) result(bond)
     class(MPO),intent(IN) :: aMPO
     integer,intent(IN) :: site
     integer :: aDirection
     if(aMPO%Initialized) then
        !next line should be (aMPO%TensorCollection(1))%GetSpin()
         bond = aMPO%BondList(site,aDirection)
     else
         call ThrowException('GetMPOBond','MPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPOBond

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function Apply_MPO_To_MPS(anMPO,anMPS) result(this)
        class(MPO),intent(IN) :: anMPO
        class(MPS),intent(IN) :: anMPS
        type(MPS) :: this
        type(MPSTensor) ::localTensor
        integer :: Length, site

        if(anMPO%IsInitialized().and.anMPS%IsInitialized()) then
            if(anMPO%getSize().eq.anMPS%GetSize()) then
                this=new_MPS(anMPS)
                do site=1,anMPS%getSize()
                    call this%SetTensorAt(site, (anMPO%TensorCollection(site)) .applyTo. (anMPS.TensorAt.site) )
                enddo
            else
                call ThrowException('Apply MPO to MPS','MPO and MPS of different size',NoErrorCode,CriticalError)
            endif
        else
            call ThrowException('Apply MPO to MPS','MPO or MPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Apply_MPO_To_MPS


    function Apply_MPS_To_MPO(anMPS,anMPO) result(this)
        class(MPS),intent(IN) :: anMPS
        class(MPO),intent(IN) :: anMPO
        type(MPS) :: this
        type(MPSTensor) ::localTensor
        integer :: Length, site

        if(anMPO%IsInitialized().and.anMPS%IsInitialized()) then
            if(anMPO%getSize().eq.anMPS%GetSize()) then
                this=new_MPS(anMPS)
                do site=1,anMPS%getSize()
                    call this%SetTensorAt(site, (anMPS.TensorAt.site) .applyTo. (anMPO%TensorCollection(site))  )
                enddo
            else
                call ThrowException('Apply MPS to MPO','MPO and MPS of different size',NoErrorCode,CriticalError)
            endif
        else
            call ThrowException('Apply MPS to MPO','MPO or MPS not initialized',NoErrorCode,CriticalError)
        endif

    end function Apply_MPS_To_MPO



 end module MPO_Class

