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

Module PEPO_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use PEPSTensor_Class
  use PEPOTensor_Class
  use PEPS_Class

!  private
  implicit none

!###############################
!#####  The class main object
!###############################
  type,public :: PEPO
     private
     integer :: XLength,YLength,spin,bond
     integer, allocatable :: bondList(:,:)
      type(PEPOTensor), allocatable :: TensorCollection(:,:)
     logical :: Initialized=.false.
   contains
     procedure,public :: GetTensorAt => GetPEPOTensorAtSite
     procedure,public :: SetTensorAt => SetPEPOTensorAtSite
     procedure,public :: delete => delete_PEPO
     procedure,public :: GetSize => GetPEPOSize
     procedure,public :: GetSpin => GetPEPOSpin
     procedure,public :: GetBond => GetPEPOBond
     procedure,public :: IsInitialized => Is_PEPO_Initialized
  end type PEPO

  interface new_PEPO
    module procedure new_PEPO_Random,new_PEPO_fromPEPO
  end interface

  interface assignment (=)
     module procedure new_PEPO_fromAssignment
  end interface

  interface operator (.applyPEPOTo.)
    module procedure Apply_PEPO_To_PEPS
  end interface

    contains
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  function new_PEPO_Random(Xlength,YLength,spin,bond) result (this)
    integer,intent(IN) :: Xlength,YLength,bond,spin
    type(PEPO) :: this
    integer :: n,m

    allocate(this%TensorCollection(0:Xlength+1,0:Ylength+1))
    !Outside of boundary terms are unit tensors
    this%TensorCollection(0,:)=new_PEPOTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%TensorCollection(XLength+1,:)=new_PEPOTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%TensorCollection(:,0)=new_PEPOTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%TensorCollection(:,YLength+1)=new_PEPOTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    !Corner tensors have two bonds equal to one
     this%TensorCollection(1,1)=new_PEPOTensor(spin,integerONE,bond,bond,integerONE)
     this%TensorCollection(XLength,1)=new_PEPOTensor(spin,bond,integerONE,bond,integerONE)
     this%TensorCollection(1,YLength)=new_PEPOTensor(spin,integerONE,bond,integerONE,bond)
     this%TensorCollection(XLength,YLength)=new_PEPOTensor(spin,bond,integerONE,integerONE,bond)
    !Boundary terms have one bond of dimension one
    do n=2,XLength-1
         this%TensorCollection(n,1)=new_PEPOTensor(spin,bond,bond,bond,integerONE)
         this%TensorCollection(n,Ylength)=new_PEPOTensor(spin,bond,bond,integerONE,bond)
    enddo
    do m=2,YLength-1
         this%TensorCollection(1,m)=new_PEPOTensor(spin,integerONE,bond,bond,bond)
         this%TensorCollection(Xlength,m)=new_PEPOTensor(spin,bond,integerONE,bond,bond)
    enddo
    !Bulk terms are proper PEPO Tensors
    do n=2,Xlength-1
        do m=2,YLength-1
             this%TensorCollection(n,m)=new_PEPOTensor(spin,bond,bond,bond,bond)
        enddo
    enddo
    this%XLength=Xlength
    this%YLength=Ylength
    this%Spin=spin
    this%bond=bond
    Allocate(this%BondList(0:XLength+1,0:YLength+1))
    this%bondList(1:XLength,1:YLength)=bond
    this%bondList(0,:)=1
    this%bondList(XLength+1,:)=1
    this%bondList(:,0)=1
    this%bondList(:,YLength+1)=1
    this%Initialized=.true.

  end function new_PEPO_Random

!##################################################################
   function new_PEPO_fromPEPO (aPEPO) result (this)
     class(PEPO),intent(in) :: aPEPO
     type(PEPO) :: this
     integer  :: n,m,error,Xlength,Ylength,spin,bond

     if(.not.aPEPO%Initialized) then
         call ThrowException('new_PEPO_fromPEPO','PEPO not initialized',NoErrorCode,CriticalError)
     endif
     spin=aPEPO%spin
     bond=aPEPO%bond
     Xlength=aPEPO%XLength
     Ylength=aPEPO%YLength
     if (this%Initialized) error=this%delete()
      allocate(this%TensorCollection(0:Xlength+1,0:Ylength+1))
     do n=0,Xlength+1
        do m=0,Xlength+1
             this%TensorCollection(n,m)=new_PEPOTensor(aPEPO%TensorCollection(n,m))
        enddo
     enddo
     this%XLength=Xlength
     this%YLength=Ylength
     this%Spin=spin
     this%Bond=bond
     allocate(this%BondList(0:XLength+1,0:YLength+1))
     this%BondList=aPEPO%BondList
     this%initialized=.true.

   end function new_PEPO_fromPEPO

   subroutine new_PEPO_fromAssignment(lhs,rhs)
     class(PEPO),intent(out) :: lhs
     type(PEPO),intent(in) :: rhs
     integer  :: n,m,error,spin,bond,XLength,YLength

     if(.not.rhs%Initialized) then
         call ThrowException('new_PEPO_fromAssignment','PEPO not initialized',NoErrorCode,CriticalError)
     endif
     spin=rhs%spin
     bond=rhs%bond
     Xlength=rhs%XLength
     Ylength=rhs%YLength
     if (lhs%Initialized) error=lhs%delete()
     allocate(lhs%TensorCollection(0:Xlength+1,0:Ylength+1))
     do n=0,Xlength+1
        do m=0,Xlength+1
             lhs%TensorCollection(n,m)=new_PEPOTensor(rhs%TensorCollection(n,m))
        enddo
     enddo
     lhs%XLength=Xlength
     lhs%YLength=Ylength
     lhs%Spin=spin
     lhs%Bond=bond
     allocate(lhs%BondList(0:XLength+1,0:YLength+1))
     lhs%BondList=rhs%BondList
     lhs%initialized=.true.

   end subroutine new_PEPO_fromAssignment

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   integer function delete_PEPO(this) result(error)
      class(PEPO),intent(INOUT) :: this
      integer :: n,m

     if(.not.this%initialized) then
         call ThrowException('delete_PEPO','Tensor is not initalized',error,CriticalError)
     endif
     error=Normal
     do n=0,this%XLength+1
        do m=0,this%YLength+1
            error=this%TensorCollection(n,m)%delete()
            if (error.ne.Normal) then
                call ThrowException('delete_PEPO','Some error while deleting tensors !',error,CriticalError)
            endif
        enddo
     enddo
     this%Xlength=0
     this%Ylength=0
     deallocate(this%BondList)
     this%Initialized=.false.

    end function delete_PEPO


   logical function Is_PEPO_initialized(this) result(AmIInitialized)
    class(PEPO) :: this

    AmIInitialized=this%Initialized

   end function Is_PEPO_initialized

!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetPEPOTensorAtSite(aPEPO, Xposition,Yposition) result(aPEPOTensor)
     class(PEPO),intent(IN) :: aPEPO
     integer,intent(IN) :: Xposition,Yposition
     type(PEPOTensor) :: aPEPOTensor
     integer n,m

     if(aPEPO%Initialized) then
        if(Xposition.ge.1 .and. Xposition.le.aPEPO%Xlength .and. Yposition.ge.1 .and. Yposition.le.aPEPO%Ylength ) then
              aPEPOTensor=aPEPO%TensorCollection(Xposition,Yposition)
        else
             call ThrowException('GetPEPOTensorAtSite','Site is wrong index',XPosition,CriticalError)
        endif
     else
         call ThrowException('GetPEPOTensorAtSite','PEPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPOTensorAtSite
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine SetPEPOTensorAtSite(thisPEPO, Xposition,Yposition,aPEPOTensor)
     class(PEPO),intent(INOUT) :: thisPEPO
     integer,intent(IN) :: Xposition,Yposition
     type(PEPOTensor) :: aPEPOTensor

     if(thisPEPO%Initialized) then
        if(Xposition.ge.1 .and. Xposition.le.thisPEPO%Xlength .and. Yposition.ge.1 .and. Yposition.le.thisPEPO%Ylength ) then
              thisPEPO%TensorCollection(Xposition,Yposition)=aPEPOTensor
        else
             call ThrowException('GetPEPOTensorAtSite','Site is wrong index',XPosition,CriticalError)
        endif
     else
         call ThrowException('GetPEPOTensorAtSite','PEPO not initialized',NoErrorCode,CriticalError)
     endif
   end Subroutine SetPEPOTensorAtSite
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetPEPOSize(aPEPO) result(TheSize)
     class(PEPO),intent(IN) :: aPEPO
     integer :: TheSize(2)

     if(aPEPO%Initialized) then
         TheSize = [aPEPO%XLength ,aPEPO%YLength]
     else
         call ThrowException('GetPEPOSize','PEPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPOSize


   integer function GetPEPOSpin(aPEPO) result(spin)
     class(PEPO),intent(IN) :: aPEPO
     if(aPEPO%Initialized) then
         spin = aPEPO%spin
     else
         call ThrowException('GetPEPOSpin','PEPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPOSpin

   integer function GetPEPOBond(aPEPO) result(bond)
     class(PEPO),intent(IN) :: aPEPO
     if(aPEPO%Initialized) then
         bond = aPEPO%bond
     else
         call ThrowException('GetPEPOBond','PEPO not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPOBond
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    function Apply_PEPO_To_PEPS(aPEPO,aPEPS) result(this)
        class(PEPO),intent(IN) :: aPEPO
        class(PEPS),intent(IN) :: aPEPS
        type(PEPS) :: this
        type(PEPSTensor) ::localTensor
        integer :: sizeOfPEPS(2),sizeOfPEPO(2), siteX,siteY

        if(aPEPO%IsInitialized().and.aPEPS%IsInitialized()) then
            sizeOfPEPS=aPEPS%GetSize()
            sizeOfPEPO=aPEPO%GetSize()
            if(sizeOfPEPS .equalvector. sizeOfPEPO) then
                this=new_PEPS(aPEPS)
                do siteX=1,sizeOfPEPS(1)
                  do siteY=1,sizeOfPEPS(2)
                    localTensor = aPEPO%TensorCollection(siteX,siteY) .applyTo. GetPEPSTensorAtSite(aPEPS,siteX,siteY)
                    call SetPEPSTensorAtSite(this, siteX,siteY, localTensor )
                  enddo
                enddo
            else
                call ThrowException('Apply PEPO to PEPS','PEPO and PEPS of different size',NoErrorCode,CriticalError)
            endif
        else
            call ThrowException('Apply PEPO to PEPS','PEPO or PEPS not initialized',NoErrorCode,CriticalError)
        endif
    end function Apply_PEPO_To_PEPS


 end module PEPO_Class

