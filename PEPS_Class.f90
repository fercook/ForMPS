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

Module PEPS_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use PEPSTensor_Class

  implicit none

!###############################
!#####  The class main object
!###############################
  type,public :: PEPS
     private
     integer :: XLength,YLength,spin,bond
     integer, allocatable :: bondList(:,:,:)
     type(PEPSTensor), allocatable :: TensorCollection(:,:)
     logical :: Initialized=.false.
     logical, public, allocatable :: HasTensorChangedAt(:,:)
   contains
     procedure,public :: ScaleBy => ScalePEPSByFactor
     procedure,public :: GetTensorAt => GetPEPSTensorAtSite
     procedure,public :: SetTensorAt => SetPEPSTensorAtSite
     procedure,public :: delete => delete_PEPS
     procedure,public :: GetSize => GetPEPSSize
     procedure,public :: GetSpin => GetPEPSSpin
     procedure,public :: GetMaxBond => GetMaxPEPSBond
     procedure,public :: GetBondAt => GetPEPSBond
     procedure,public :: SetMaxBond => SetMaxPEPSBond
     procedure,public :: IsInitialized => Is_PEPS_Initialized
     procedure,public :: CheckPointState => SetTensorsToUnchanged
  end type PEPS

  interface new_PEPS
    module procedure new_PEPS_Random,new_PEPS_fromPEPS
  end interface

  interface assignment (=)
     module procedure new_PEPS_fromAssignment
  end interface

    contains

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!
!!    Diagram of PEPS positions:  XLength=5  YLength=4
!!
!!
!!     /\  Y
!!     |
!!     |
!!     |        X
!!     ---------->
!!
!!
!!(1,YLen)    (XLen,YLen)
!!    o---o---o---o
!!    |   |   |   |
!!    o---o---o---o
!!    |   |   |   |
!!    o---o---o---o
!!    |   |   |   |
!!    o---o---o---o
!!    |   |   |   |
!!    o---o---o---o
!!  (1,1)       (XLen,1)
!!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function new_PEPS_Random(XLength,YLength,spin,bond) result (this)
    integer,intent(IN) :: XLength,YLength,bond,spin
    type(PEPS) :: this
    integer :: n,m

    allocate(this%TensorCollection(0:XLength+1,0:YLength+1))
    Allocate(this%BondList(0:XLength+1,0:YLength+1,LEFT:DOWN))
    allocate(this%HasTensorChangedAt(0:XLength+1,0:YLength+1))

    !Outside of boundary terms are unit tensors
    this%TensorCollection(0,:)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%bondList(0,:,:)=integerONE
    this%TensorCollection(XLength+1,:)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%bondList(XLength+1,:,:)=integerONE
    this%TensorCollection(:,0)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%bondList(:,0,:)=integerONE
    this%TensorCollection(:,YLength+1)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%bondList(:,YLength+1,:)=integerONE
    !Corner tensors have two bonds equal to one
    this%TensorCollection(1,1)=new_PEPSTensor(spin,integerONE,bond,bond,integerONE)
    this%bondList(1,1,:)= [ integerONE, bond, bond, integerONE ]
    this%TensorCollection(XLength,1)=new_PEPSTensor(spin,bond,integerONE,bond,integerONE)
    this%bondList(XLength,1,:)= [ bond,integerONE,bond,integerONE ]
    this%TensorCollection(1,YLength)=new_PEPSTensor(spin,integerONE,bond,integerONE,bond)
    this%bondList(1,YLength,:)= [ integerONE,bond,integerONE,bond ]
    this%TensorCollection(XLength,YLength)=new_PEPSTensor(spin,bond,integerONE,integerONE,bond)
    this%bondList(XLength,YLength,:)= [ bond,integerONE,integerONE,bond ]
    !Boundary terms have one bond of dimension one
    do n=2,XLength-1
        this%TensorCollection(n,1)=new_PEPSTensor(spin,bond,bond,bond,integerONE)
        this%bondList(n,1,:)= [ bond,bond,bond,integerONE ]
        this%TensorCollection(n,Ylength)=new_PEPSTensor(spin,bond,bond,integerONE,bond)
        this%bondList(n,YLength,:)= [ bond,bond,integerONE,bond ]
    enddo
    do n=2,YLength-1
        this%TensorCollection(1,n)=new_PEPSTensor(spin,integerONE,bond,bond,bond)
        this%bondList(1,n,:)= [ integerONE,bond,bond,bond ]
        this%TensorCollection(Xlength,n)=new_PEPSTensor(spin,bond,integerONE,bond,bond)
        this%bondList(Xlength,n,:)= [ bond,integerONE,bond,bond ]
    enddo
    !Bulk terms are proper PEPS Tensors
    do m=2,YLength-1
        do n=2,Xlength-1
            this%TensorCollection(n,m)=new_PEPSTensor(spin,bond,bond,bond,bond)
            this%bondList(n,m,:)= [ bond,bond,bond,bond ]
        enddo
    enddo
    this%XLength=Xlength
    this%YLength=Ylength
    this%Spin=spin
    this%bond=bond
    this%HasTensorChangedAt=.true.
    this%Initialized=.true.

  end function new_PEPS_Random

!##################################################################

   function new_PEPS_fromPEPS (aPEPS) result (this)
     class(PEPS),intent(in) :: aPEPS
     type(PEPS) :: this
     integer  :: n,m,error,Xlength,Ylength,spin,bond

     if(.not.aPEPS%Initialized) then
         call ThrowException('new_PEPS_fromPEPS','PEPS not initialized',NoErrorCode,CriticalError)
     endif
     spin=aPEPS%spin
     bond=aPEPS%bond
     Xlength=aPEPS%XLength
     Ylength=aPEPS%YLength
     if (this%Initialized) error=this%delete()
     allocate(this%TensorCollection(0:Xlength+1,0:Ylength+1))
     do m=0,Xlength+1
        do n=0,Xlength+1
            this%TensorCollection(n,m)=new_PEPSTensor(aPEPS%TensorCollection(n,m))
        enddo
     enddo
     this%XLength=Xlength
     this%YLength=Ylength
     this%Spin=spin
     this%Bond=bond
     allocate(this%BondList(0:XLength+1,0:YLength+1,LEFT:DOWN))
     this%BondList=aPEPS%BondList
     allocate(this%HasTensorChangedAt(0:XLength+1,0:YLength+1))
     this%HasTensorChangedAt=.true.
     this%initialized=.true.

   end function new_PEPS_fromPEPS

   subroutine new_PEPS_fromAssignment(lhs,rhs)
     class(PEPS),intent(out) :: lhs
     type(PEPS),intent(in) :: rhs
     integer  :: n,m,error,length,spin,bond,XLength,YLength

     if(.not.rhs%initialized) then
         call ThrowException('new_PEPS_fromAssignment','PEPS not initialized',NoErrorCode,CriticalError)
     endif
     Xlength=rhs%Xlength
     Ylength=rhs%Ylength
     spin=rhs%spin
     bond=rhs%bond
     if (lhs%initialized) error=lhs%delete()
     allocate(lhs%TensorCollection(0:Xlength+1,0:Ylength+1))
     do m=0,Ylength+1
        do n=0,Xlength+1
            lhs%TensorCollection(n,m)=new_PEPSTensor(rhs%TensorCollection(n,m))
        enddo
     enddo
     lhs%XLength=Xlength
     lhs%YLength=Ylength
     lhs%Spin=spin
     lhs%Bond=bond
     allocate(lhs%BondList(0:XLength+1,0:YLength+1,LEFT:DOWN))
     lhs%BondList=rhs%BondList
     allocate(lhs%HasTensorChangedAt(0:XLength+1,0:YLength+1))
     lhs%HasTensorChangedAt=.true.
     lhs%initialized=.true.

   end subroutine new_PEPS_fromAssignment

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   integer function delete_PEPS(this) result(error)
      class(PEPS),intent(INOUT) :: this
      integer :: n,m

     if(.not.this%initialized) then
         call ThrowException('delete_PEPS','Tensor is not initalized',error,CriticalError)
     endif

     do m=0,this%YLength+1
        do n=0,this%XLength+1
            error=this%TensorCollection(n,m)%delete()
            if (error.ne.Normal) then
                call ThrowException('delete_PEPS','Some error while deleting tensors !',error,CriticalError)
            endif
        enddo
     enddo
     this%Xlength=0
     this%Ylength=0
     deallocate(this%BondList)
     deallocate(this%HasTensorChangedAt)

     this%Initialized=.false.

    end function delete_PEPS


   logical function Is_PEPS_initialized(this) result(AmIInitialized)
    class(PEPS) :: this

    AmIInitialized=this%Initialized

   end function Is_PEPS_initialized

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine ScalePEPSByFactor(aPEPS,factor)
    class(PEPS),intent(INOUT) :: aPEPS
    complex(8), intent(IN) :: factor
    integer :: n,m

    if(aPEPS%Initialized) then
        do m=1,aPEPS%Ylength
        do n=1,aPEPS%Xlength
            aPEPS%TensorCollection(n,m) = factor* (aPEPS%TensorCollection(n,m))
        enddo
        enddo
    else
        call ThrowException('ScalePEPSByFactor','PEPS not initialized',NoErrorCode,CriticalError)
    endif

  end subroutine ScalePEPSByFactor

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetPEPSTensorAtSite(aPEPS,Xposition,Yposition) result(aPEPSTensor)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: Xposition,Yposition
     type(PEPSTensor) :: aPEPSTensor

     if(aPEPS%Initialized) then
        if(Xposition.ge.1 .and. Xposition.le.aPEPS%Xlength .and. Yposition.ge.1 .and. Yposition.le.aPEPS%Ylength ) then
             aPEPSTensor=aPEPS%TensorCollection(Xposition,Yposition)
        else
             call ThrowException('GetPEPSTensorAtSite','Site is wrong index',XPosition,CriticalError)
        endif
     else
         call ThrowException('GetPEPSTensorAtSite','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSTensorAtSite

   subroutine SetPEPSTensorAtSite(thisPEPS,Xposition,Yposition,aPEPSTensor)
     class(PEPS),intent(INOUT) :: thisPEPS
     integer,intent(IN) :: Xposition,Yposition
     class(PEPSTensor),intent(IN) :: aPEPSTensor

     if(thisPEPS%Initialized.and.aPEPSTensor%IsInitialized()) then
        if(Xposition.ge.1 .and. Xposition.le.thisPEPS%Xlength .and. Yposition.ge.1 .and. Yposition.le.thisPEPS%Ylength ) then
             thisPEPS%TensorCollection(Xposition,Yposition)=aPEPSTensor
             !Now update Bond list
             thisPEPS%BondList(Xposition,Yposition,:)=[ aPEPSTensor%getDLeft(), aPEPSTensor%getDRight(), aPEPSTensor%getDUp(), aPEPSTensor%getDDown() ]
             thisPEPS%bond=max(thisPEPS%bond,maxval(thisPEPS%BondList(Xposition,Yposition,:)) )
             thisPEPS%HasTensorChangedAt(Xposition,Yposition)=.true.
        else
             call ThrowException('SetPEPSTensorAtSite','Site is wrong index',Xposition,CriticalError)
        endif
     else
         call ThrowException('SetPEPSTensorAtSite','PEPS or Tensor not initialized',NoErrorCode,CriticalError)
     endif
   end Subroutine SetPEPSTensorAtSite

   function GetPEPSSize(aPEPS) result(TheSize)
     class(PEPS),intent(IN) :: aPEPS
     integer :: TheSize(2)

     if(aPEPS%Initialized) then
         TheSize = [ aPEPS%Xlength, aPEPS%YLength ]
     else
         call ThrowException('GetPEPSLength','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSSize


   integer function GetPEPSSpin(aPEPS) result(spin)
     class(PEPS),intent(IN) :: aPEPS
     if(aPEPS%Initialized) then
         spin = aPEPS%spin
     else
         call ThrowException('GetPEPSSpin','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSSpin

   integer function GetMAXPEPSBond(aPEPS) result(bond)
     class(PEPS),intent(IN) :: aPEPS
     if(aPEPS%Initialized) then
         bond = aPEPS%bond
     else
         call ThrowException('GetPEPSBond','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMaxPEPSBond

   integer function GetPEPSBond(aPEPS,siteX,siteY,aDirection) result(bond)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: siteX,siteY,aDirection

     if(aPEPS%Initialized) then
         bond = aPEPS%BondList(siteX,siteY,aDirection)
     else
         call ThrowException('GetPEPSBond','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSBond

   subroutine SetTensorsToUnchanged(this,Xpos,Ypos)
      class(PEPS),intent(INOUT) :: this
      integer,intent(IN) :: Xpos
      integer,intent(IN),optional :: Ypos

     if(this%Initialized) then
        if(Xpos.eq.ALLTENSORS) then
           this%HasTensorChangedAt=.false.
        else if (present(Ypos) ) then
           this%HasTensorChangedAt(Xpos,Ypos)=.false.
        else !Xpos is interpreted as a row...quite lame
           this%HasTensorChangedAt(:,Xpos)=.false.
        endif
     else
         call ThrowException('SetTensorsToUnchanged','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end subroutine SetTensorsToUnchanged

   function SetMAXPEPSBond(thisPEPS,newMaxBond) result(newPEPS)
     class(PEPS),intent(IN) :: thisPEPS
     integer,intent(IN) :: newMaxBond
     type(PEPS) :: newPEPS
     integer :: X,Y
     type(PEPSTensor) :: coreTensor,tempTensor
     type(Tensor2) :: UMatrices(4)

     if(aPEPS%Initialized) then
         if (newMaxBond.ge.thisPEPS%GetMaxBond()) then
            newPEPS=thisPEPS
         else
            newPEPS=new_PEPS(thisPEPS%XLength,thisPEPS%YLength,thisPEPS%spin,newMaxBond)
            do Y=1,thisPEPS%YLength
                do X=1,thisPEPS%XLength
                    call thisPEPS%TensorCollection(X,Y)%HOSVD(coreTensor,Umatrices,newMaxBond)
                    tempTensor=
                    call newPEPS%SetTensorAt(X,Y,
         endif
     else
         call ThrowException('SetPEPSBond','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function SetMaxPEPSBond

 end module PEPS_Class


