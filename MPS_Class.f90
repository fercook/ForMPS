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

Module MPS_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use MPSTensor_Class

  private

  public :: new_MPS,delete_MPS,LeftCanonizeAtSite,RightCanonizeAtSite,Normalize

!###############################
!#####  The class main object
!###############################
  type,public :: MPS
     private
     integer :: length
     type(MPSTensor), allocatable :: TensorCollection(:)
     logical :: IsInitialized=.false.
     integer :: CanonizedAt=0
   contains
     procedure,public :: LeftCanonizeAtSite => LeftCanonizeMPS_AtSite
     procedure,public :: RightCanonizeAtSite => RightCanonizeMPS_AtSite
     procedure,public :: Canonize => CanonizeMPS
     procedure,public :: Normalize => NormalizeMPS
     procedure,public :: TensorAt => GetMPSTensorAtSite
     procedure,public :: delete => delete_MPS
  end type MPS

  interface new_MPS
    module procedure new_MPS_Random,new_MPS_fromMPS
  end interface

  interface assignment (=)
     module procedure new_MPS_fromAssignment
  end interface

    contains
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function new_MPS_Random(length,spin,bond) result (this)
    integer,intent(IN) :: length,bond,spin
    type(MPS) :: this
    integer :: n

    allocate(this%TensorCollection(0:length+1))
    this%TensorCollection(0)=new_MPSTensor(spin,integerONE,integerONE,ONE)
    this%TensorCollection(length+1)=new_MPSTensor(spin,integerONE,integerONE,ONE)
    this%TensorCollection(1)=new_MPSTensor(spin,integerONE,bond)
    this%TensorCollection(length)=new_MPSTensor(spin,bond,integerONE)
    do n=2,length-1
        this%TensorCollection(n)=new_MPSTensor(spin,bond,bond)
    enddo
    this%Length=length
    this%IsInitialized=.true.

  end function new_MPS_Random

!##################################################################
   function new_MPS_fromMPS (aMPS) result (this)
     class(MPS),intent(in) :: aMPS
     type(MPS) :: this
     integer  :: n,error

     if(.not.aMPS%IsInitialized) then
         call ThrowException('new_MPS_fromMPS','MPS not initialized',NoErrorCode,CriticalError)
     endif
     length=aMPS%length
     if (this%IsInitialized) error=this%delete()
     allocate(this%TensorCollection(0:length+1))
     do n=0,length+1
         this%TensorCollection(n)=new_MPSTensor(aMPS%TensorCollection(n))
     enddo
     this%Length=length
     this%IsInitialized=.true.

   end function new_MPS_fromMPS



   subroutine new_MPS_fromAssignment(lhs,rhs)
     class(MPS),intent(out) :: lhs
     type(MPS),intent(in) :: rhs
     integer  :: n,error

     if(.not.rhs%IsInitialized) then
         call ThrowException('new_MPS_fromAssignment','MPS not initialized',NoErrorCode,CriticalError)
     endif
     length=rhs%length
     if (lhs%IsInitialized) error=lhs%delete()
     allocate(lhs%TensorCollection(0:length+1))
     do n=0,length+1
         lhs%TensorCollection(n)=new_MPSTensor(rhs%TensorCollection(n))
     enddo
     lhs%Length=length
     lhs%IsInitialized=.true.

   end subroutine new_MPS_fromAssignment

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   integer function delete_MPS(this) result(error)
      class(MPS),intent(INOUT) :: this
      integer :: n

     if(.not.this%IsInitialized) then
         call ThrowException('delete_MPS','Tensor is not initalized',error,CriticalError)
     endif
     n=1
     error=Normal
     do while (error.eq.Normal.and.n.le.this%length)
        error=this%TensorCollection(n)%delete()
        n=n+1
     enddo
     if (error.ne.Normal) then
         call ThrowException('delete_MPS','Some error while deleting tensors !',error,CriticalError)
     endif

    end function delete_MPS


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine LeftCanonizeMPS_AtSite(aMPS,site)
    class(MPS),intent(INOUT) :: aMPS
    integer,intent(IN) :: site
    type(Tensor2) :: canonizingMatrix

    if(.not.aMPS%IsInitialized) then
        call ThrowException('LeftCanonizeMPS_AtSite','MPS not initialized',NoErrorCode,CriticalError)
    endif
    if(site.lt.1.or.site.gt.aMPS%length) then
        call ThrowException('LeftCanonizeMPS_AtSite','1<=Site<=Lenght',site,CriticalError)
    endif
    canonizingMatrix=LeftCanonize(aMPS%TensorCollection(site))
    aMPS%TensorCollection(site-1)=MPSTensor_times_matrix(aMPS%TensorCollection(site-1),canonizingMatrix)
  end subroutine LeftCanonizeMPS_AtSite

  subroutine RightCanonizeMPS_AtSite(aMPS,site)
    class(MPS),intent(INOUT) :: aMPS
    integer,intent(IN) :: site
    type(Tensor2) :: canonizingMatrix

    if(.not.aMPS%IsInitialized) then
        call ThrowException('LeftCanonizeMPS_AtSite','MPS not initialized',NoErrorCode,CriticalError)
    endif
    if(site.lt.1.or.site.gt.aMPS%length) then
        call ThrowException('LeftCanonizeMPS_AtSite','1<=Site<=Lenght',site,CriticalError)
    endif
    canonizingMatrix=RightCanonize(aMPS%TensorCollection(site))
    aMPS%TensorCollection(site+1)=matrix_times_MPSTensor(canonizingMatrix,aMPS%TensorCollection(site+1))
  end subroutine RightCanonizeMPS_AtSite
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine CanonizeMPS(aMPS)
    class(MPS),intent(INOUT) :: aMPS
    integer :: n

    if(aMPS%IsInitialized) then
        !First canonize to the right
        do n=1,aMPS%length
            call aMPS%RightCanonizeAtSite(n)
        enddo
        !Now canonize to the left and end at first site
        do n=aMPS%length,1,-1
            call aMPS%LeftCanonizeAtSite(n)
        enddo
    else
        call ThrowException('CanonizeMPS','MPS not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine CanonizeMPS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine NormalizeMPS(aMPS,factor)
    class(MPS),intent(INOUT) :: aMPS
    complex(8), intent(IN) :: factor
    integer :: n

    if(aMPS%IsInitialized) then
        do n=1,aMPS%length
            aMPS%TensorCollection(n)=factor*(aMPS%TensorCollection(n))
        enddo
    else
        call ThrowException('NormalizeMPS','MPS not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine NormalizeMPS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetMPSTensorAtSite(aMPS,site) result(aMPSTensor)
     class(MPS),intent(INOUT) :: aMPS
     integer,intent(IN) :: site
     type(MPSTensor) :: aMPSTensor

     if(aMPS%IsInitialized) then
        if(site.ge.1.or.site.le.length) then
             aMPSTensor=aMPS%TensorCollection(n)
        else
             call ThrowException('GetMPSTensorAtSite','Site is wrong index',site,CriticalError)
        endif
     else
         call ThrowException('GetMPSTensorAtSite','MPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPSTensorAtSite

 end module MPS_Class
