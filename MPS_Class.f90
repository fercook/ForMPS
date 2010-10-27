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

  implicit none
!  private

!  public :: new_MPS,delete_MPS,LeftCanonizeAtSite,RightCanonizeAtSite,ScaleBy

!###############################
!#####  The class main object
!###############################
  type,public :: MPS
     private
     integer :: length,spin,bond
     complex(8) :: Norm=ONE
     integer, allocatable :: bondList(:,:)
     type(MPSTensor), allocatable :: TensorCollection(:)
     logical :: Initialized=.false.
   contains
     procedure,public :: LeftCanonizeAtSite => LeftCanonizeMPS_AtSite
     procedure,public :: RightCanonizeAtSite => RightCanonizeMPS_AtSite
     procedure,public :: Canonize => CanonizeMPS
     procedure,public :: SetNorm => setMPSNorm
     procedure,public :: GetNorm => GetMPSNorm
     procedure,public :: GetTensorAt => GetMPSTensorAtSite
     procedure,public :: SetTensorAt => SetMPSTensorAtSite
     procedure,public :: delete => delete_MPS
     procedure,public :: GetSize => GetMPSLength
     procedure,public :: GetSpin => GetMPSSpin
     procedure,public :: GetBondAt => GetMPSBond
     procedure,public :: GetMaxBond => GetMaxMPSBond
     procedure,public :: IsInitialized => Is_MPS_Initialized
  end type MPS

  interface new_MPS
    module procedure new_MPS_Random,new_MPS_fromMPS,new_MPS_Template
  end interface

  interface assignment (=)
     module procedure new_MPS_fromAssignment
  end interface

  interface operator (.TensorAt.)
    module procedure GetMPSTensorAtSite
  end interface

    contains
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function new_MPS_Random(length,spin,bond) result (this)
    integer,intent(IN) :: length,bond,spin
    type(MPS) :: this
    integer :: n

    allocate(this%TensorCollection(0:length+1))
    Allocate(this%BondList(0:Length+1,LEFT:RIGHT))
    this%TensorCollection(0)=new_MPSTensor(spin,integerONE,integerONE,ONE)
    this%bondList(0,:)=integerONE
    this%bondList(Length+1,:)=integerONE
    this%TensorCollection(length+1)=new_MPSTensor(spin,integerONE,integerONE,ONE)
    this%TensorCollection(1)=new_MPSTensor(spin,integerONE,bond)
    this%bondList(1,:)=[ integerONE, bond ]
    this%TensorCollection(length)=new_MPSTensor(spin,bond,integerONE)
    this%bondList(length,:)=[ bond,integerONE ]
    do n=2,length-1
        this%TensorCollection(n)=new_MPSTensor(spin,bond,bond)
    enddo
    this%bondList(2:Length-1,:)=bond
    this%Length=length
    this%Spin=spin
    this%bond=bond
    this%Initialized=.true.

  end function new_MPS_Random

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function new_MPS_Template(length) result (this)
    integer,intent(IN) :: length
    type(MPS) :: this
    integer :: n,spin=1

    allocate(this%TensorCollection(0:length+1))
    do n=0,length+1
        this%TensorCollection(0)=new_MPSTensor(spin,integerONE,integerONE,ONE)
    enddo
    this%Length=length
    this%Spin=spin
    this%bond=integerONE
    Allocate(this%BondList(0:Length+1,LEFT:RIGHT))
    this%bondList=integerONE
    this%Initialized=.true.

  end function new_MPS_Template

!##################################################################
   function new_MPS_fromMPS (aMPS) result (this)
     class(MPS),intent(in) :: aMPS
     type(MPS) :: this
     integer  :: n,error,length,spin,bond

     if(.not.aMPS%Initialized) then
         call ThrowException('new_MPS_fromMPS','MPS not initialized',NoErrorCode,CriticalError)
     endif
     length=aMPS%length
     spin=aMPS%spin
     bond=aMPS%bond
     if (this%Initialized) error=this%delete()
     allocate(this%TensorCollection(0:length+1))
     do n=0,length+1
         this%TensorCollection(n)=new_MPSTensor(aMPS%TensorCollection(n))
     enddo
     this%Length=length
     this%Spin=spin
     this%Bond=bond
     allocate(this%BondList(0:Length+1,LEFT:RIGHT))
     this%BondList=aMPS%BondList
     this%initialized=.true.

   end function new_MPS_fromMPS

   subroutine new_MPS_fromAssignment(lhs,rhs)
     class(MPS),intent(out) :: lhs
     type(MPS),intent(in) :: rhs
     integer  :: n,error,length,spin,bond

     if(.not.rhs%initialized) then
         call ThrowException('new_MPS_fromAssignment','MPS not initialized',NoErrorCode,CriticalError)
     endif
     length=rhs%length
     spin=rhs%spin
     bond=rhs%bond
     if (lhs%initialized) error=lhs%delete()
     allocate(lhs%TensorCollection(0:length+1))
     do n=0,length+1
         lhs%TensorCollection(n)=new_MPSTensor(rhs%TensorCollection(n))
     enddo
     lhs%Length=length
     lhs%Spin=spin
     lhs%Bond=bond
     allocate(lhs%BondList(0:Length+1,LEFT:RIGHT))
     lhs%BondList=rhs%BondList
     lhs%initialized=.true.

   end subroutine new_MPS_fromAssignment

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   integer function delete_MPS(this) result(error)
      class(MPS),intent(INOUT) :: this
      integer :: n

     if(.not.this%initialized) then
         call ThrowException('delete_MPS','Tensor is not initalized',error,CriticalError)
     endif
     n=0
     error=Normal
     do while (error.eq.Normal.and.n.le.this%length+1)
        error=this%TensorCollection(n)%delete()
        n=n+1
     enddo
     if (error.ne.Normal) then
         call ThrowException('delete_MPS','Some error while deleting tensors !',error,CriticalError)
     endif
     this%length=0
     deallocate(this%BondList)
     this%Initialized=.false.

    end function delete_MPS


   logical function Is_MPS_initialized(this) result(AmIInitialized)
    class(MPS) :: this

    AmIInitialized=this%Initialized

   end function Is_MPS_initialized

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine LeftCanonizeMPS_AtSite(aMPS,site)
    class(MPS),intent(INOUT) :: aMPS
    integer,intent(IN) :: site
    type(Tensor2) :: canonizingMatrix

    if(.not.aMPS%Initialized) then
        call ThrowException('LeftCanonizeMPS_AtSite','MPS not initialized',NoErrorCode,CriticalError)
    endif
    if(site.lt.1.or.site.gt.aMPS%length) then
        call ThrowException('LeftCanonizeMPS_AtSite','1<=Site<=Lenght',site,CriticalError)
    endif
    canonizingMatrix=LeftCanonize(aMPS%TensorCollection(site))
    if (1.eq.site) then
        call aMPS%SetNorm(aMPS%Norm*(canonizingMatrix%Norm()))
!        aMPS%TensorCollection(0) = new_MPSTensor(aMPS%spin,integerONE,integerONE,ONE)
    else
        aMPS%TensorCollection(site-1)=MPSTensor_times_matrix(aMPS%TensorCollection(site-1),canonizingMatrix)
    endif
    !Now fix the bond list
    aMPS%BondList(site,:)=[ aMPS%TensorCollection(site)%getDLeft(), aMPS%TensorCollection(site)%getDRight() ]
    aMPS%BondList(site-1,:)=[ aMPS%TensorCollection(site-1)%getDLeft(), aMPS%TensorCollection(site-1)%getDRight() ]
    aMPS%bond=max(aMPS%bond,maxval(aMPS%BondList(site,:)),maxval(aMPS%BondList(site-1,:)))
  end subroutine LeftCanonizeMPS_AtSite

  subroutine RightCanonizeMPS_AtSite(aMPS,site)
    class(MPS),intent(INOUT) :: aMPS
    integer,intent(IN) :: site
    type(Tensor2) :: canonizingMatrix

    if(.not.aMPS%Initialized) then
        call ThrowException('RightCanonizeMPS_AtSite','MPS not initialized',NoErrorCode,CriticalError)
    endif
    if(site.lt.1.or.site.gt.aMPS%length) then
        call ThrowException('RightCanonizeMPS_AtSite','1<=Site<=Lenght',site,CriticalError)
    endif
    canonizingMatrix=RightCanonize(aMPS%TensorCollection(site))
    if (aMPS%length.eq.site) then
        call aMPS%SetNorm(aMPS%Norm*(canonizingMatrix%Norm()))
!        aMPS%TensorCollection(aMPS%length+1) = new_MPSTensor(aMPS%spin,integerONE,integerONE,ONE)
    else
        aMPS%TensorCollection(site+1)=matrix_times_MPSTensor(canonizingMatrix,aMPS%TensorCollection(site+1))
    endif
    !now fix bond list
    aMPS%BondList(site,:)=[ aMPS%TensorCollection(site)%getDLeft(), aMPS%TensorCollection(site)%getDRight() ]
    aMPS%BondList(site+1,:)=[ aMPS%TensorCollection(site+1)%getDLeft(), aMPS%TensorCollection(site+1)%getDRight() ]
    aMPS%bond=max(aMPS%bond,maxval(aMPS%BondList(site,:)),maxval(aMPS%BondList(site+1,:)))
  end subroutine RightCanonizeMPS_AtSite
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine CanonizeMPS(aMPS)
    class(MPS),intent(INOUT) :: aMPS
    integer :: n

    if(aMPS%Initialized) then
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

  subroutine SetMPSNorm(aMPS,newNorm)
    class(MPS),intent(INOUT) :: aMPS
    complex(8), intent(IN) :: newNorm

    if(aMPS%Initialized) then
        aMPS%Norm=newNorm
    else
        call ThrowException('SetMPSNorm','MPS not initialized',NoErrorCode,CriticalError)
    endif
  end subroutine SetMPSNorm

  function GetMPSNorm(aMPS)
    class(MPS),intent(INOUT) :: aMPS
    complex(8) :: GetMPSNorm
    if(aMPS%Initialized) then
        GetMPSNorm=aMPS%Norm
    else
        call ThrowException('GetMPSNorm','MPS not initialized',NoErrorCode,CriticalError)
    endif
  end function GetMPSNorm
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetMPSTensorAtSite(aMPS,site) result(aMPSTensor)
     class(MPS),intent(IN) :: aMPS
     integer,intent(IN) :: site
     type(MPSTensor) :: aMPSTensor

     if(aMPS%Initialized) then
        if(site.ge.1.and.site.le.aMPS%length) then
             aMPSTensor=aMPS%TensorCollection(site)
        else
             call ThrowException('GetMPSTensorAtSite','Site is wrong index',site,CriticalError)
        endif
     else
         call ThrowException('GetMPSTensorAtSite','MPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPSTensorAtSite

   subroutine SetMPSTensorAtSite(thisMPS,site,aMPSTensor)
     class(MPS),intent(INOUT) :: thisMPS
     integer,intent(IN) :: site
     class(MPSTensor),intent(IN) :: aMPSTensor

     if(thisMPS%Initialized.and.aMPSTensor%IsInitialized()) then
        if(site.ge.1 .and. site.le.thisMPS%length) then
             thisMPS%TensorCollection(site)=aMPSTensor
             !Now update Bond list
             thisMPS%BondList(site,:)=[ thisMPS%TensorCollection(site)%getDLeft(), thisMPS%TensorCollection(site)%getDRight() ]
             thisMPS%bond=max(thisMPS%bond,maxval(thisMPS%BondList(site,:)))
        else
             call ThrowException('SetMPSTensorAtSite','Site is wrong index',site,CriticalError)
        endif
     else
         call ThrowException('SetMPSTensorAtSite','MPS or Tensor not initialized',NoErrorCode,CriticalError)
     endif
   end Subroutine SetMPSTensorAtSite

   integer function GetMPSLength(aMPS) result(length)
     class(MPS),intent(IN) :: aMPS
     if(aMPS%Initialized) then
         length = aMPS%length
     else
         call ThrowException('GetMPSLength','MPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPSLength


   integer function GetMPSSpin(aMPS) result(spin)
     class(MPS),intent(IN) :: aMPS
     if(aMPS%Initialized) then
         spin = aMPS%spin
     else
         call ThrowException('GetMPSSpin','MPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPSSpin

   integer function GetMaxMPSBond(aMPS) result(bond)
     class(MPS),intent(IN) :: aMPS
     if(aMPS%Initialized) then
        !next line should be (aMPS%TensorCollection(1))%GetSpin()
         bond = aMPS%bond
     else
         call ThrowException('GetMaxMPSBond','MPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMaxMPSBond

   integer function GetMPSBond(aMPS,site,aDirection) result(bond)
     class(MPS),intent(IN) :: aMPS
     integer,intent(IN) :: site,aDirection
     if(aMPS%Initialized) then
        !next line should be (aMPS%TensorCollection(1))%GetSpin()
         bond = aMPS%BondList(site,aDirection)
     else
         call ThrowException('GetMPSBond','MPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMPSBond

 end module MPS_Class



! module m_strategy_pattern
!implicit none
!
!abstract interface
!    !! A generic interface to a subroutine accepting array of integers
!    subroutine generic_function(numbers)
!        integer, dimension(:), intent(in) :: numbers
!    end subroutine
!end interface
!
!type :: Button
!    character(len=20) :: label
!    procedure(generic_function), pointer, nopass :: on_submit
!contains
!    procedure :: init
!end type Button
!
!contains
!
!    subroutine init(self, func, label)
!        class(Button), intent(inout) :: self
!        procedure(generic_function) :: func
!        character(len=*) :: label
!        self%on_submit => func      !! Procedure pointer
!        self%label = label
!    end subroutine init
!
!    subroutine summation(array)
!        integer, dimension(:), intent(in) :: array
!        integer :: total
!        total = sum(array)
!        write(*,*) total
!    end subroutine summation
!
!    subroutine join(array)
!        integer, dimension(:), intent(in) :: array
!        write(*,*) array        !! Just write out the whole array
!    end subroutine join
!
!end module m_strategy_pattern
!
! The following program demonstrates the usage of the module
!program test_strategy
!use m_strategy_pattern
!implicit none
!
!    type(Button) :: button1, button2
!    integer :: i
!
!    call button1%init(summation, "Add them")
!    call button2%init(join, "Join them")
!
!    call button1%on_submit([(i, i=1,10)])   !! Displays 55
!    call button2%on_submit([(i, i=1,10)])   !! Prints out the array
!
!end program test_strategy
