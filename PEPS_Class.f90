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
  use MPS_Class
  use MPSTensor_Class
  use MPSAlgorithms_Class
!
!	use Multiplicator_Class
!	use MPSAlgorithms_Class
  implicit none


  integer,parameter :: MaxBondDefaultForPEPSCanonization = 10

!###############################
!#####  The class main object
!###############################
  type,public :: PEPS
     private
     integer :: XLength,YLength
     type(PEPSTensor), allocatable :: TensorCollection(:,:)
     logical :: Initialized=.false.
   contains
     procedure,public :: ScaleBy => ScalePEPSByFactor
     procedure,public :: GetTensorAt => GetPEPSTensorAtSite
     procedure,public :: SetTensorAt => SetPEPSTensorAtSite
     procedure,public :: delete => delete_PEPS
     procedure,public :: GetSize => GetPEPSSize
     procedure,public :: GetSpin => GetPEPSSpin
     procedure,public :: GetMaxBond => GetMaxPEPSBond
     procedure,public :: GetBondAt => GetPEPSBond
     procedure,public :: ReduceBond => ReduceMaxPEPSBond
     procedure,public :: IsInitialized => Is_PEPS_Initialized
     procedure,public :: PrintBondDimensions => PrintPEPSBondDimensionsMap
     procedure,public :: IsPEPSWellFormed => Integrity_Check_of_bonds
     procedure,public :: CanonizeAt => CanonizePEPSAtSite
     procedure,public :: GetRowAsMPS => GetPEPSRowAsMPS
     procedure,public :: GetColAsMPS => GetPEPSColAsMPS
  end type PEPS

  interface new_PEPS
    module procedure new_PEPS_Random,new_PEPS_fromPEPS
  end interface

  interface assignment (=)
     module procedure new_PEPS_fromAssignment
  end interface

  interface getRowAsMPS
   module procedure GetPEPSRowAsMPS, GetCoreRowAsMPS
  end interface

  interface getColAsMPS
   module procedure GetPEPSColAsMPS, GetCoreColAsMPS
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
    integer :: n,m,error
    complex(8) :: sqrtOfNorm

    if (this%Initialized) error=this%delete()

    allocate(this%TensorCollection(0:XLength+1,0:YLength+1))

    !Insert some attempt at normalization
    sqrtOfNorm=ONE/(spin*bond**4)   !ONE/sqrt(0.5d0**(1.0d0/( 1.0d0*(XLength-1)*(YLength-1)*spin*bond**4+ ((XLength-1)+(YLength-1))*2.0d0*spin*bond**3+4.0d0*spin*bond**2)) )

    !Outside of boundary terms are unit tensors
    this%TensorCollection(0,:)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%TensorCollection(XLength+1,:)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%TensorCollection(:,0)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    this%TensorCollection(:,YLength+1)=new_PEPSTensor(spin,integerONE,integerONE,integerONE,integerONE,ONE)
    !Corner tensors have two bonds equal to one
    this%TensorCollection(1,1)=sqrtOfNorm*new_PEPSTensor(spin,integerONE,bond,bond,integerONE)
    this%TensorCollection(XLength,1)=sqrtOfNorm*new_PEPSTensor(spin,bond,integerONE,bond,integerONE)
    this%TensorCollection(1,YLength)=sqrtOfNorm*new_PEPSTensor(spin,integerONE,bond,integerONE,bond)
    this%TensorCollection(XLength,YLength)=sqrtOfNorm*new_PEPSTensor(spin,bond,integerONE,integerONE,bond)
    !Boundary terms have one bond of dimension one
    do n=2,XLength-1
        this%TensorCollection(n,1)=sqrtOfNorm*new_PEPSTensor(spin,bond,bond,bond,integerONE)
        this%TensorCollection(n,Ylength)=sqrtOfNorm*new_PEPSTensor(spin,bond,bond,integerONE,bond)
    enddo
    do n=2,YLength-1
        this%TensorCollection(1,n)=sqrtOfNorm*new_PEPSTensor(spin,integerONE,bond,bond,bond)
        this%TensorCollection(Xlength,n)=sqrtOfNorm*new_PEPSTensor(spin,bond,integerONE,bond,bond)
    enddo
    !Bulk terms are proper PEPS Tensors
    do m=2,YLength-1
        do n=2,Xlength-1
            this%TensorCollection(n,m)=sqrtOfNorm*new_PEPSTensor(spin,bond,bond,bond,bond)
        enddo
    enddo
    this%XLength=Xlength
    this%YLength=Ylength

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
     if (lhs%initialized) error=lhs%delete()
     allocate(lhs%TensorCollection(0:Xlength+1,0:Ylength+1))
     do m=0,Ylength+1
        do n=0,Xlength+1
            lhs%TensorCollection(n,m)=new_PEPSTensor(rhs%TensorCollection(n,m))
        enddo
     enddo
     lhs%XLength=Xlength
     lhs%YLength=Ylength

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

     deallocate(this%TensorCollection,stat=error)

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
    integer :: x,y
    type(PEPSTensor) :: tempTensor

    if(aPEPS%Initialized) then
        do y=1,aPEPS%Ylength
        do x=1,aPEPS%Xlength
            tempTensor=factor * (aPEPS%TensorCollection(x,y))
            call aPEPS%SetTensorAt(x,y, tempTensor )
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
        else
             call ThrowException('SetPEPSTensorAtSite','Site is wrong index',Xposition,CriticalError)
        endif
     else
         call ThrowException('SetPEPSTensorAtSite','PEPS or Tensor not initialized',NoErrorCode,CriticalError)
     endif
   end Subroutine SetPEPSTensorAtSite



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


   function GetPEPSRowAsMPS(aPEPS,whichRow,ColRange) result(aMPS)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: whichRow,ColRange(2)
     type(MPS) :: aMPS
     type(MPSTensor) :: tempTensor
     integer :: x

     if(aPEPS%Initialized) then
        if(IsPositionInsidePEPS(aPEPS,ColRange(1),whichRow) .and. IsPositionInsidePEPS(aPEPS,ColRange(2),whichRow) ) then
             aMPS=new_MPS(ColRange(2)-ColRange(1)+1)
             do x=ColRange(1),ColRange(2)
               tempTensor=aPEPS%TensorCollection(x,whichRow)%asMPSTensor(LEFT,RIGHT)
               call aMPS%SetTensorAt(x-ColRange(1)+1,tempTensor)
             enddo
        else
             call ThrowException('GetPEPSRowAsMPS','Site is wrong index',whichRow,CriticalError)
        endif
     else
         call ThrowException('GetPEPSRowAsMPS','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSRowAsMPS




   function GetPEPSColAsMPS(aPEPS,whichCol,RowRange) result(aMPS)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: whichCol,RowRange(2)
     type(MPS) :: aMPS
     type(MPSTensor) :: tempTensor
     integer :: y

     if(aPEPS%Initialized) then
        if(IsPositionInsidePEPS(aPEPS,whichCol,RowRange(1)) .and. IsPositionInsidePEPS(aPEPS,whichCol,RowRange(2)) ) then
             aMPS=new_MPS(RowRange(2)-RowRange(1)+1)
             do y=RowRange(1),RowRange(2)
               tempTensor=aPEPS%TensorCollection(whichCol,y)%asMPSTensor(DOWN,UP)
               call aMPS%SetTensorAt(y-RowRange(1)+1,tempTensor)
             enddo
        else
             call ThrowException('GetPEPSColAsMPS','Site is wrong index',whichCol,CriticalError)
        endif
     else
         call ThrowException('GetPEPSColAsMPS','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSColAsMPS

!!!!  TODO : String of PEPS as MPS, with corners and all.

   function GetCoreColAsMPS(aPEPS,whichCol,LeftMatrices,RightMatrices) result(aMPS)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: whichCol
     type(Tensor2),optional :: LeftMatrices(:),RightMatrices(:) !Make sure this vector is correct length
     type(MPS) :: aMPS
     type(PEPSTensor) :: tempTensor
     integer :: y

     if(aPEPS%Initialized) then
        if(IsPositionInsidePEPS(aPEPS,whichCol,1) ) then
             aMPS=new_MPS(aPEPS%Ylength)
             do y=1,aPEPS%YLength
               if (present(LeftMatrices) .and. present(RightMatrices) ) then
                  tempTensor=ApplyMatrixToBond(ApplyMatrixToBond(aPEPS%TensorCollection(whichCol,y),LeftMatrices(y),LEFT),RightMatrices(y),RIGHT)
               else
                  tempTensor=aPEPS%TensorCollection(whichCol,y)
               endif
               call aMPS%SetTensorAt(y,tempTensor%AsMPSTensor(DOWN,UP) )
             enddo
        else
             call ThrowException('GetCoreColAsMPS','Site is wrong index',whichCol,CriticalError)
        endif
     else
         call ThrowException('GetCoreColAsMPS','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetCoreColAsMPS


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   function GetCoreRowAsMPS(aPEPS,whichRow,BelowMatrices,AboveMatrices) result(aMPS)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: whichRow
     type(Tensor2),optional :: BelowMatrices(:),AboveMatrices(:) !Make sure this vector is correct length
     type(MPS) :: aMPS
     type(PEPSTensor) :: tempPEPS
     integer :: x

     if(aPEPS%Initialized) then
        if(IsPositionInsidePEPS(aPEPS,1,whichRow) ) then
             aMPS=new_MPS(aPEPS%Xlength)
             do x=1,aPEPS%XLength
               if (present(BelowMatrices) .and. present(AboveMatrices)) then
                  tempPEPS=ApplyMatrixToBond(ApplyMatrixToBond(aPEPS%TensorCollection(x,whichRow),BelowMatrices(x),DOWN),AboveMatrices(x),UP)
               else
                  tempPEPS=aPEPS%TensorCollection(x,whichRow)
               endif
               call aMPS%SetTensorAt(x,tempPEPS%AsMPSTensor(LEFT,RIGHT) )
             enddo
        else
             call ThrowException('GetCoreRowAsMPS','Site is wrong index',whichRow,CriticalError)
        endif
     else
         call ThrowException('GetCoreRowAsMPS','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetCoreRowAsMPS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



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
         spin = 2 !Magic thingy here, TODO
     else
         call ThrowException('GetPEPSSpin','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSSpin








   integer function GetMAXPEPSBond(aPEPS) result(bond)
     class(PEPS),intent(IN) :: aPEPS
     integer :: X,Y

     bond=0
     if(aPEPS%Initialized) then
       do Y=1,aPEPS%YLength
         do X=1,aPEPS%XLength
           ! call aPEPS%TensorCollection(X,Y)%PrintDimensions()
           bond= max( bond, maxval(aPEPS%TensorCollection(X,Y)%getBonds()) )
         enddo
       enddo
     else
         call ThrowException('GetPEPSBond','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetMaxPEPSBond







   integer function GetPEPSBond(aPEPS,siteX,siteY,aDirection) result(bond)
     class(PEPS),intent(IN) :: aPEPS
     integer,intent(IN) :: siteX,siteY,aDirection

     if(aPEPS%Initialized) then
        select case (aDirection)
            case (LEFT)
                bond=aPEPS%TensorCollection(siteX,siteY)%getDLeft()
            case (RIGHT)
                bond=aPEPS%TensorCollection(siteX,siteY)%getDRight()
            case (UP)
                bond=aPEPS%TensorCollection(siteX,siteY)%getDUp()
            case (DOWN)
                bond=aPEPS%TensorCollection(siteX,siteY)%getDDown()
            case default
                 call ThrowException('GetPEPSBond','direction must be LEFT/RIGHT/UP/DOWN',aDirection,CriticalError)
         end select
     else
         call ThrowException('GetPEPSBond','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function GetPEPSBond








!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




   function ReduceMAXPEPSBond(thisPEPS,newMaxBond) result(newPEPS)
     class(PEPS),intent(IN) :: thisPEPS
     integer,intent(IN) :: newMaxBond
     type(PEPS) :: newPEPS
     integer :: X,Y
     type(PEPSTensor) :: coreTensor,tempTensor
     type(Tensor2) :: UMatrices(RIGHT:UP)

     if(thisPEPS%Initialized) then
        newPEPS=thisPEPS
        if (newMaxBond.lt.thisPEPS%GetMaxBond()) then
            do Y=1,thisPEPS%YLength
                do X=1,thisPEPS%XLength
                    call newPEPS%TensorCollection(X,Y)%HOSVD(coreTensor,Umatrices,maxBondDim=[newMaxBond,newMaxBond],whichDirections=[RIGHT,UP])
                    !Umatrices are ordered as LEFT/RIGHT/UP/DOWN
!                    tempTensor=ApplyMatrixToBond(ApplyMatrixToBond(coreTensor,Umatrices(LEFT),LEFT),Umatrices(DOWN),DOWN)
                    call newPEPS%SetTensorAt(X,Y,CoreTensor)
                    !Now update the neighboring peps, careful with boundaries
                    if (X.lt.thisPEPS%XLength) then
                        tempTensor=ApplyMatrixToBond(newPEPS%TensorCollection(X+1,Y),TensorTranspose(Umatrices(RIGHT)),LEFT)
                        call newPEPS%SetTensorAt(X+1,Y,tempTensor)
                    endif
                    if (Y.lt.thisPEPS%YLength) then
                        tempTensor=ApplyMatrixToBond(newPEPS%TensorCollection(X,Y+1),TensorTranspose(Umatrices(UP)),DOWN)
                        call newPEPS%SetTensorAt(X,Y+1,tempTensor)
                    endif
                enddo
            enddo
            if (.not.newPEPS%IsPEPSWellFormed()) then
              call ThrowException('SetMAXPEPSBond','Bonds do not match -- check for errors :)',NoErrorCode,CriticalError)
            endif
         endif
     else
         call ThrowException('SetMAXPEPSBond','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function ReduceMaxPEPSBond


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  subroutine CanonizePEPSAtSite(aPEPS,siteX,siteY,direction,MaxLongitudinalBond,MaxTransverseBond,returnTheNorm)
   class(PEPS),intent(INOUT) :: aPEPS
   integer,intent(IN) :: siteX,siteY
   integer,intent(IN) :: direction
   real(8),intent(OUT),optional :: returnTheNorm
   integer :: MaxLongitudinalBond,MaxTransverseBond
   integer :: x,y
   real(8) :: NormOfTensor
   type(PEPSTensor) :: tempTensor

     if(aPEPS%Initialized) then

      if (direction.eq.HORIZONTAL) then
         do x=1,siteX-1
            if (verbose) print *,'Canonizing PEPS to right, col:',x
            call DecouplePEPSVerticallyAtCol(aPEPS,x,RIGHT,MaxLongitudinalBond,MaxTransverseBond)
         enddo
         do x=aPEPS%XLength,siteX+1,-1
            if (verbose) print *,'Canonizing PEPS to left, col:',x
            call DecouplePEPSVerticallyAtCol(aPEPS,x,LEFT,MaxLongitudinalBond,MaxTransverseBond)
         enddo
      !   call CanonizeCoreColumnAtSite(aPEPS,siteX,siteY,MaxTransverseBond**2)
         do y=1,siteY-1
            if (verbose) print *,'Canonizing PEPS up, row:',y
            call CanonizeCoreSiteAndPushNorm(aPEPS,siteX,y,UP,MaxTransverseBond**2)
         enddo
         do y=aPEPS%YLength,siteY+1,-1
            if (verbose) print *,'Canonizing PEPS down, row:',y
            call CanonizeCoreSiteAndPushNorm(aPEPS,siteX,y,DOWN,MaxTransverseBond**2)
         enddo
       else if (direction.eq.VERTICAL) then
         do y=1,siteY-1
            call DecouplePEPSHorizontallyAtRow(aPEPS,y,UP,MaxLongitudinalBond,MaxTransverseBond)
         enddo
         do y=aPEPS%YLength,siteY+1,-1
            call DecouplePEPSHorizontallyAtRow(aPEPS,y,DOWN,MaxLongitudinalBond,MaxTransverseBond)
         enddo
     !    call CanonizeCoreRowAtSite(aPEPS,siteX,siteY,MaxTransverseBond**2)
         do x=1,siteX-1
            call CanonizeCoreSiteAndPushNorm(aPEPS,x,siteY,RIGHT,MaxTransverseBond**2)
         enddo
         do x=aPEPS%XLength,siteX+1,-1
            call CanonizeCoreSiteAndPushNorm(aPEPS,x,siteY,LEFT,MaxTransverseBond**2)
         enddo
       else
         call ThrowException('CanonizePEPSAtSite','Direction must be HORIZONTAL or VERTICAL',NoErrorCode,CriticalError)
       endif
      !Finally normalize the centered tensor so that norm of PEPS is ONE
       NormOfTensor=aPEPS%TensorCollection(siteX,siteY)%Norm()
       if(present(returnTheNorm)) returnTheNorm=NormOfTensor
!       tempTensor=(ONE/sqrt(NormOfTensor)) * (aPEPS%TensorCollection(siteX,siteY))
!       call aPEPS%SetTensorAt(siteX,siteY, tempTensor )
!       if (debug) print *,' -x- PEPS Core site Normalized by:',NormOfTensor
     else
         call ThrowException('CanonizePEPSAtSite','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end subroutine CanonizePEPSAtSite


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine DecouplePEPSHorizontallyAtRow(aPEPS,WhichRow,directionToPush,MaxLongitudinalBond,MaxTransverseBond)
      class(PEPS),intent(INOUT) :: aPEPS
      integer,intent(IN) :: WhichRow,directionToPush,MaxLongitudinalBond,MaxTransverseBond
      integer :: x

      do x=1,aPEPS%XLength
          call FullCanonizePEPSSite(aPEPS,x,WhichRow,directionToPush,MaxLongitudinalBond,MaxTransverseBond)
      enddo

   end subroutine DecouplePEPSHorizontallyAtRow

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine DecouplePEPSVerticallyAtCol(aPEPS,WhichCol,directionToPush,MaxLongitudinalBond,MaxTransverseBond)
      class(PEPS),intent(INOUT) :: aPEPS
      integer,intent(IN) :: WhichCol,directionToPush,MaxLongitudinalBond,MaxTransverseBond
      integer :: y

      do y=1,aPEPS%YLength
          call FullCanonizePEPSSite(aPEPS,WhichCol,y,directionToPush,MaxLongitudinalBond,MaxTransverseBond)
      enddo

   end subroutine DecouplePEPSVerticallyAtCol


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine FullCanonizePEPSSite(aPEPS,Xcol,Yrow,directionToPush,MaxLongitudinalBond,MaxTransverseBond)
      class(PEPS),intent(INOUT) :: aPEPS
      integer,intent(IN) :: XCol,Yrow,directionToPush,MaxLongitudinalBond,MaxTransverseBond
      integer :: y,spin
      type(PEPSTensor) :: tempPEPSTensor
      type(Tensor4) :: tempTensor4,pushTensor
      type(Tensor3) :: tempTensor3
      type(Tensor2) :: UMatrices(LEFT:DOWN)
      integer :: BondSizes(LEFT:DOWN),BondLimits(LEFT:DOWN)
      integer :: deltaX,deltaY,TransversePush,TransverseDoNotPush

      select case (directionToPush)
         case (RIGHT)
            deltaX=1
            deltaY=0
            TransversePush=UP
            TransverseDoNotPush=DOWN
         case (LEFT)
            deltaX=-1
            deltaY=0
            TransversePush=UP
            TransverseDoNotPush=DOWN
         case (UP)
            deltaX=0
            deltaY=1
            TransversePush=RIGHT
            TransverseDoNotPush=LEFT
         case (DOWN)
            deltaX=0
            deltaY=-1
            TransversePush=RIGHT
            TransverseDoNotPush=LEFT
      end select

      !First join spin dimension of tensor with the leg opposite the push direction
      tempTensor4=aPEPS%TensorCollection(Xcol,Yrow)%JoinSpinWith(DirectionOppositeTo(directionToPush))
      !Need to prepare the truncation of bonds
      BondSizes=aPEPS%TensorCollection(Xcol,Yrow)%GetBonds()
      spin=aPEPS%TensorCollection(Xcol,Yrow)%GetSpin()
      !keep track of growing bonds near borders
      BondLimits(directionToPush)=BondSizes(directionToPush)
      BondLimits(DirectionOppositeTo(directionToPush))= min(MaxLongitudinalBond, spin*BondSizes(DirectionOppositeTo(directionToPush)) )
      !Transverse dir gets special treatment
      BondLimits(TransverseDoNotPush)=BondSizes(TransverseDoNotPush)
      BondLimits(TransversePush)=MaxTransverseBond

      call tempTensor4%SVD(pushTensor,Umatrices,BondLimits)

      !Recover the do-not push matrix into the push tensor
      pushTensor=nModeProduct(Umatrices(TransverseDoNotPush),pushTensor,[TransverseDoNotPush])
      !And now push the transverse matrix
      if ( IsPositionInsidePEPS(aPEPS,Xcol+abs(deltaY),Yrow+abs(deltaX)) ) then
        tempPEPSTensor=ApplyMatrixToBond(aPEPS%TensorCollection(Xcol+abs(deltaY),Yrow+abs(deltaX)), &
           & TensorTranspose(Umatrices(TransversePush)),directionOppositeTo(TransversePush))
        call aPEPS%SetTensorAt(Xcol+abs(deltaY),Yrow+abs(deltaX),tempPEPSTensor)
      endif

      !and now push longitudinal matrix first, tensor comes later (this reduces cost because matrix is truncated)
      if (IsPositionInsidePEPS(aPEPS,Xcol+deltaX,Yrow+deltaY)) then
         tempPEPSTensor=ApplyMatrixToBond(aPEPS%TensorCollection(Xcol+deltaX,Yrow+deltaY), &
            & TensorTranspose(Umatrices(directionToPush)),DirectionOppositeTo(directionToPush))
         tempPEPSTensor=ApplyMPOToBond(tempPEPSTensor,pushTensor,DirectionOppositeTo(directionToPush))
         call aPEPS%SetTensorAt(Xcol+deltaX,Yrow+deltaY,tempPEPSTensor)
      endif

      !Finally reconstruct a new PEPSTensor on site from the remaining Umatrix
      select case (directionToPush)
         case (LEFT)
            tempTensor3 = TensorReshape( TensorTranspose(Umatrices(RIGHT)), & !Its transposed because second index always corresponds to core tensor
               &  [ BondLimits(RIGHT), BondSizes(RIGHT), spin ] )
            tempPEPSTensor=new_PEPSTensor(tempTensor3, 3, HORIZONTAL)
         case (RIGHT)
            tempTensor3 = TensorReshape( Umatrices(LEFT), &
               &  [ BondSizes(LEFT), spin, BondLimits(LEFT) ] )
            tempPEPSTensor=new_PEPSTensor(tempTensor3, 2, HORIZONTAL)
         case (UP)
            tempTensor3 = TensorReshape( TensorTranspose(Umatrices(DOWN)), & !Second index of U is always the one corresponding to core tensor
               &  [ BondLimits(DOWN), BondSizes(DOWN), spin ] )
            tempPEPSTensor=new_PEPSTensor(tempTensor3, 3, VERTICAL)
         case (DOWN)
            tempTensor3 = TensorReshape( Umatrices(UP), &
               &  [ BondSizes(UP), spin, BondLimits(UP) ] )
            tempPEPSTensor=new_PEPSTensor(tempTensor3, 2, VERTICAL)
      end select

      call aPEPS%SetTensorAt(Xcol,Yrow,tempPEPSTensor)

  end subroutine FullCanonizePEPSSite

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine CanonizeCoreSiteAndPushNorm(aPEPS,siteX,siteY,directionToPush,MaxBond)
      class(PEPS),intent(INOUT) :: aPEPS
      integer,intent(IN) :: siteX,siteY,directionToPush,MaxBond
      type(PEPSTensor) :: coreTensor,tempTensor
      type(Tensor2) :: UMatrices(1),NormMatrix(1),SqrtOfNormMatrix,InverseSqrtNorm
      integer :: deltaR(2),RankTruncation(1),direction

      RankTruncation=MaxBond
      call aPEPS%TensorCollection(siteX,siteY)%HOSVD(coreTensor,Umatrices,RankTruncation,[directionToPush],SqrtOfNormMatrix)
      !NormMatrix(1)=TakeDiagonalPart(coreTensor%CollapseAllIndicesBut(directionToPush))
      !SqrtOfNormMatrix=TensorSqrt(NormMatrix(1))

      !Finally push the U and then the norm matrix to the neighbor tensor
      select case (directionToPush)
         case (UP)
            deltaR=[0,1]
         case (DOWN)
            deltaR=[0,-1]
         case (LEFT)
            deltaR=[-1,0]
         case (RIGHT)
            deltaR=[1,0]
      end select
      if (IsPositionInsidePEPS(aPEPS,siteX+deltaR(1),siteY+deltaR(2))) then
         !renormalize the tensor with the inverse sqrt of the norm
         coreTensor=ApplyMatrixToBond(coreTensor,PseudoInverseDiagonal(SqrtOfNormMatrix),directionToPush)
         call aPEPS%SetTensorAt(siteX,siteY,coreTensor)
         !And now push the U matrix and then the sqrt of the norm
         tempTensor=aPEPS%TensorCollection(siteX+deltaR(1),siteY+deltaR(2))
         tempTensor=ApplyMatrixToBond(tempTensor,TensorTranspose(Umatrices(1)),DirectionOppositeTo(directionToPush))
         tempTensor=ApplyMatrixToBond(tempTensor,SqrtOfNormMatrix,DirectionOppositeTo(directionToPush))
         call aPEPS%SetTensorAt(siteX+deltaR(1),siteY+deltaR(2),tempTensor)
      else
         !Store the tensor but without adding the norm
         coreTensor=ApplyMatrixToBond(coreTensor,Umatrices(1),directionToPush)
         call aPEPS%SetTensorAt(siteX,siteY,coreTensor)
      endif

   end subroutine CanonizeCoreSiteAndPushNorm

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine CanonizeCoreColumnAtSite(aPEPS,whichCol,siteY,MaxBond)
      class(PEPS),intent(INOUT) :: aPEPS
      integer,intent(IN) :: siteY,whichCol,MaxBond
      type(PEPSTensor) :: tempTensor
      type(MPSTensor) :: tempMPS
      type(MPS) :: CoreMPS,reducedMPS
      integer :: y,dims(4)
      real(8) :: Norma

      !First extract the core as an MPS
      if (Debug) print *,'Extracting Core...'
      CoreMPS=GetCoreColAsMPS(aPEPS,whichCol)
      if (Debug) print *,'Core extracted, normalizing it'
      call Normalize_MPS(CoreMPS)
      if (Debug) print *,'Core normalized, reducing it'
      reducedMPS=Approximate(CoreMPS,MaxBond)
!      if (Debug) print *,'Core reduced, canonizing it'
!      call CoreMPS%Canonize()
!      if (Debug) print *,'Core canonized, putting it in normal form at requested site...'
!      do y=1,siteY-1
!         call reducedMPS%RightCanonizeAtSite(y)
!      enddo
!   !Might need this
!      !Finally normalize the center tensor
!      tempMPS=CoreMPS%GetTensorAt(siteX)
!      Norma=Norm(tempMPS)
!      tempMPS=(ONE/sqrt(Norma)) * tempMPS
!      call CoreMPS%SetTensorAt(siteX,tempMPS)
      !Finally, rewrite the MPS as a core column
      if (Debug) print *,'Restoring MPS to core in PEPS'
      do y=1,aPEPS%YLength
         if (Debug) print *,'    --- row:',y
         dims=aPEPS%TensorCollection(whichCol,y)%GetBonds()
         tempMPS=reducedMPS%GetTensorAt(y)
         tempTensor=new_PEPSTensor(tempMPS,3,VERTICAL,[dims(1),dims(2)])
         call aPEPS%SetTensorAt(whichCol,y,tempTensor)
      enddo

!      !Canonize in a unorthodox way, first left part
!      do y=siteY-1,1-1
!         call CoreMPS%LeftCanonizeAtSite(y)
!      enddo
!      call CoreMPS%SetTensorAt(0,new_MPSTensor(integerONE,integerONE,integerONE,ONE))
!      do y=1,siteY-1
!         call CoreMPS%RightCanonizeAtSite(y)
!      enddo
!      !And now right part
!      do y=siteY+1,aPEPS%Ylength
!         call CoreMPS%RightCanonizeAtSite(y)
!      enddo
!      call CoreMPS%SetTensorAt(aPEPS%Ylength,new_MPSTensor(integerONE,integerONE,integerONE,ONE))
!      do y=aPEPS%Ylength,siteY+1,-1
!         call CoreMPS%LeftCanonizeAtSite(y)
!      enddo
!      !Finally normalize the center tensor
!      tempMPS=CoreMPS%GetTensorAt(siteY)
!      Norma=Norm(tempMPS)
!      tempMPS=(ONE/sqrt(Norma)) * tempMPS
!      call CoreMPS%SetTensorAt(siteY,tempMPS)
!      !Finally, rewrite the MPS as a core column
!      do y=1,aPEPS%YLength
!         dims=aPEPS%TensorCollection(whichCol,y)%GetBonds()
!         tempTensor=new_PEPSTensor(CoreMPS%GetTensorAt(y),3,VERTICAL,[dims(1),dims(2)])
!         call aPEPS%SetTensorAt(whichCol,y,tempTensor)
!      enddo

   end subroutine CanonizeCoreColumnAtSite

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine CanonizeCoreRowAtSite(aPEPS,siteX,whichRow,MaxBond)
      class(PEPS),intent(INOUT) :: aPEPS
      integer,intent(IN) :: siteX,whichRow,MaxBond
      type(PEPSTensor) :: tempTensor
      type(MPSTensor) :: tempMPS
      type(MPS) :: CoreMPS,reducedMPS
      integer :: x,dims(4)
      real(8) :: Norma

      if (Debug) print *,'Extracting Row Core...'
      !First extract the core as an MPS
      CoreMPS=GetCoreRowAsMPS(aPEPS,whichRow)
      if (Debug) print *,'Core extracted, normalizing it'
      call Normalize_MPS(CoreMPS)
      if (Debug) print *,'Core normalized, reducing it'
      reducedMPS=Approximate(CoreMPS,MaxBond)
      if (Debug) print *,'Core normalized, canonizing it'
!      call CoreMPS%Canonize()
!      do x=1,siteX-1
!         call reducedMPS%RightCanonizeAtSite(x)
!      enddo
!   !Might need this
!      !Finally normalize the center tensor
!      tempMPS=CoreMPS%GetTensorAt(siteX)
!      Norma=Norm(tempMPS)
!      tempMPS=(ONE/sqrt(Norma)) * tempMPS
!      call CoreMPS%SetTensorAt(siteX,tempMPS)
      !Finally, rewrite the MPS as a core column
      do x=1,aPEPS%XLength
         dims=aPEPS%TensorCollection(x,whichRow)%GetBonds()
         tempMPS=reducedMPS%GetTensorAt(x)
         tempTensor=new_PEPSTensor(tempMPS,3,HORIZONTAL,[dims(3),dims(4)])
         call aPEPS%SetTensorAt(x,whichRow,tempTensor)
      enddo

   end subroutine CanonizeCoreRowAtSite

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    logical function Integrity_Check_of_bonds(thisPEPS) result(isAllOK)
     class(PEPS),intent(IN) :: thisPEPS
     integer :: X,Y
     logical :: TempCheck

     isAllOK=.true.
     if(thisPEPS%Initialized) then
         do Y=1,thisPEPS%YLength
            do X=1,thisPEPS%XLength
               isAllOK= isAllOK.and.thisPEPS%TensorCollection(X,Y)%IsInitialized()
               if (verbose.and. .not.thisPEPS%TensorCollection(X,Y)%IsInitialized()) then
                   call ThrowException('Integrity_Check_of_bonds','Tensor not initialized',(Y-1)*thisPEPS%XLength+X,Warning)
                   print *,'Coordinate for previous error is:',X,Y
               endif
               TempCheck=thisPEPS%TensorCollection(X,Y)%GetDLeft().eq.thisPEPS%TensorCollection(X-1,Y)%GetDRight()
               isAllOK= isAllOK.and.TempCheck
               if (verbose.and. .not.TempCheck) then
                   call ThrowException('Integrity_Check_of_bonds','Left bond does not match',(Y-1)*thisPEPS%XLength+X,Warning)
                   print *,'Coordinate for previous error is:',X,Y
                   print *,'Bonds are: ',thisPEPS%TensorCollection(X-1,Y)%GetDRight(),thisPEPS%TensorCollection(X,Y)%GetDLeft()
               endif
               TempCheck=thisPEPS%TensorCollection(X,Y)%GetDRight().eq.thisPEPS%TensorCollection(X+1,Y)%GetDLeft()
               isAllOK= isAllOK.and.TempCheck
               if (verbose.and. .not.TempCheck) then
                   call ThrowException('Integrity_Check_of_bonds','Right bond does not match',(Y-1)*thisPEPS%XLength+X,Warning)
                   print *,'Coordinate for previous error is:',X,Y
                   print *,'Bonds are: ',thisPEPS%TensorCollection(X,Y)%GetDRight(),thisPEPS%TensorCollection(X+1,Y)%GetDLeft()
               endif
               TempCheck=thisPEPS%TensorCollection(X,Y)%GetDUp().eq.thisPEPS%TensorCollection(X,Y+1)%GetDDown()
               isAllOK= isAllOK.and.TempCheck
               if (verbose.and. .not.TempCheck) then
                   call ThrowException('Integrity_Check_of_bonds','Up bond does not match',(Y-1)*thisPEPS%XLength+X,Warning)
                   print *,'Coordinate for previous error is:',X,Y
                   print *,'Bonds are: ',thisPEPS%TensorCollection(X,Y)%GetDUp(),thisPEPS%TensorCollection(X,Y+1)%GetDDown()
               endif
               TempCheck=thisPEPS%TensorCollection(X,Y)%GetDDown().eq.thisPEPS%TensorCollection(X,Y-1)%GetDUp()
               isAllOK= isAllOK.and.TempCheck
               if (verbose.and. .not.TempCheck) then
                   call ThrowException('Integrity_Check_of_bonds','Down bond does not match',(Y-1)*thisPEPS%XLength+X,Warning)
                   print *,'Coordinate for previous error is:',X,Y
                   print *,'Bonds are: ',thisPEPS%TensorCollection(X,Y)%GetDDown(),thisPEPS%TensorCollection(X,Y-1)%GetDUp()
               endif
            enddo
         enddo
     else
         call ThrowException('Integrity_Check_of_bonds','PEPS not initialized',NoErrorCode,CriticalError)
     endif
   end function Integrity_Check_of_bonds


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine PrintPEPSBondDimensionsMap(aPEPS,message)
    class(PEPS),intent(IN) :: aPEPS
    character*(*),optional :: message
    integer :: x,y,length
    integer,allocatable :: dims(:)
!    character(LEN=20),parameter :: upaboveFormat='(A,I2)'

    if (present(message)) then
      print *,message
    endif
    print *,'PEPS dims:'
    print *,'========='

    length=2*aPEPS%XLength
    allocate (dims(length))

    do y=aPEPS%Ylength,1,-1
        do x=1,aPEPS%XLength
            dims(x)=aPEPS%TensorCollection(x,y)%GetDUp()
        enddo
        write(*,'(2x,<length>(I2,6x))'), (dims(x),x=1,aPEPS%XLength)
        do x=1,aPEPS%XLength
            dims(2*x-1)=aPEPS%TensorCollection(x,y)%GetDLeft()
            dims(2*x)=aPEPS%TensorCollection(x,y)%GetDRight()
        enddo
        write(*,'(<length>(I2," H",I2," -"))') (dims(x),x=1,2*aPEPS%XLength)
        !'(I2,"100x ",I2,"-")'), dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)
        do x=1,aPEPS%XLength
            dims(x)=aPEPS%TensorCollection(x,y)%GetDDown()
        enddo
        write(*,'(2x,<length>(I2,6x))'), (dims(x),x=1,aPEPS%XLength)
        write(*,'(2x,<length>(A2,6x))'), ('|',x=1,aPEPS%XLength)
    enddo

end subroutine PrintPEPSBondDimensionsMap

logical function IsPositionInsidePEPS(aPEPS,siteX,siteY)
   class(PEPS),intent(IN) :: aPEPS
   integer,intent(IN) :: siteX,siteY

   IsPositionInsidePEPS=siteX.ge.1
   IsPositionInsidePEPS=IsPositionInsidePEPS.and. (siteX.le.aPEPS%XLength)
   IsPositionInsidePEPS=IsPositionInsidePEPS.and. (siteY.ge.1)
   IsPositionInsidePEPS=IsPositionInsidePEPS.and. (siteY.le.aPEPS%YLength)

end function IsPositionInsidePEPS
 end module PEPS_Class


