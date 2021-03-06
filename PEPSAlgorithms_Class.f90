
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

module PEPSAlgorithms_Class

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
  use Multiplicator2D_Class

  implicit none

    real(8) :: PSEUDOINVERSETOLERANCE = 1.0d-8

    interface Approximate
        module procedure Approximate_PEPS
    end interface

    interface Normalize
        module procedure Normalize_PEPS
    end interface

    interface Overlap
        module procedure Overlap_PEPS, Canonical_PEPS_Overlap
    end interface

    interface ExpectationValue
        module procedure Expectation_Value_PEPS_PEPO, Expectation_Value_Canonical_PEPS_MPO
    end interface

  contains

!*****************************************************************************
  function Approximate_PEPS(bigPEPS, newBondDimension, returnOverlap) result(smallPEPS)
    type(PEPS), intent(INOUT) :: bigPEPS
    integer,intent(IN) :: newBondDimension
    real(8),optional,intent(OUT) :: returnOverlap
    type(PEPS) :: smallPEPS
    type(PEPSTensor) :: localTensor
    type(Multiplicator2D) :: smallMultiplicator,bigMultiplicator
    integer :: PEPSSize(2),localPEPSDims(5),PEPSSpin
    real(8) :: overlap
    integer :: siteX,siteY,sweep
    type(Tensor4) :: aTempTensor
    logical,allocatable,target :: HasPEPSChangedAt(:,:)

    if(bigPEPS%IsInitialized()) then

        PEPSsize=bigPEPS%GetSize()
        PEPSSpin=bigPEPS%GetSpin()
        if (newBondDimension.ge.bigPEPS%GetMaxBond() ) then
            smallPEPS = new_PEPS(bigPEPS)
            return
        else
            print *,'inside AP'
            call Normalize(bigPEPS)
            print *,'normalized B_PEPS'
            smallPEPS = ReduceMAXPEPSBond(bigPEPS,newBondDimension) !new_PEPS(PEPSSize(1),PEPSSize(2),PEPSSpin,newBondDimension) !
            print *,'reduced bond'
            call Normalize(smallPEPS)
            print *,'normalized S_PEPS'
            !if( smallPEPS
        endif

        allocate(HasPEPSChangedAt(PEPSSize(1),PEPSSize(2)))
        HasPEPSChangedAt=.true.

        smallMultiplicator = new_Multiplicator2D(smallPEPS,MatrixToTrackChanges=HasPEPSChangedAt)
        bigMultiplicator = new_Multiplicator2D(bigPEPS,smallPEPS,MatrixToTrackChanges=HasPEPSChangedAt)

        print *,'Max sweeping:',MaxSweeps
        sweep=0
        !Start sweeping
        !overlap=Abs(bigMultiplicator%FullContraction())**2
        do while (sweep .le. MaxSweeps ) !.and. (1.d0-overlap).gt.ApproximationTolerance )
            sweep=sweep+1
            !Sweep from 1,1 first to the right and then up
            do siteY=1,PEPSSize(2)
	            do siteX=1,PEPSSize(1)
	                print *,'Optimizing site :',siteX,siteY
	                localTensor = smallPEPS%GetTensorAt(siteX,siteY)
	                localPEPSDims = localTensor%GetDimensions()
                    localTensor = ComputePEPSOptimum ( smallMultiplicator%LeftAt(siteX-1,siteY), smallMultiplicator%RightAt(siteX+1,siteY), &
                        & smallMultiplicator%AboveAt(siteX,siteY+1), smallMultiplicator%BelowAt(siteX,siteY-1), &
                        & bigMultiplicator%LeftAt(siteX-1,siteY), bigMultiplicator%RightAt(siteX+1,siteY), &
                        & bigMultiplicator%AboveAt(siteX,siteY+1), bigMultiplicator%BelowAt(siteX,siteY-1), &
                        & bigPEPS%GetTensorAt(siteX,siteY), localPEPSDims )
                        call smallPEPS%SetTensorAt(siteX,siteY,localTensor)  ; HasPEPSChangedAt(siteX,siteY)=.true.
                enddo
                call smallMultiplicator%Reset(RIGHT)
                call bigMultiplicator%Reset(RIGHT)
            enddo
            !At the end of the sweep, the  multiplicators are reset
            call smallMultiplicator%Reset(UP)
            call bigMultiplicator%Reset(UP)
	        !Sweep in opposite direction
            do siteY=PEPSSize(2),1,-1
                do siteX=PEPSSize(1),1,-1
                    localTensor = smallPEPS%GetTensorAt(siteX,siteY)
                    localPEPSDims = localTensor%GetDimensions()
                    localTensor = ComputePEPSOptimum ( smallMultiplicator%LeftAt(siteX-1,siteY), smallMultiplicator%RightAt(siteX+1,siteY), &
                        & smallMultiplicator%AboveAt(siteX,siteY+1), smallMultiplicator%BelowAt(siteX,siteY-1), &
                        & bigMultiplicator%LeftAt(siteX-1,siteY), bigMultiplicator%RightAt(siteX+1,siteY), &
                        & bigMultiplicator%AboveAt(siteX,siteY+1), bigMultiplicator%BelowAt(siteX,siteY-1), &
                        & bigPEPS%GetTensorAt(siteX,siteY) , localPEPSDims)
                        call smallPEPS%SetTensorAt(siteX,siteY,localTensor) ; HasPEPSChangedAt(siteX,siteY)=.true.
                enddo
                call smallMultiplicator%Reset(LEFT)
                call bigMultiplicator%Reset(LEFT)
            enddo
            !At the end of the sweep, the multiplicators are reset
            call smallMultiplicator%Reset(DOWN)
            call bigMultiplicator%Reset(DOWN)
            !overlap=Abs(bigMultiplicator%FullContraction())**2
        enddo
        call smallMultiplicator%Delete()
        call bigMultiplicator%Delete()
        if (present(returnOverlap)) then
            returnOverlap=overlap
        endif
        If (Verbose) print *,'PEPS-Approximate/Total sweeps performed: ',sweep
    else
        call ThrowException('PEPS-approximate algorithm','PEPS not initialized',NoErrorCode,CriticalError)
    endif

  end function Approximate_PEPS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function ComputePEPSOptimum ( smallE_Left, smallE_Right, smallE_Up, smallE_Down, &
                             &  bigE_Left, bigE_Right, bigE_Up, bigE_Down, aPEPSTensor, newDims) result(newTensor)
    type(Tensor4),intent(IN) :: smallE_Left, smallE_Right, smallE_Up, smallE_Down
    type(Tensor4),intent(IN) :: bigE_Left, bigE_Right, bigE_Up, bigE_Down
    type(PEPSTensor),intent(IN) :: aPEPSTensor
    integer,intent(IN) :: newDims(5)
    type(PEPSTensor) :: newTensor
    type(Tensor2) :: smallMatrix,bigMatrixTimesVector,aTempMatrix

    print *,'About to contract 4 environments'
    smallMatrix=TensorTranspose( TensorTrace(  &
                & ((smallE_left.xplus.smallE_Up).xplus.smallE_Right).xplus. TensorTranspose(smallE_Down, [2,1,3,4] ), &
                &                           THIRDANDFOURTH )      )

    call smallMatrix%PrintDimensions('Small environment dims:')

    bigMatrixTimesVector=TensorTranspose( TensorTrace( &
                & ((bigE_left.xplus.bigE_Up).xplus.bigE_Right).xplus. TensorTranspose( bigE_Down, [2,1,3,4] ), &
                                            THIRDANDFOURTH )      ) .x. aPEPSTensor%CompactBonds()

    print *,'Environments computed'

!    call bigMatrixTimesVector%PrintDimensions('Big environment dims:')

!    call smallMatrix%Print('Small matrix data:')
!    call bigMatrixTimesVector%Print('Big matrix*vector data:')

    newTensor= new_PEPSTensor(  &
       & SplitIndexOf(   &
       & SolveLinearProblem(smallMatrix, bigMatrixTimesVector, PSEUDOINVERSETOLERANCE ), newDims )  &
       &                      )
!        call newTensor%Print('New Tensor data:')
  end function ComputePEPSOptimum

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function Overlap_PEPS(onePEPS, anotherPEPS) result(theOverlap)
      class(PEPS),intent(IN) :: onePEPS
      class(PEPS),intent(IN),optional :: anotherPEPS
      type(Multiplicator2D) :: theEnvironment
      complex(8) :: theOverlap
      integer :: dims(2)
      logical,allocatable,target :: HasPEPSChangedAt(:,:)

      dims=onePEPS%GetSize()
      allocate(HasPEPSChangedAt(dims(1),dims(2)))
      HasPEPSChangedAt=.true.
      if (present(anotherPEPS)) then
          theEnvironment=new_Multiplicator2D(onePEPS, anotherPEPS,MatrixToTrackChanges=HasPEPSChangedAt)
      else
          theEnvironment=new_Multiplicator2D(onePEPS,MatrixToTrackChanges=HasPEPSChangedAt)
      endif
      theOverlap = Overlap_PEPSAboveBelow(theEnvironment)
      call theEnvironment%Delete()

  end function Overlap_PEPS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function Expectation_Value_PEPS_PEPO(onePEPS, aPEPO, anotherPEPS) result(theOverlap)
      class(PEPS),intent(IN) :: onePEPS
      class(PEPO),intent(IN) :: aPEPO
      class(PEPS),intent(IN),optional :: anotherPEPS
      type(Multiplicator2D) :: theEnvironment
      complex(8) :: theOverlap
      integer :: dims(2)
      logical,allocatable,target :: HasPEPSChangedAt(:,:)

      dims=onePEPS%GetSize()
      allocate(HasPEPSChangedAt(dims(1),dims(2)))
      HasPEPSChangedAt=.true.
      if (present(anotherPEPS)) then
          theEnvironment=new_Multiplicator2D(onePEPS, anotherPEPS,aPEPO,MatrixToTrackChanges=HasPEPSChangedAt)
      else
          theEnvironment=new_Multiplicator2D(onePEPS, PEPO_C=aPEPO, MatrixToTrackChanges=HasPEPSChangedAt)
      endif
      theOverlap = Overlap_PEPSAboveBelow(theEnvironment)
      call theEnvironment%Delete()

  end function Expectation_Value_PEPS_PEPO

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine Normalize_PEPS(aPEPS)
      class(PEPS),intent(INOUT) :: aPEPS
      real(8) :: theNorm
      integer :: TotalNumberOfTensors

      theNorm = abs(Overlap_PEPS(aPEPS))   !!!Notice I should not use **2
      TotalNumberOfTensors=product(aPEPS%GetSize())
      call aPEPS%ScaleBy(ONE/(theNorm**(0.5d0/TotalNumberOfTensors)))

  end subroutine Normalize_PEPS




!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  function Canonical_PEPS_Overlap(onePEPS, anotherPEPS, CorePosition, CoreDirection) result(theOverlap)
      class(PEPS),intent(IN),target :: onePEPS
      class(PEPS),intent(IN),optional,target :: anotherPEPS
      type(PEPS),pointer :: belowPEPS
      integer,intent(IN) :: CorePosition, CoreDirection
      complex(8) :: theOverlap
      type(Tensor2),allocatable :: RightMatrices(:),LeftMatrices(:)
      type(Multiplicator) :: aMultiplicator
      type(MPS) :: MPSAbove,MPSBelow,CoreMPSBelow,CoreMPSAbove
      integer :: x,y,dims(2)
      integer :: iTemp(2)

      if (present(anotherPEPS)) then
         belowPEPS => anotherPEPS
      else
         belowPEPS => onePEPS
      endif

      dims=onePEPS%GetSize()
      if (dims.equalvector.belowPEPS%GetSize()) then
        select case (CoreDirection)
         case (VERTICAL)
            allocate (RightMatrices(dims(2)),LeftMatrices(dims(2)))
            do y=1,dims(2)
               if (CorePosition.gt.1) then
                  MPSAbove=onePEPS%GetRowAsMPS(y,[1,CorePosition-1])
                  MPSBelow=belowPEPS%GetRowAsMPS(y,[1,CorePosition-1])
                  aMultiplicator=new_Multiplicator(MPSAbove,MPSBelow)
                  LeftMatrices(y)=LeftAtSite(aMultiplicator,CorePosition-1)
                  !iTemp=LeftMatrices(y)%GetDimensions()
               else
                  LeftMatrices(y)=Identity(1)
               endif
               if (CorePosition.lt.dims(1)) then
                  MPSAbove=onePEPS%GetRowAsMPS(y,[CorePosition+1,dims(1)])
                  MPSBelow=belowPEPS%GetRowAsMPS(y,[CorePosition+1,dims(1)])
                  aMultiplicator=new_Multiplicator(MPSAbove,MPSBelow)
                  RightMatrices(y)=RightAtSite(aMultiplicator,1)
               else
                  RightMatrices(y)=Identity(1)
               endif
            enddo
            CoreMPSBelow=GetColAsMPS(belowPEPS,CorePosition)
            CoreMPSAbove=GetColAsMPS(onePEPS,CorePosition,LeftMatrices,RightMatrices)
         case (HORIZONTAL)
            allocate (RightMatrices(dims(1)),LeftMatrices(dims(1)))
            do x=1,dims(1)
               if (CorePosition.gt.1) then
                  MPSAbove=onePEPS%GetColAsMPS(x,[1,CorePosition-1])
                  MPSBelow=belowPEPS%GetColAsMPS(x,[1,CorePosition-1])
                  aMultiplicator=new_Multiplicator(MPSAbove,MPSBelow)
                  LeftMatrices(x)=LeftAtSite(aMultiplicator,CorePosition-1)
                  if (Debug) Print *,'Left matrices x=',x,Norm(LeftMatrices(x)-Identity(onePEPS%GetBondAt(x,CorePosition,DOWN)))
               else
                  LeftMatrices(x)=Identity(1)
               endif
               if (CorePosition.lt.dims(2)) then
                  MPSAbove=onePEPS%GetColAsMPS(x,[CorePosition+1,dims(2)])
                  MPSBelow=belowPEPS%GetColAsMPS(x,[CorePosition+1,dims(2)])
                  aMultiplicator=new_Multiplicator(MPSAbove,MPSBelow)
                  RightMatrices(x)=RightAtSite(aMultiplicator,1)
                  if (Debug) Print *,'Right matrix x=',x,Norm(LeftMatrices(x)-Identity(onePEPS%GetBondAt(x,CorePosition,UP)))
               else
                  RightMatrices(x)=Identity(1)
               endif
            enddo
            CoreMPSBelow=GetRowAsMPS(belowPEPS,CorePosition)
            CoreMPSAbove=GetRowAsMPS(onePEPS,CorePosition,LeftMatrices,RightMatrices)
          end select

           theOverlap=Overlap(CoreMPSAbove,CoreMPSBelow)

        else
            call ThrowException('PEPS-Canonical_PEPS_Overlap','PEPS not initialized',NoErrorCode,CriticalError)
        endif


  end function Canonical_PEPS_Overlap


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!!!  CAREFUL FOR NOW THIS ROUTINE ASSUMES THAT PEPS IS CANONIZED IN COL/ROW OF CHOICE
  function Expectation_Value_Canonical_PEPS_MPO(onePEPS, aMPO, anotherPEPS,colORrowPosition,direction,corePosition) result(theResult)
      class(PEPS),intent(IN),target :: onePEPS
      class(MPO),intent(IN) :: aMPO
      class(PEPS),intent(IN),optional,target :: anotherPEPS
      integer,intent(IN) :: colORrowPosition,direction,corePosition
      type(PEPS),pointer :: belowPEPS
      type(MPS) :: MPSAbove,MPSBelow
      type(PEPSTensor) :: tempPEPSTensor
      type(PEPOTensor) :: tempPEPOTensor
      complex(8) :: theResult
      integer :: dims(2)

      dims=onePEPS%GetSize()

      if (present(anotherPEPS)) then
         belowPEPS => anotherPEPS
      else
         belowPEPS => onePEPS
      endif

      dims=onePEPS%GetSize()
      if (dims.equalvector.belowPEPS%GetSize()) then
        select case (direction)
         case (VERTICAL)
            MPSAbove=onePEPS%GetColAsMPS(colORrowPosition,[1,dims(2)])
            MPSBelow=belowPEPS%GetColAsMPS(colORrowPosition,[1,dims(2)])
            tempPEPOTensor=new_PEPOTensor(aMPO%GetTensorAt(CorePosition),3,4,VERTICAL)
            tempPEPSTensor=tempPEPOTensor.applyTo. (onePEPS%GetTensorAt(CorePosition, colORrowPosition))
            call MPSAbove%SetTensorAt(CorePosition,tempPEPSTensor%AsMPSTensor(DOWN,UP))
         case (HORIZONTAL)
            MPSAbove=onePEPS%GetRowAsMPS(colORrowPosition,[1,dims(1)])
            MPSBelow=belowPEPS%GetRowAsMPS(colORrowPosition,[1,dims(1)])
            tempPEPOTensor=new_PEPOTensor(aMPO%GetTensorAt(CorePosition),3,4,HORIZONTAL)
            tempPEPSTensor=tempPEPOTensor.applyTo.(onePEPS%GetTensorAt(CorePosition, colORrowPosition))
            call MPSAbove%SetTensorAt(CorePosition,tempPEPSTensor%AsMPSTensor(LEFT,RIGHT))
          end select

           theResult=Overlap(MPSAbove,MPSBelow)

        else
            call ThrowException('Expectation_Value_Canonical_PEPS_MPO','PEPS not initialized',NoErrorCode,CriticalError)
        endif

  end function Expectation_Value_Canonical_PEPS_MPO

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


end module PEPSAlgorithms_Class
