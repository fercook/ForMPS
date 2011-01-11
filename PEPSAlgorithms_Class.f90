
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

    real(8) :: PSEUDOINVERSETOLERANCE = 1.0d-4

    interface Approximate
        module procedure Approximate_PEPS
    end interface

    interface Normalize
        module procedure Normalize_PEPS
    end interface

    interface Overlap
        module procedure Overlap_PEPS
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

    if(bigPEPS%IsInitialized()) then

        PEPSsize=bigPEPS%GetSize()
        PEPSSpin=bigPEPS%GetSpin()
        if (newBondDimension.ge.bigPEPS%GetMaxBond() ) then
            smallPEPS = new_PEPS(bigPEPS)
            return
        else
            smallPEPS= new_PEPS(PEPSSize(1),PEPSSize(2),PEPSSpin,newBondDimension)
        endif

        smallMultiplicator = new_Multiplicator2D(smallPEPS)
        bigMultiplicator = new_Multiplicator2D(bigPEPS,smallPEPS)

        sweep=0
        !Start sweeping
        !overlap=Abs(bigMultiplicator%FullContraction())**2
        do while (sweep .le. MaxSweeps ) !.and. (1.d0-overlap).gt.ApproximationTolerance )
            sweep=sweep+1
            !Sweep from 1,1 first to the right and then up
            do siteY=1,PEPSSize(2)
	            do siteX=1,PEPSSize(1)
	                print *,'Optimizing site ',siteX,',',siteY
	                localTensor = smallPEPS%GetTensorAt(siteX,siteY)
	                localPEPSDims = localTensor%GetDimensions()
	                call localTensor%Print('Previous tensor data')
	                aTempTensor=smallMultiplicator%LeftAt(siteX-1,siteY)
	                call aTempTensor%Print('Left environment data')
	                call aTempTensor%PrintDimensions('xxxxxxxx-------   0,1 dims LEFT small')
	                aTempTensor=smallMultiplicator%RightAt(siteX+1,siteY)
	                call aTempTensor%Print('Right environment data')
	                call aTempTensor%PrintDimensions('xxxxxxxx-------   2,1 dims RIGHT small')
                    aTempTensor= smallMultiplicator%AboveAt(siteX,siteY+1)
                    call aTempTensor%Print('Above environment data')
                    call aTempTensor%PrintDimensions('xxxxxxxx-------   1,2 dims ABOVE small')
                    aTempTensor=smallMultiplicator%BelowAt(siteX,siteY-1)
                    call aTempTensor%Print('Below environment data')
                    call aTempTensor%PrintDimensions('xxxxxxxx-------   1,0 dims BELOW small')
                    aTempTensor=bigMultiplicator%LeftAt(siteX-1,siteY)
                    call aTempTensor%Print('BIG Left environment data')
                    call aTempTensor%PrintDimensions('xxxxxxxx-------   0,1 dims LEFT BIG')
                    aTempTensor=bigMultiplicator%RightAt(siteX+1,siteY)
                    call aTempTensor%Print('BIG Right environment data')
                    call aTempTensor%PrintDimensions('xxxxxxxx-------   2,1 dims RIGHT BIG')
                    aTempTensor=bigMultiplicator%AboveAt(siteX,siteY+1)
                    call aTempTensor%Print('BIG Above environment data')
                    call aTempTensor%PrintDimensions('xxxxxxxx-------   1,2 dims ABOVE BIG')
                    aTempTensor=bigMultiplicator%BelowAt(siteX,siteY-1)
                    call aTempTensor%Print('BIG Below environment data')
                    call aTempTensor%PrintDimensions('xxxxxxxx-------   1,0 dims BELOW BIG')

                    localTensor = ComputePEPSOptimum ( smallMultiplicator%LeftAt(siteX-1,siteY), smallMultiplicator%RightAt(siteX+1,siteY), &
                        & smallMultiplicator%AboveAt(siteX,siteY+1), smallMultiplicator%BelowAt(siteX,siteY-1), &
                        & bigMultiplicator%LeftAt(siteX-1,siteY), bigMultiplicator%RightAt(siteX+1,siteY), &
                        & bigMultiplicator%AboveAt(siteX,siteY+1), bigMultiplicator%BelowAt(siteX,siteY-1), &
                        & bigPEPS%GetTensorAt(siteX,siteY), localPEPSDims )
                        call smallPEPS%SetTensorAt(siteX,siteY,localTensor)
                    call localTensor%Print('NEW Tensor data')
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
                        call smallPEPS%SetTensorAt(siteX,siteY,localTensor)
                enddo
                call smallMultiplicator%Reset(LEFT)
                call bigMultiplicator%Reset(LEFT)
            enddo
            !At the end of the sweep, the multiplicators are reset
            call smallMultiplicator%Reset(DOWN)
            call bigMultiplicator%Reset(DOWN)
            !overlap=Abs(bigMultiplicator%FullContraction())**2
        enddo
        print *,'---------------- Out of big loop ------------'
        call smallMultiplicator%Delete()
        call bigMultiplicator%Delete()
        if (present(returnOverlap)) then
            returnOverlap=overlap
        endif
        If (Verbose) print *,'PEPS-Approximate/Total sweeps performed: ',sweep
    else
        call ThrowException('Overlap algorithm','PEPS not initialized',NoErrorCode,CriticalError)
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

    call bigMatrixTimesVector%PrintDimensions('Big environment dims:')

    call smallMatrix%Print('Small matrix data:')
    call bigMatrixTimesVector%Print('Big matrix*vector data:')

    newTensor= new_PEPSTensor(  &
       & SplitIndexOf(   &
       & SolveLinearProblem(smallMatrix, bigMatrixTimesVector, PSEUDOINVERSETOLERANCE ), newDims )  &
       &                      )
        call newTensor%Print('New Tensor data:')
  end function ComputePEPSOptimum

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function Overlap_PEPS(onePEPS, anotherPEPS) result(theOverlap)
      class(PEPS),intent(INOUT) :: onePEPS
      class(PEPS),intent(INOUT),optional :: anotherPEPS
      type(Multiplicator2D) :: theEnvironment
      complex(8) :: theOverlap

      if (present(anotherPEPS)) then
          theEnvironment=new_Multiplicator2D(onePEPS, anotherPEPS)
      else
          theEnvironment=new_Multiplicator2D(onePEPS)
      endif
      theOverlap = Overlap_PEPSAboveBelow(theEnvironment)

  end function Overlap_PEPS

  subroutine Normalize_PEPS(aPEPS)
      class(PEPS),intent(INOUT) :: aPEPS
      type(Multiplicator2D) :: theEnvironment
      real(8) :: theNorm
      integer :: TotalNumberOfTensors

      theEnvironment=new_Multiplicator2D(aPEPS)
      theNorm = abs(Overlap_PEPSAboveBelow(theEnvironment))   !!!Notice I should not use **2

      TotalNumberOfTensors=product(aPEPS%GetSize())

      call aPEPS%ScaleBy(ONE/(theNorm**(0.5d0/TotalNumberOfTensors)))

  end subroutine Normalize_PEPS

end module PEPSAlgorithms_Class
