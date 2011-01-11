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

module MPSAlgorithms_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class
  use MPSTensor_Class
  use MPS_Class
  use Multiplicator_Class

  implicit none

    integer :: MaxSweeps = 10
    integer :: ApproximationTolerance = 1.0D-9

    interface Approximate
        module procedure Approximate_MPS
    end interface
  contains

!*****************************************************************************
  complex(8) function Overlap(mps1,mps2)
    type(MPS),intent(IN) :: mps1, mps2
    type(Multiplicator) :: aMultiplicator

    if(mps1%IsInitialized().and.mps2%IsInitialized()) then
        aMultiplicator = new_Multiplicator(mps1,mps2)
        Overlap=TensorTrace( (LeftAtSite(aMultiplicator,0))*(RightAtSite(aMultiplicator,1)) )
    else
        call ThrowException('Overlap algorithm','MPS not initialized',NoErrorCode,CriticalError)
    endif
    call aMultiplicator%Delete()

  end function Overlap


!*****************************************************************************
  function Approximate_MPS(bigMPS, newBondDimension,returnOverlap) result(smallMPS)
    !TODO: Notice that bigMPS must be canonized on site 1 before entry
    type(MPS), intent(IN) :: bigMPS
    integer,intent(IN) :: newBondDimension
    real(8),optional,intent(OUT) :: returnOverlap
    type(MPS) :: smallMPS
    type(MPSTensor) :: localTensor
    type(Multiplicator) :: smallMultiplicator,bigMultiplicator
    integer :: MPSLength,MPSSpin,sweeps
    real(8) :: overlap
    integer :: sweep,site
    type(Tensor2) :: TempTensor
    type(MPSTensor) :: TempMPS

    if(bigMPS%IsInitialized()) then

        MPSLength=bigMPS%GetSize()
        MPSSpin=bigMPS%GetSpin()
        if (newBondDimension.ge.bigMPS%GetMaxBond() ) then
            smallMPS = new_MPS(bigMPS)
            return
        else
            smallMPS= new_MPS(MPSLength,MPSSpin,newBondDimension)
            call smallMPS%Canonize()
        endif

        smallMultiplicator = new_Multiplicator(smallMPS)
        bigMultiplicator = new_Multiplicator(bigMPS,smallMPS)
        overlap=abs(TensorTrace( (LeftAtSite(bigMultiplicator,0))*(RightAtSite(bigMultiplicator,1)) ) )**2

        sweep=0
        !Start sweeping
        do while (sweep .le. MaxSweeps .and. (1.d0-overlap).gt.ApproximationTolerance )
            sweep=sweep+1
            !Sweep to the right
            do site=1,MPSLength
                localTensor=(LeftAtSite(bigMultiplicator,site-1)) .times. &
                    & (bigMPS%GetTensorAt(site)) .times. (RightAtSite(bigMultiplicator,site+1))
                call smallMPS%SetTensorAt(site,localTensor)
                call smallMPS%RightCanonizeAtSite(site)
            enddo
            !At the end of the sweep, the right multiplicators are reset
            call smallMultiplicator%Reset(RIGHT)
            call bigMultiplicator%Reset(RIGHT)
            !Sweep to the left
            do site=MPSLength,1,-1
                localTensor=(LeftAtSite(bigMultiplicator,site-1)) .times. &
                    & (bigMPS%GetTensorAt(site)) .times. (RightAtSite(bigMultiplicator,site+1))
                call smallMPS%SetTensorAt(site,localTensor)
                call smallMPS%LeftCanonizeAtSite(site)
            enddo
            call smallMultiplicator%Reset(LEFT)
            call bigMultiplicator%Reset(LEFT)
            overlap=abs(TensorTrace( (LeftAtSite(bigMultiplicator,0))*(RightAtSite(bigMultiplicator,1)) ) )**2
        enddo
        call smallMultiplicator%Delete()
        call bigMultiplicator%Delete()
        if (present(returnOverlap)) then
            returnOverlap=overlap
        endif
        If (Verbose) print *,'Total sweeps performed: ',sweep
    else
        call ThrowException('Overlap algorithm','MPS not initialized',NoErrorCode,CriticalError)
    endif

  end function Approximate_MPS

end module MPSAlgorithms_Class
