!!   Copyright 2011 Fernando M. Cucchietti
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
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


program Ising_tester
    use PEPO_Class
    use PEPS_Class
    use Constants
    use ErrorHandling
    use IsingHelper

    implicit none

    integer,parameter :: Xsize=4,Ysize=4
    integer,parameter :: BondDimension = 2, SpinDim = 2, LongitudinalBond=6,TransverseBond=6
    integer :: numberOfSteps,numberOfStepHalvings,error
    integer :: problemCount
    type(PEPS) :: thePEPS
    type(PEPO) :: currentPEPO
    real(8) :: tau,hField,energy,ElapsedTime,preEnergy,bigEnergy
    real(8) :: theta
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localState(:,:,:,:,:)

    integer :: n,m,step,halving

    print *,'Hello, initializing parameters...'
    !Set up initial values
    hField=0.2d0
    numberOfSteps=400
    numberOfStepHalvings=4
    tau=0.05d0
    elapsedTime=0.0d0

!    !LETS DO SOME TESTING FIRST
!    thePEPS=new_PEPS(Xsize,Ysize,SpinDim,integerONE)
!    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
!    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
!    localState=ZERO
!    theta=0.31d0
!    localState(1,1,1,1,1)=Cos(theta/2.0d0)
!    localState(2,1,1,1,1)=Sin(theta/2.0d0)
!    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
!    do n=1,Xsize
!      do m=1,Xsize
!        call thePEPS%SetTensorAt(n,m,localPEPS)
!      enddo
!    enddo
!    do n=0,50
!        energy = ComputeIsingEnergy(thePEPS,Xsize,Ysize,10.0d0*n/50)
!        print *,10.0d0*n/50,energy
!    enddo
!    error= thePEPS%delete()

    print *,'Creating new PEPS object...'
    thePEPS=new_PEPS(Xsize,Ysize,SpinDim,BondDimension)
    call thePEPS%CanonizeAt(3,1,VERTICAL,LongitudinalBond,TransverseBond)
    call thePEPS%PrintBondDimensions('Initial PEPS diagram:')
    energy = ComputeIsingEnergy(thePEPS,Xsize,Ysize,hField)
    print *,'Starting energy (random):',energy
    print *,'Starting up engines...'
 !   do halving=1,numberOfStepHalvings
     !   print *,'Halving: ',halving

  !      do n=1,10
        !First step is different to accomplish 2nd order evol operator
!        call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau/2.0d0,2.0d0*hField)
!        call evolvePEPS(thePEPS,currentPEPO,BondDimension)

        !print *,'Preparing general PEPO'
        call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau,hField)
        write (*,'(A6,3x,A12,3x,A12,3x,A12)') 'Time','prevEnergy','newBigEnergy','TruncEnergy'
        problemCount=0
        do step=1,numberOfSteps !-1 !One less because 2nd order Suzuki is used
            !print *,step,' of',numberOfSteps
            preEnergy=energy
            call evolvePEPSCanonical(thePEPS,currentPEPO,LongitudinalBond,TransverseBond)
            ElapsedTime = ElapsedTime +tau
            energy = ComputeIsingEnergy(thePEPS,Xsize,Ysize,hField)
            if(  (energy-preEnergy).gt.1.0d-4) problemCount=problemCount+1
            if (problemCount.gt.20) then
                tau=tau/2.0d0
                call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau,hField)
            endif
            write (*,'(f6.3,3x,f12.7,3x,f12.7,3x,f12.7)') ElapsedTime,preEnergy,bigEnergy,energy
        enddo
        !Final step is special too
!        call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau/2.0d0,0.0d0)
!        call evolvePEPS(thePEPS,currentPEPO,BondDimension)

        !Finally halve the time and double the steps
   !     enddo
  !      tau=tau/2.0d0
  !      numberOfSteps=2*numberOfSteps
 !   enddo

end program Ising_tester

