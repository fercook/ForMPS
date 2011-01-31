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
    integer,parameter :: BondDimension = 2, SpinDim = 2
    integer :: numberOfSteps,numberOfStepHalvings
    type(PEPS) :: thePEPS
    type(PEPO) :: currentPEPO
    real(8) :: tau,hField,energy,ElapsedTime

    integer :: n,m,step,halving

    print *,'Hello, initializing parameters...'
    !Set up initial values
    hField=10.0d0
    numberOfSteps=5
    numberOfStepHalvings=4
    tau=0.16d0
    elapsedTime=0.0d0

    print *,'Creating new PEPS object...'
    thePEPS=new_PEPS(Xsize,Ysize,SpinDim,BondDimension)

    print *,'Starting up engines...'
    do halving=1,numberOfStepHalvings
        print *,'Halving: ',halving

        do n=1,10
        !First step is different to accomplish 2nd order evol operator
        call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau/2.0d0,2.0d0*hField)
        call evolvePEPS(thePEPS,currentPEPO,BondDimension)

        !print *,'Preparing general PEPO'
        call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau,hField)
        do step=1,numberOfSteps-1 !One less because 2nd order Suzuki is used
            print *,'Tau:',tau,' -x- Step ',step,' of ',numberOfSteps
            call evolvePEPS(thePEPS,currentPEPO,BondDimension)
        enddo
        !Final step is special too
        call prepareIsingPEPO(currentPEPO,Xsize,Ysize,tau/2.0d0,0.0d0)
        call evolvePEPS(thePEPS,currentPEPO,BondDimension)
        ElapsedTime = ElapsedTime + tau*numberOfSteps

        energy = ComputeIsingEnergy(thePEPS,Xsize,Ysize,hField)
        print *,ElapsedTime,energy
        !Finally halve the time and double the steps
        enddo
        tau=tau/2.0d0
        numberOfSteps=2*numberOfSteps
    enddo

end program Ising_tester

