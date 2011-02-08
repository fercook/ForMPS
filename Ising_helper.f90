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
module IsingHelper

    use PEPSAlgorithms_Class
!    use Tensor_Class
    use PEPOTensor_Class
    use Constants
    use ErrorHandling

implicit none

interface fillOutMatrix
    module procedure fillOutMatrix2D
end interface

contains

subroutine prepareIsingPEPO(aPEPO,Xsize,Ysize,tau,field)
    type(PEPO),intent(INOUT) :: aPEPO
    real(8),intent(IN) :: tau,field
    integer,intent(IN) :: Xsize,Ysize
    type(PEPOTensor) :: localTensor
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8) :: identity(2,2),pauliZ(2,2),pauliX(2,2),A(2,2,2),AL(2,1,2),AR(2,2,1)
    complex(8) :: spinPart(2,2,2,2)
    integer :: n,m
    integer,parameter :: BondDim=2,SpinDim=2
    real(8) :: cosht,sinht,coshtH,sinhtH

    !All these constants should be prepared elsewhere...
    cosht=Cosh(tau);   sinht=Sinh(tau)
    coshtH=Cosh(tau*field/2.0d0);  sinhtH=Sinh(tau*field/2.0d0)

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE

    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE

    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    A=ZERO
        A(1,1,1)=cosht;    A(1,2,2)=sinht
        A(2,1,2)=sinht;    A(2,2,1)=cosht

    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=cosht;  AL(1,1,2)=ZERO
        AL(2,1,1)=ZERO;   AL(2,1,2)=sinht

        AR(1,1,1)=ONE;    AR(1,2,1)=ZERO
        AR(2,1,1)=ZERO;   AR(2,2,1)=ONE

    spinPart(1,1,:,:)=cosht*identity-sinhtH*pauliX
    spinPart(1,2,:,:)=cosht*pauliZ-sinhtH*matmul(pauliZ,pauliX)
    spinPart(2,1,:,:)=cosht*pauliZ-sinhtH*matmul(pauliZ,pauliX)
    spinPart(2,2,:,:)=cosht*identity-sinhtH*pauliX

    aPEPO=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    if (allocated(localMatrix)) deallocate(localMatrix)
    !First do boundary terms
    !left-bottom corner
    allocate(localMatrix(2,2,integerONE,BondDim,BondDim,integerONE))
    call fillOutMatrix(localMatrix,AL,AL,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    call aPEPO%SetTensorAt(1,1,localTensor)
    !right-bottom corner
    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,integerONE,BondDim,integerONE))
    call fillOutMatrix(localMatrix,AR,AL,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    call aPEPO%SetTensorAt(Xsize,1,localTensor)
    !left-top corner
    deallocate(localMatrix)
    allocate(localMatrix(2,2,integerONE,BondDim,integerONE,BondDim))
    call fillOutMatrix(localMatrix,AL,AR,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    call aPEPO%SetTensorAt(1,Ysize,localTensor)
    !right-top corner
    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,integerONE,integerONE,BondDim))
    call fillOutMatrix(localMatrix,AR,AR,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    call aPEPO%SetTensorAt(Xsize,Ysize,localTensor)


    !Bottom line Y=1
    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,BondDim,BondDim,integerONE))
    call fillOutMatrix(localMatrix,A,AL,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    do n=2,Xsize-1
        call aPEPO%SetTensorAt(n,1,localTensor)
    enddo
    !Top line Y=Ymax
    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,BondDim,integerONE,BondDim))
    call fillOutMatrix(localMatrix,A,AR,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    do n=2,Xsize-1
        call aPEPO%SetTensorAt(n,Ysize,localTensor)
    enddo
    !Left line X=1
    deallocate(localMatrix)
    allocate(localMatrix(2,2,integerONE,BondDim,BondDim,BondDim))
    call fillOutMatrix(localMatrix,AL,A,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    do n=2,Ysize-1
        call aPEPO%SetTensorAt(1,n,localTensor)
    enddo
    !Right line X=Xmax
    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,integerONE,BondDim,BondDim))
    call fillOutMatrix(localMatrix,AR,A,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    do n=2,Ysize-1
        call aPEPO%SetTensorAt(Xsize,n,localTensor)
    enddo

    !Bulk
    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,BondDim,BondDim,BondDim))
    call fillOutMatrix(localMatrix,A,A,spinPart)
    localTensor=new_PEPOTensor(localMatrix)
    do m=2,Ysize-1
      do n=2,Xsize-1
        call aPEPO%SetTensorAt(n,m,localTensor)
      enddo
    enddo

end subroutine prepareIsingPEPO

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine fillOutMatrix2D(matrix,Xtensor,Ytensor,operators)
    complex(8),intent(OUT) :: matrix(:,:,:,:,:,:)
    complex(8),intent(IN) :: Xtensor(:,:,:),Ytensor(:,:,:)
    complex(8),intent(IN) :: operators(:,:,:,:)
    integer :: MatrixDims(6),XTensorDims(3),YTensorDims(3),OpDims(4)
    integer :: i,j,u,d,l,r,n,m

    MatrixDims=shape(matrix)
    XTensorDims=shape(Xtensor)
    YTensorDims=shape(Ytensor)
    OpDims=shape(operators)

    matrix=ZERO
    do i=1,MatrixDims(1)
     do j=1,MatrixDims(2)
      do l=1,XTensorDims(2)
       do r=1,XTensorDims(3)
        do d=1,YTensorDims(2) !Left dim of YTensor is down
         do u=1,YTensorDims(3)!Right is up
          do n=1,XTensorDims(1)
           do m=1,YTensorDims(1)
             matrix(i,j,l,r,u,d)=matrix(i,j,l,r,u,d)+XTensor(n,l,r)*YTensor(m,d,u)*operators(n,m,i,j)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo

end subroutine fillOutMatrix2D


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine prepareIsingHamiltonianROW(theH,Xsize,Ysize,field,row)
    type(PEPO),intent(INOUT) :: theH
    real(8),intent(IN) :: field
    integer,intent(IN) :: Xsize,Ysize,row
    type(PEPOTensor) :: localTensor
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    integer,parameter :: BondDim=3,OperatorDim=3,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(2,2),pauliZ(2,2),pauliX(2,2)
    complex(8) :: A(OperatorDim,BondDim,BondDim),AL(OperatorDim,1,BondDim),AR(OperatorDim,BondDim,1)
    complex(8) :: AIdentity(1,1,1)
    complex(8) :: spinPart(OperatorDim,2,2)


    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ; spinPart(3,:,:)=PauliX
    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,3,3)=ONE
        A(2,1,2)=ONE;       A(2,2,3)=-ONE
        A(3,1,3)=field
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;  AL(3,1,3)=field;
        AR(1,3,1)=ONE;      AR(2,2,1)=-ONE;  AR(3,1,1)=field;
    !Prepare an identity PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=identity
    localTensor=new_PEPOTensor(localMatrix)
    do m=1,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    !Now set the row to have the Ising Hamiltonian
    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE, BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do n=1,1
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,n,m,1,1)=localMatrix(:,:,n,m,1,1)+AL(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(1,row,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,BondDim,integerONE,integerONE,integerONE))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,1
        do k=1,OperatorDim
          localMatrix(:,:,n,m,1,1)=localMatrix(:,:,n,m,1,1)+AR(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(Xsize,row,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,BondDim,BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,n,m,1,1)=localMatrix(:,:,n,m,1,1)+A(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    do n=2,Xsize-1
        call theH%SetTensorAt(n,row,localTensor)
    enddo

end subroutine prepareIsingHamiltonianROW


subroutine prepareIsingHamiltonianCOL(theH,Xsize,Ysize,field,col)
    type(PEPO),intent(INOUT) :: theH
    real(8),intent(IN) :: field
    integer,intent(IN) :: Xsize,Ysize,col
    type(PEPOTensor) :: localTensor
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    integer,parameter :: BondDim=3,OperatorDim=3,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(2,2),pauliZ(2,2),pauliX(2,2)
    complex(8) :: A(OperatorDim,BondDim,BondDim),AL(OperatorDim,1,BondDim),AR(OperatorDim,BondDim,1)
    complex(8) :: spinPart(OperatorDim,2,2)


    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ; spinPart(3,:,:)=PauliX
    A=ZERO
        A(1,1,1)=ONE;       A(1,3,3)=ONE
        A(2,1,2)=ONE;       A(2,2,3)=-ONE
        A(3,1,3)=field
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;  AL(3,1,3)=field;
        AR(1,3,1)=ONE;      AR(2,2,1)=-ONE;  AR(3,1,1)=field;
    !Prepare an identity PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=identity
    localTensor=new_PEPOTensor(localMatrix)
    do m=1,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    !Now set the col to have the Ising Hamiltonian
    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE,integerONE,BondDim,integerONE))
    localMatrix=ZERO
    do n=1,1
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,1,1,m,n)=localMatrix(:,:,1,1,m,n)+AL(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(col,1,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE,integerONE,integerONE,BondDim))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,1
        do k=1,OperatorDim
          localMatrix(:,:,1,1,m,n)=localMatrix(:,:,1,1,m,n)+AR(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(col,Ysize,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE,integerONE,BondDim,BondDim))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,1,1,m,n)=localMatrix(:,:,1,1,m,n)+A(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    do n=2,Ysize-1
        call theH%SetTensorAt(col,n,localTensor)
    enddo

end subroutine prepareIsingHamiltonianCOL

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine evolvePEPS(aPEPS,PEPOEvolver,BondDimension,preEnergy)
    type(PEPS),intent(INOUT) :: aPEPS
    type(PEPO),intent(IN) :: PEPOEvolver
    real(8),intent(OUT),optional :: preEnergy
    type(PEPS) :: tempPEPS
    integer,intent(IN) :: BondDimension
    integer :: error

    tempPEPS=PEPOEvolver.applyPEPOTo.aPEPS
    if (present(preEnergy)) preEnergy=real(ComputeIsingEnergy(tempPEPS,4,4,0.2d0)/overlap_PEPS(tempPEPS,tempPEPS))
    aPEPS=ReduceMAXPEPSBond(tempPEPS,BondDimension)
    error = tempPEPS%delete()
    call Normalize(aPEPS)

end subroutine evolvePEPS

subroutine evolvePEPSApprox(aPEPS,PEPOEvolver,BondDimension,preEnergy)
    type(PEPS),intent(INOUT) :: aPEPS
    type(PEPO),intent(IN) :: PEPOEvolver
    real(8),intent(OUT),optional :: preEnergy
    type(PEPS) :: tempPEPS
    integer,intent(IN) :: BondDimension
    integer :: error

    tempPEPS=PEPOEvolver.applyPEPOTo.aPEPS
    if (present(preEnergy)) preEnergy=real(ComputeIsingEnergy(tempPEPS,4,4,0.2d0)/overlap_PEPS(tempPEPS,tempPEPS))
    aPEPS=Approximate_PEPS(tempPEPS,BondDimension)
    error = tempPEPS%delete()
    call Normalize(aPEPS)

end subroutine evolvePEPSApprox

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

real(8) function ComputeIsingEnergy(aPEPS,Xsize,Ysize,field) result(energy)
    type(PEPS),intent(IN) :: aPEPS
    integer,intent(IN) :: Xsize,Ysize
    real(8) :: field
    integer :: row,col
    type(PEPO) :: lineHamiltonian

    energy=0.0d0
    do row=1,Ysize
        call prepareIsingHamiltonianROW(lineHamiltonian,Xsize,Ysize,field/2.0d0,row)
        energy=energy+ExpectationValue(aPEPS,lineHamiltonian)
    enddo
    do col=1,Xsize
        call prepareIsingHamiltonianCOL(lineHamiltonian,Xsize,Ysize,field/2.0d0,col)
        energy=energy+ExpectationValue(aPEPS,lineHamiltonian)
    enddo

end function ComputeIsingEnergy


end module IsingHelper
