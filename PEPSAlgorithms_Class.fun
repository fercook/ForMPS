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

test_suite PEPSAlgorithms_Class

setup
  !Set testing mode
  MaxErrorAllowed=Warning
end setup

teardown

end teardown

!test MPS_MemoryDrawn
!
!  type(Tensor3) :: atensor
!type(Tensor1) :: avector
!type(Tensor2) :: amatrix
!
!  print *,'Before requesting memory'
!  aVector=New_Tensor(16*256*256)
!  print *,'V'
!  aMatrix=New_Tensor(16*256,256)
!  print *,'M'
!  aTensor=New_Tensor(256,256,16)
!  print *,'After requesting memory'
!
!end test


test PEPS_Canonization_Routines
    type(PEPS) :: aPEPS,smallPEPS
    type(PEPSTensor) :: aTensor
    integer :: dims(4)
    complex(8) :: overlapBS

    aPEPS=new_PEPS(4,4,2,2)
    call Normalize(aPEPS)
    print *,'Norm of big overlap is',Overlap_PEPS(aPEPS,aPEPS)
    smallPEPS=aPEPS
    call smallPEPS%CanonizeAt(3,2,HORIZONTAL,6,6)
    overlapBS=Overlap_PEPS(smallPEPS,smallPEPS)
    print *,'Norm of canonical state is',overlapBS
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-9)
    call smallPEPS%PrintbondDimensions()
    assert_true(smallPEPS%IsPEPSWellFormed())
    overlapBS=Overlap_PEPS(aPEPS,smallPEPS)
    print *, 'overlap BTW peps and Canonized PEPS (EXACT): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-9)

    aPEPS=new_PEPS(4,4,2,3)
    call Normalize(aPEPS)
    smallPEPS=aPEPS
    print *,'About to canonize...'
    call smallPEPS%CanonizeAt(3,2,HORIZONTAL,6,6)
    overlapBS=Overlap_PEPS(aPEPS,smallPEPS)
    print *, 'overlap BTW peps and Canonized PEPS (H small bond): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-6)

end test


test OverlapAlgorithm
  type(PEPS) :: aPEPS
  integer :: length=4,width=4,spin=2,bond=2, error
  complex(8) :: overlap12

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(4.0d0)**(1.0d0/2.0d0))
  print *,'ABOUT TO COMPUTE OVERLAP'
  overlap12 = Overlap_PEPS(aPEPS)
  print *,'PREVIOUS OVERLAP',overlap12
  assert_false(abs(overlap12)**2.eq.1.0d0)

  print *,'About to NORMALIZE ----------'
  call Normalize(aPEPS)
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  print *,'OVERLAP AFTER NORMALIZATION',overlap12
!!  assert_equal_within(abs(overlap12)**2,1.0d0,1.0d-8)

  assert_false(WasThereError())

  error= aPEPS%Delete()

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(2.0d0)**(1.0d0/2.0d0))
  overlap12 = Overlap_PEPS(aPEPS)
  print *,'SECOND OVERLAP',overlap12

  print *,'About to NORMALIZE ----------'
  call Normalize(aPEPS)
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  print *,'OVERLAP AFTER NORMALIZATION',overlap12

end test

test Progressive_truncation
  type(PEPS) :: aPEPS,smallPEPS
  integer :: length=4,width=4,spin=2,bond=3, error
  complex(8) :: overlap12
  integer :: smallbond

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(4.0d0)**(1.0d0/2.0d0))
  call Normalize(aPEPS)
  do smallbond=bond-1,1,-1
      print *,'Reducing Bond ',bond,' to ',smallbond
      smallPEPS=ReduceMAXPEPSBond(aPEPS,smallbond)
      call Normalize(smallPEPS)
      overlap12 = Overlap_PEPS(aPEPS,smallPEPS)
      print *,'Bond: ',smallbond,', overlap: ',abs(overlap12)**2
      assert_true(abs(overlap12)**2.gt.0.9)
  enddo

  assert_false(WasThereError())

end test




test ExptValueProdHamilt
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=0.d0
    real(8),parameter :: theta=0.5d0
    integer,parameter :: Xsize=4,Ysize=4,row=2
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=1,OperatorDim=2,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim,BondDim,BondDim),AL(OperatorDim,integerONE,BondDim),AR(OperatorDim,BondDim,integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart(OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    !Prepare a product PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=pauliZ
    localTensor=new_PEPOTensor(localMatrix)
    do m=1,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)

    assert_equal_within(energy,cos(theta)**16,1.0d-3)

end test


test ExptValueLocalSumHamtn
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=1.0d0
    real(8),parameter :: theta=4.0d0
    integer,parameter :: Xsize=4,Ysize=4,row=2
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=2,OperatorDim=2,SpinDim=2
    integer :: n,m,k,i,j,p,dims(8)
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim,BondDim,BondDim),AL(OperatorDim,integerONE,BondDim),AR(OperatorDim,BondDim,integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart(OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ

    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,2,2)=ONE
        A(2,1,2)=ONE;
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;
        AR(1,2,1)=ONE;      AR(2,1,1)=ONE;

    !Prepare a product PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=identity
    localTensor=new_PEPOTensor(localMatrix)
    do m=2,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    deallocate(localMatrix)
    allocate(localMatrix(2,2,integerONE,BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do m=1,BondDim
     do i=1,SpinDim
     do j=1,SpinDim
      do k=1,OperatorDim
        localMatrix(i,j,1,m,1,1)=localMatrix(i,j,1,m,1,1)+AL(k,1,m)*spinPart(k,i,j)
      enddo
     enddo
     enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(1,1,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,integerONE,integerONE,integerONE))
    localMatrix=ZERO
    do m=1,BondDim
     do i=1,SpinDim
     do j=1,SpinDim
      do k=1,OperatorDim
        localMatrix(i,j,m,1,1,1)=localMatrix(i,j,m,1,1,1)+AR(k,m,1)*spinPart(k,i,j)
      enddo
     enddo
     enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(4,1,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do m=1,BondDim
    do p=1,BondDim
     do i=1,SpinDim
     do j=1,SpinDim
      do k=1,OperatorDim
        localMatrix(i,j,m,p,1,1)=localMatrix(i,j,m,p,1,1)+A(k,m,p)*spinPart(k,i,j)
      enddo
     enddo
     enddo
    enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(2,1,localTensor)
    call theH%SetTensorAt(3,1,localTensor)

    !call theH%PrintBondDimensions()

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)

    assert_equal_within(energy,4*cos(theta),1.0d-3)

end test




test ExptValueHamiltonian
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=0.17d0
    real(8),parameter :: theta=0.71d0
    integer,parameter :: Xsize=4,Ysize=4,row=1
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=3,OperatorDim=3,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim, BondDim, BondDim),AL( OperatorDim, integerONE, BondDim),AR( OperatorDim, BondDim, integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart( OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ; spinPart(3,:,:)=PauliX
    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,3,3)=ONE
        A(2,1,2)=ONE;       A(2,2,3)=ONE
        A(3,1,3)=field
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;  AL(3,1,3)=field;
        AR(1,3,1)=ONE;      AR(2,2,1)=ONE;  AR(3,1,1)=field;
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

    !Now set the second row to have the Ising like Hamiltonian
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
    call theH%SetTensorAt(2,row,localTensor)
    call theH%SetTensorAt(3,row,localTensor)

    !call theH%PrintBondDimensions()

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)

    assert_equal_within(energy,3*cos(theta)**2+4*field*sin(theta),1.0d-3)

end test


test ExptValueHamiltonianCol
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=0.17d0
    real(8),parameter :: theta=0.71d0
    integer,parameter :: Xsize=4,Ysize=4,col=1
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=3,OperatorDim=3,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim, BondDim, BondDim),AL( OperatorDim, integerONE, BondDim),AR( OperatorDim, BondDim, integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart( OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ; spinPart(3,:,:)=PauliX
    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,3,3)=ONE
        A(2,1,2)=ONE;       A(2,2,3)=ONE
        A(3,1,3)=field
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;  AL(3,1,3)=field;
        AR(1,3,1)=ONE;      AR(2,2,1)=ONE;  AR(3,1,1)=field;
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

    !Now set the first col to have the Ising like Hamiltonian
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
    call theH%SetTensorAt(col,2,localTensor)
    call theH%SetTensorAt(col,3,localTensor)

    !call theH%PrintBondDimensions()

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)
    assert_equal_within(energy,3*cos(theta)**2+4*field*sin(theta),1.0d-3)

end test


test Canonical_Overlap
    type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
    complex(8) :: overlapBS,overlapNormal

    onePEPS=new_PEPS(4,4,2,2)
    smallPEPS1=onePEPS
    call smallPEPS1%CanonizeAt(2,2,HORIZONTAL,6,6)
    overlapBS=Overlap(smallPEPS1,CorePosition=2, CoreDirection=VERTICAL)
    print *, 'overlap BTW Canonized PEPSs (EXACT): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    twoPEPS=new_PEPS(4,4,2,2)
    smallPEPS2=twoPEPS
    call smallPEPS2%CanonizeAt(2,2,HORIZONTAL,6,6)
    overlapNormal=Overlap(onePEPS,twoPEPS)
    overlapBS=Overlap(smallPEPS1,smallPEPS2,CorePosition=2, CoreDirection=VERTICAL)
    assert_equal_within(abs(overlapBS-overlapNormal),0.0d0,1.0d-10)

end test



end test_suite
