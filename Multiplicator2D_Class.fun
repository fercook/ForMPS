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


test_suite Multiplicator2D_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test up_down_boundaries_of_mult
  type(PEPS) :: aPEPS
  type(Multiplicator2D) :: theEnvironment
  integer,parameter :: length=4,width=4,spin=2,bond=3
  type(Tensor4) :: theTemplate
  type(Tensor4) :: aTensor
  logical :: matrixOfChanges(0:width+1,0:length+1)

  matrixOfChanges=.true.
  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS,MatrixToTrackChanges=matrixOfChanges)
  call theEnvironment%SetMaxApproxBond(20)
  theTemplate=new_Tensor(integerONE,integerONE,integerONE,integerONE,ONE)
  aTensor=theEnvironment%BelowAt(2,0)
  assert_equal_within(aTensor.absdiff.theTemplate,0.0d0,1.0d-8)
  aTensor=theEnvironment%AboveAt(3,5)
  assert_equal_within(aTensor.absdiff.theTemplate,0.0d0,1.0d-8)
  assert_false(WasThereError())

end test

test SetCurrentRow
  type(PEPS) :: aPEPS
  type(Multiplicator2D) :: theEnvironment
  integer,parameter :: length=4,width=4,spin=2,bond=3
  type(MPOTensor) :: aTensor,correctTensor
  type(PEPSTensor) :: aPEPSt
  integer :: dims(4),newdims(4)
  logical,target :: matrixOfChanges(0:width+1,0:length+1)

  matrixOfChanges = .true.
  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS,MatrixToTrackChanges=matrixOfChanges)
  call theEnvironment%SetMaxApproxBond(20)
  call theEnvironment%PrepareRowAsMPO(1)
  call theEnvironment%PrepareRowAsMPO(2)
  aTensor=theEnvironment%RowsAsMPO(1)%GetTensorAt(1)
  dims=[1,bond**2,bond**2,1]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)
  aTensor=theEnvironment%RowsAsMPO(2)%GetTensorAt(3)
  dims=[bond**2,bond**2,bond**2,bond**2]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)
  print *,'Set row test finished'
end test


test UpperlowerMPSCreation
  type(PEPS) :: aPEPS
  type(Multiplicator2D) :: theEnvironment
  integer,parameter :: length=4,width=4,spin=2,bond=2
  type(MPO) :: anMPO
  type(MPOTensor) :: anMPOTensor
  type(MPS) :: anMPS,approxMPS
  type(MPSTensor) :: oneTensor,anotherTensor
  logical,target :: matrixOfChanges(0:width+1,0:length+1)

  matrixOfChanges=.true.
  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS,MatrixToTrackChanges=matrixOfChanges)
  call theEnvironment%SetMaxApproxBond(20)

  anMPS=theEnvironment%MPS_Above(length+1)
  call theEnvironment%PrepareRowAsMPO(length)
  anMPO=new_MPO(length)
  anMPO=new_MPO(theEnvironment%RowsAsMPO(length))

  anMPS=anMPO.ApplyMPOTo.anMPS
  call anMPS%Canonize()
  approxMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,length)
  oneTensor=anMPS%GetTensorAt(2)
  anotherTensor=approxMPS%GetTensorAt(2)
  assert_equal_within(oneTensor.absdiff.anotherTensor,0.0d0,1.0d-8)

  anMPS=theEnvironment%MPS_Above(length)
  call theEnvironment%PrepareRowAsMPO(length-1)
  anMPO=theEnvironment%RowsAsMPO(length-1)
  anMPS=anMPO.ApplyMPOTo.anMPS
  call anMPS%Canonize()
  approxMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,length-1)
  oneTensor=anMPS%GetTensorAt(2)
  anotherTensor=approxMPS%GetTensorAt(2)
  assert_equal_within(oneTensor.absdiff.anotherTensor,0.0d0,1.0d-8)

  anMPS=theEnvironment%MPS_Below(0)
  call theEnvironment%PrepareRowAsMPO(1)
  anMPO=theEnvironment%RowsAsMPO(1)
  anMPS=anMPS.ApplyMPOTo.anMPO
  call anMPS%Canonize()
  approxMPS=Multiplicator2D_RowAsLowerMPS(theEnvironment,1)
  oneTensor=anMPS%GetTensorAt(2)
  anotherTensor=approxMPS%GetTensorAt(2)
  call oneTensor%PrintDimensions('Exact tensor dims:')
  call anotherTensor%PrintDimensions('Approx tensor dims:')
  assert_equal_within(oneTensor.absdiff.anotherTensor,0.0d0,1.0d-8)

end test

test lateralEnvironments     !DIMENSIONAL ONLY TESTING OF LEFT AND RIGHT
  type(PEPS) :: aPEPS
  type(Multiplicator2D) :: theEnvironment
  integer,parameter :: length=4,width=4,spin=2,bond=2
  type(Tensor4) :: aTensor
  type(MPO) :: anMPO
  type(MPS) :: aboveMPS,belowMPS
  type(MPSTensor) :: aboveMPSTensor,belowMPSTensor
  type(MPOTensor) :: anMPOTensor
  integer :: dims(4),newdims(4),dims3(3),newdims3(3)
  logical,target :: matrixOfChanges(0:width+1,0:length+1)

  matrixOfChanges=.true.
  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS,MatrixToTrackChanges=matrixOfChanges)
  call theEnvironment%SetMaxApproxBond(20)

  call theEnvironment%PrepareRowAsMPO(1)
  anMPO=theEnvironment%RowsAsMPO(1)
  aboveMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,4)
  aboveMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,3)
  aboveMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,2)
  belowMPS=Multiplicator2D_RowAsLowerMPS(theEnvironment,0)
  anMPOTensor = anMPO%GetTensorAt(1)
  aboveMPSTensor = aboveMPS%GetTensorAt(1)
  belowMPSTensor = belowMPS%GetTensorAt(1)
  dims=[1,bond**2,bond**2,1]
  newdims=anMPOTensor%GetDimensions()
  assert_true(dims.equalvector.newdims)
  dims3=[1,bond**2,bond**2]
  newdims3=aboveMPSTensor%GetDimensions()
  assert_true(dims3.equalvector.newdims3)
  dims3=[1,1,1]
  newdims3=belowMPSTensor%GetDimensions()
  assert_true(dims3.equalvector.newdims3)

  aTensor= theEnvironment%LeftAt(2,3)
  dims=[bond**4,bond**2,bond,bond]
  newdims=aTensor%GetDimensions()
  assert_true(dims.equalvector.newdims)
  aTensor= theEnvironment%RightAt(3,3)
  dims=[bond**2,bond**4,bond,bond]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)

  aTensor= theEnvironment%LeftAt(2,2)
  dims=[bond**2,bond**4,bond,bond]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)
  aTensor= theEnvironment%RightAt(3,2)
  dims=[bond**4,bond**2,bond,bond]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)

  aTensor= theEnvironment%LeftAt(2,1)
  dims=[1,bond**4,bond,bond]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)
  aTensor= theEnvironment%RightAt(3,1)
  dims=[bond**4,1,bond,bond]
  newdims=aTensor%GetDimensions()
  assert_true(newdims.equalvector.dims)

end test

test Multip_With2_PEPS
  type(PEPS) :: aPEPS, aBigPEPS
  integer,parameter :: height=4,width=4,spin=2,bondSmall=2,bondBig=3
  integer:: error
  real(8) :: overlap12
  type(Multiplicator2D) :: theEnvironment
  type(Tensor4) :: aTensor
  integer :: dims(4),newdims(4)
  logical :: matrixOfChanges(0:width+1,0:height+1)

  aPEPS=new_PEPS(width,height,spin,bondSmall)
  aBigPEPS=new_PEPS(width,height,spin,bondBig)
  matrixOfChanges=.true.
  theEnvironment = new_Multiplicator2D(aBigPEPS,aPEPS,MatrixToTrackChanges=matrixOfChanges)
  call theEnvironment%SetMaxApproxBond(20)

  aTensor= theEnvironment%LeftAt(1,1)
  call aTensor%PrintDimensions('1-1 Left Dimensions')
  aTensor= theEnvironment%LeftAt(0,1)
  call aTensor%PrintDimensions('0-1 Left Dimensions')
  aTensor= theEnvironment%BelowAt(1,0)
  call aTensor%PrintDimensions('1-0 Below Dimensions')
  aTensor= theEnvironment%RightAt(2,1)
  call aTensor%PrintDimensions('2-1 Right Dimensions')
!
  aTensor= theEnvironment%LeftAt(2,2)
  call aTensor%PrintDimensions('Left Dimensions')
  dims=[bondSMall*BondBig,min((bondSMall*BondBig)**2,theEnvironment%GetMaxApproxBond()),bondBig,bondSmall]
  newdims=aTensor%GetDimensions()
  assert_true(dims.equalvector.newdims)
  aTensor= theEnvironment%RightAt(3,2)
  call aTensor%PrintDimensions('Right Dimensions')
  dims=[min((bondSMall*BondBig)**2,theEnvironment%GetMaxApproxBond()),bondSMall*BondBig,bondBig,bondSmall]
  newdims=aTensor%GetDimensions()
  assert_true(dims.equalvector.newdims)

  assert_false(WasThereError())

end test

end test_suite
