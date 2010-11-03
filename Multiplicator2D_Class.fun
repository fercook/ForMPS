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
  integer :: length=4,width=4,spin=2,bond=3
  type(Tensor4) :: theTemplate
  type(Tensor4) :: aTensor

  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS)
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
  integer :: length=4,width=4,spin=2,bond=3
  type(MPOTensor) :: aTensor,correctTensor
  type(PEPSTensor) :: aPEPSt

  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS)
  call theEnvironment%PrepareRowAsMPO(1)
  call theEnvironment%PrepareRowAsMPO(2)
  aTensor=theEnvironment%RowsAsMPO(1)%GetTensorAt(1)
  assert_true(aTensor%GetDimensions().equalvector.[1,bond**2,bond**2,1])
  aTensor=theEnvironment%RowsAsMPO(2)%GetTensorAt(3)
  assert_true(aTensor%GetDimensions().equalvector.[bond**2,bond**2,bond**2,bond**2])

end test


test UpperlowerMPSCreation
  type(PEPS) :: aPEPS
  type(Multiplicator2D) :: theEnvironment
  integer :: length=4,width=4,spin=2,bond=2
  type(MPO) :: anMPO
  type(MPS) :: anMPS,approxMPS
  type(MPSTensor) :: oneTensor,anotherTensor

  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS)

  anMPS=theEnvironment%MPS_Above(length+1)
  call theEnvironment%PrepareRowAsMPO(length)
  anMPO=theEnvironment%RowsAsMPO(length)
  anMPS=anMPO.ApplyMPOTo.anMPS
  call anMPS%Canonize()
  call anMPS%SetNorm(ONE)
  approxMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,length)
  oneTensor=anMPS%GetTensorAt(2)
  anotherTensor=approxMPS%GetTensorAt(2)
  assert_equal_within(oneTensor.absdiff.anotherTensor,0.0d0,1.0d-8)

  anMPS=theEnvironment%MPS_Above(length)
  call theEnvironment%PrepareRowAsMPO(length-1)
  anMPO=theEnvironment%RowsAsMPO(length-1)
  anMPS=anMPO.ApplyMPOTo.anMPS
  call anMPS%Canonize()
  call anMPS%SetNorm(ONE)
  approxMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,length-1)
  oneTensor=anMPS%GetTensorAt(2)
  anotherTensor=approxMPS%GetTensorAt(2)
  assert_equal_within(oneTensor.absdiff.anotherTensor,0.0d0,1.0d-8)

  anMPS=theEnvironment%MPS_Below(0)
  call theEnvironment%PrepareRowAsMPO(1)
  anMPO=theEnvironment%RowsAsMPO(1)
  anMPS=anMPS.ApplyMPOTo.anMPO
  call anMPS%Canonize()
  call anMPS%SetNorm(ONE)
  approxMPS=Multiplicator2D_RowAsLowerMPS(theEnvironment,1)
  oneTensor=anMPS%GetTensorAt(2)
  anotherTensor=approxMPS%GetTensorAt(2)
  call oneTensor%PrintDimensions('Exact tensor dims:')
  call anotherTensor%PrintDimensions('Approx tensor dims:')
  assert_equal_within(oneTensor.absdiff.anotherTensor,0.0d0,1.0d-8)

end test

test lateralEnvironments
  type(PEPS) :: aPEPS
  type(Multiplicator2D) :: theEnvironment
  integer :: length=4,width=4,spin=2,bond=2
  type(Tensor4) :: aTensor
  type(MPO) :: anMPO
  type(MPS) :: aboveMPS,belowMPS
  type(MPSTensor) :: aboveMPSTensor,belowMPSTensor
  type(MPOTensor) :: anMPOTensor

  aPEPS=new_PEPS(width,length,spin,bond)
  theEnvironment=new_Multiplicator2D(aPEPS)
!  aTensor = theEnvironment%LeftAt(1,1)

  call theEnvironment%PrepareRowAsMPO(1)
  anMPO=theEnvironment%RowsAsMPO(1)
  print *,'--------------- Up to 4'
  aboveMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,4)
  print *,'--------------- Up to 3'
  aboveMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,3)
  print *,'--------------- Up to 2'
  aboveMPS=Multiplicator2D_RowAsUpperMPS(theEnvironment,2)
  print *,'And here'
  belowMPS=Multiplicator2D_RowAsLowerMPS(theEnvironment,0)
  print *,'Now continuing'
  anMPOTensor = anMPO%GetTensorAt(1)
  aboveMPSTensor = aboveMPS%GetTensorAt(1)
  belowMPSTensor = belowMPS%GetTensorAt(1)
  call anMPOTensor%PrintDimensions('MPO dims')
  call aboveMPSTensor%PrintDimensions('MPS above dims')
  call belowMPSTensor%PrintDimensions('MPS below dims')

    print *,LEFT,' to ',RIGHT

  aTensor= theEnvironment%LeftAt(2,3)
  call aTensor%PrintDimensions()

end test

end test_suite
