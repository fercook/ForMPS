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

test_suite PEPSTensor_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion
  type(PEPSTensor) :: A,B
  integer :: spinT=2,BondL=3,BondR=4,bondU=2,bondD=3

  A=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
  assert_false(WasThereError())
  assert_equal(A%delete(),Normal)
  call LowerFlag()
  B=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
!  call B%Print('TESTING PRINTING')
  assert_equal(B%delete(),Normal)
  assert_false(WasThereError())
end test

!

test accesors_work
  type(PEPSTensor) :: A
  integer :: spinT=2,BondL=3,BondR=4,bondU=2,bondD=3

  A=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
  assert_equal(A%getspin(),SpinT)
  assert_equal(A%getDLeft(),BondL)
  assert_equal(A%getDRight(),BondR)
  assert_equal(A%getDDown(),BondD)
  assert_equal(A%getDUp(),BondU)
  assert_equal(A%delete(),Normal)
  assert_false(WasThereError())
end test

test CreateATransposedTensor
    type(PEPSTensor) :: A
    type(Tensor5) :: aTensor
    integer :: spinT=2, bondL=3, bondR=4,bondU=2,bondD=3
    integer :: dims(5)

    aTensor=new_Tensor(bondR,spinT,bondL,bondD,bondU)
    A=new_PEPSTensor(aTensor,SecondDimension,ThirdDimension,FirstDimension,FifthDimension,FourthDimension)
    dims=A%GetDimensions()
    assert_true( dims .equalvector. [bondL, bondR, bondU, bondD, spinT] )
  assert_false(WasThereError())
end test

test HOSVD_of_PEPS
    type(PEPSTensor) :: aTensor
    type(PEPSTensor) :: theCore, reconstructedPEPS
    type(Tensor2) :: Umatrices(4)
    integer :: spinT=2, bondL=3, bondR=3,bondU=3,bondD=3

    aTensor=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
    call aTensor%HOSVD(theCore,Umatrices)
    reconstructedPEPS=nModeProduct(Umatrices(4), &
                        & nModeProduct(Umatrices(3),  &
                          &  nModeProduct(Umatrices(2), &
                            &  nModeProduct(Umatrices(1),theCore, &
                            &  FIRST), &
                          &  SECOND), &
                         & THIRD), &
                      & FOURTH)
    assert_equal_within(reconstructedPEPS.absdiff.aTensor,0.0d0,1.0d-8 )
    assert_false(WasThereError())

end test

test Multiplication_By_MPO
   type(PEPSTensor) :: aPEPSt,newPEPSt
   type(Tensor4) :: aTensor
   integer :: spinT=2, bondL=2, bondR=2,bondU=2,bondD=2, dims(4),correctDims(4)

    aPEPSt=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)

    aTensor=new_Tensor(3,2,4,5)
    newPEPSt=ApplyMPOToBond(aPEPSt,aTensor,LEFT)
    dims=newPEPSt%getBonds()
    correctDims=[3,2,8,10]
    assert_true(dims.equalvector.correctDims)

    aTensor=new_Tensor(2,3,4,5)
    newPEPSt=ApplyMPOToBond(aPEPSt,aTensor,RIGHT)
    dims=newPEPSt%getBonds()
    correctDims=[2,3,8,10]
    assert_true(dims.equalvector.correctDims)

    aTensor=new_Tensor(3,4,5,2)
    newPEPSt=ApplyMPOToBond(aPEPSt,aTensor,UP)
    dims=newPEPSt%getBonds()
    correctDims=[6,8,5,2]
    assert_true(dims.equalvector.correctDims)

    aTensor=new_Tensor(3,4,2,5)
    newPEPSt=ApplyMPOToBond(aPEPSt,aTensor,DOWN)
    dims=newPEPSt%getBonds()
    correctDims=[6,8,2,5]
    assert_true(dims.equalvector.correctDims)
end test



test go_to_mps_and_back
  type(PEPSTensor) :: aPEPS,anotherPEPS
  type(MPSTensor) :: aMPSt
  integer :: error, spin=2,LSbond=3,RSbond=4,USbond=5,DSbond=6

  aPEPS=new_PEPSTensor(spin,LSbond,RSbond,USbond,DSbond)
  aMPSt=aPEPS%AsMPSTensor(LEFT,RIGHT)
  call aMPSt%PrintDimensions()
  anotherPEPS=new_PEPSTensor(aMPSt,3,HORIZONTAL,[USbond,DSbond])
  assert_equal_within(Norm(aPEPS-anotherPEPS),0.0d0,1.0d-15)

  aMPSt=aPEPS%AsMPSTensor(DOWN,UP)
  call aMPSt%PrintDimensions()
  anotherPEPS=new_PEPSTensor(aMPSt,3,VERTICAL,[LSbond,RSbond])
  call anotherPeps%PrintDimensions()
  assert_equal_within(Norm(aPEPS-anotherPEPS),0.0d0,1.0d-15)

end test


end test_suite

