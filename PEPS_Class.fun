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
test_suite PEPS_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test PEPS_creation_deletion
   type(PEPS) :: aPEPS
   type(PEPSTensor) :: aTensor

    Print *,'Inside PEPS'
   aPEPS=new_PEPS(4,4,2,2)
   aTensor=aPEPS%GetTensorAt(3,4)
   call aTensor%PrintDimensions('Dims of a PEPSTensor')
   print *, aTensor%GetDimensions()
   assert_true(aTensor%GetDimensions().equalvector.[2,2,1,2,2])
   assert_false(WasThereError())
   print *,'After bussiness'
   assert_equal(aPEPS%delete(),Normal)
end test


test PEPS_Reduce_Bound_Dim
    type(PEPS) :: aPEPS,smallPEPS

    aPEPS=new_PEPS(4,4,2,4)
    smallPEPS=ReduceMAXPEPSBond(aPEPS,3)
    assert_equal(smallPEPS%GetMaxBond(),3)
    assert_true(smallPEPS%IsPEPSWellFormed())
    assert_true(aPEPS%IsPEPSWellFormed())
end test



test PEPS_Canonization_Routines_Basic
    type(PEPS) :: aPEPS,smallPEPS
    type(PEPSTensor) :: aTensor
    integer :: dims(4)
    type(Tensor2) :: aMatrix,theId

    aPEPS=new_PEPS(4,4,2,5)
    smallPEPS=aPEPS
    call smallPEPS%CanonizeAt(3,2,HORIZONTAL,6,6)
    dims=[4,2,36,16]
    aTensor=smallPEPS%GetTensorAt(3,2)
    assert_true(aTensor%GetBonds().equalvector.dims)
    assert_true(smallPEPS%IsPEPSWellFormed())
    call smallPEPS%PrintBondDimensions()

    !now check canonization of core
    aTensor=smallPEPS%GetTensorAt(3,4)
    aMatrix=aTensor%CollapseAllIndicesBut(DOWN)
    dims=aTensor%GetBonds()
    theId=Identity(dims(DOWN))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

    aTensor=smallPEPS%GetTensorAt(3,3)
    aMatrix=aTensor%CollapseAllIndicesBut(DOWN)
    dims=aTensor%GetBonds()
    theId=Identity(dims(DOWN))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

    aTensor=smallPEPS%GetTensorAt(3,1)
    aMatrix=aTensor%CollapseAllIndicesBut(UP)
    dims=aTensor%GetBonds()
    theId=Identity(dims(UP))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)
end test


end test_suite
