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

test PEPS_checkpointing
   type(PEPS) :: aPEPS
   type(PEPSTensor) :: aTensor

   aPEPS=new_PEPS(4,4,2,2)
   call aPEPS%CheckPointState(ALLTENSORS)
   assert_false( ALL(aPEPS%HasTensorChangedAt))
   aTensor=new_PEPSTensor(2,4,4,4,4)
   call aPEPS%SetTensorAt(2,2,aTensor)
   assert_false( ALL(aPEPS%HasTensorChangedAt))
   assert_true( ANY(aPEPS%HasTensorChangedAt))
   assert_false(WasThereError())
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

end test_suite
