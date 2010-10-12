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
test_suite MPS_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test MPS_creation_deletion
   type(MPS) :: anMPS
   type(MPSTensor) :: aTensor

   anMPS=new_MPS(5,2,7)
   aTensor=anMPS%GetTensorAt(3)
   assert_true(aTensor%GetDimensions().equalvector.[7,7,2])
   assert_false(WasThereError())
   assert_equal(anMPS%delete(),Normal)
end test


test MPS_canonicalForm
   type(MPS) :: anMPS
   type(MPSTensor) :: aTensor

   anMPS=new_MPS(5,2,7)
   call anMPS%Canonize()
   aTensor=anMPS%GetTensorAt(3)
   assert_true(aTensor%GetDimensions().equalvector.[4,4,2])
   assert_false(WasThereError())
   assert_equal(anMPS%delete(),Normal)
end test


test MPS_memory_usage
   type(MPS) :: anMPS
   type(MPSTensor) :: aTensor

   print *,'Before allocation. Check memory now'
!   pause
   anMPS=new_MPS(200,2,100)
   print *,'Memory allocated. Check memory now'
!   pause
   assert_false(WasThereError())
   assert_equal(anMPS%delete(),Normal)
   print *,'Memory released. Check memory now'
!   pause

end test


end test_suite
