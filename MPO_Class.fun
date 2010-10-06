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

test_suite MPO_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test MPO_creation_deletion
   type(MPO) :: anMPO
   type(MPOTensor) :: aTensor

   anMPO=new_MPO(5,2,4)
   aTensor=anMPO%GetTensorAt(3)
   assert_true(aTensor%GetDimensions().equalvector.[4,4,2,2])
   assert_false(WasThereError())
   assert_equal(anMPO%delete(),Normal)
end test

test MPO_applied_to_MPS
    type(MPO) :: anMPO
    type(MPS) :: anMPS,anotherMPS
    type(MPSTensor) :: aTensor
    integer :: bondMPS=20,bondMPO=4, Length=6, spin=2

   anMPO=new_MPO(Length,spin,bondMPO)
   anMPS=new_MPS(Length,spin,bondMPS)
   anotherMPS=anMPO.ApplyMPOTo.anMPS
   aTensor=anotherMPS%GetTensorAt(Length/2)
   assert_true(aTensor%GetDimensions().equalvector.[bondMPO*bondMPS,bondMPO*bondMPS,spin])
   assert_false(WasThereError())
   assert_equal(anMPO%delete(),Normal)

end test
!
!test MPS_memory_usage
!   type(MPS) :: anMPS
!   type(MPSTensor) :: aTensor
!
!   print *,'Before allocation. Check memory now'
!!   pause
!   anMPS=new_MPS(200,2,100)
!   print *,'Memory allocated. Check memory now'
!!   pause
!   assert_false(WasThereError())
!   assert_equal(anMPS%delete(),Normal)
!   print *,'Memory released. Check memory now'
!!   pause
!
!end test
!

end test_suite

