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

test_suite MPOTensor_Class

!TODO: New tests with all possible combinations of index bonding

! use ErrorHandling
! use Constants

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion

  type(MPOTensor) :: anMPO
  type(MPSTensor) :: anMPS,anotherMPS
  integer :: error, spin=2,LSbond=2,RSbond=7,LObond=3,RObond=5

  anMPS=new_MPSTensor(spin,LSbond,RSbond)
  anMPO=new_MPOTensor(spin,LObond,RObond)
  !THIS IS THE KEY OPERATION TO TEST
  anotherMPS=anMPO.applyTo.anMPS
  !TODO: This is just a dimensionality test, the numbers could still be wrong
  assert_equal(anotherMPS%getDLeft(),(anMPS%getDLeft())*(anMPO%getDLeft()))

  error=anMPO%delete()
  assert_equal(error,Normal)
  error=anMPS%delete()
  assert_equal(error,Normal)
  error=anotherMPS%delete()
  assert_equal(error,Normal)
  assert_false(WasThereError())

end test

end test_suite
