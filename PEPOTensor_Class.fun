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

test_suite PEPOTensor_Class

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

  type(PEPOTensor) :: aPEPO
  type(PEPSTensor) :: aPEPS,anotherPEPS
  integer :: error, spin=2,LSbond=2,RSbond=3,USbond=3,DSbond=2
  integer :: LObond=3,RObond=3,UObond=3,DObond=2

  aPEPS=new_PEPSTensor(spin,LSbond,RSbond,USbond,DSbond)
  aPEPO=new_PEPOTensor(spin,LObond,RObond,UObond,DObond)
  print *,'Before Applying PEPO'
  !THIS IS THE KEY OPERATION TO TEST
  anotherPEPS=aPEPO.applyTo.aPEPS
  print *,'After Applying PEPO'
    !TODO: This is just a dimensionality test, the numbers could still be wrong
  call aPEPO%PrintDimensions("PEPO Dims")
  call aPEPS%PrintDimensions("PEPS Dims")
  call anotherPEPS%PrintDimensions("new PEPS Dims")
  assert_equal(anotherPEPS%getDLeft(),(aPEPS%getDLeft())*(aPEPO%getDLeft()))
  assert_equal(anotherPEPS%getDRight(),(aPEPS%getDRight())*(aPEPO%getDRight()))
  assert_equal(anotherPEPS%getDUp(),(aPEPS%getDUp())*(aPEPO%getDUp()))
  assert_equal(anotherPEPS%getDDown(),(aPEPS%getDDown())*(aPEPO%getDDown()))

  error=aPEPO%delete()
  assert_equal(error,Normal)
  error=aPEPS%delete()
  assert_equal(error,Normal)
  error=anotherPEPS%delete()
  assert_equal(error,Normal)
  assert_false(WasThereError())

end test

end test_suite
