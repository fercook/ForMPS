!!   Copyright 2010 Fernando M. Cucchietti
!
!    This file is part of FortranMPS
!
!    FortranMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    FortranMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    test_suite MatrixProductState_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion

  type(MatrixProductState) :: mps
  integer error
  mps=new_MatrixProductState(10,2,17)
  error=mps%delete()
  assert_equal(error,Normal)
  Print *,' And here'
  assert_false(WasThereError())

end test

test assignments_of_mps

  type(MatrixProductState) :: mps1,mps2
  integer error
  mps1=new_MatrixProductState(10,2,10)
  mps2=mps1
  assert_equal(mps1%delete(),Normal)
  assert_equal(mps2%delete(),Normal)
  assert_false(WasThereError())

end test


test canonization_of_mps
  type(MatrixProductState) :: mps
  integer site
  real(8) result
  mps=new_MatrixProductState(10,2,20)
  assert_false(mps%isCanonized())
  result=mps%RCanonize()
  result=mps%LCanonize()
  assert_equal_within(result, 1.0d0, 1.0e-10)
  assert_true(mps%isCanonized())
  result=mps%CanonizeAtSite(5)
  assert_true(mps%isCanonized())
  assert_equal_within(result, 5.0d0, 1.0e-10)
end test



end test_suite
