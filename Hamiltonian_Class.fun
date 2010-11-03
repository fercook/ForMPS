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
!    along with FortranMPS.  If not, see <http://www.gnu.org/licenses/>.

test_suite Hamiltonian_Class




setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test HamiltonianCreation
  type(Hamiltonian) :: Ham
  integer error
  Ham=new_Hamiltonian()
  error=Ham%delete()
  assert_equal(error,Normal)
  assert_false(WasThereError())
end test

test AddingInteractions
  type(Hamiltonian) :: Ham
  type(Operator) :: 

end test




end test_suite
