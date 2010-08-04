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


test_suite Multiplicator_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion
  type(MPS),target :: anMPS
  type(Multiplicator) :: Prod
  integer :: length=10,spin=2,bond=20
  complex(8) :: overlap
  type(Tensor2) :: matrix

  anMPS=new_MPS(length,spin,bond)
  Prod=new_Multiplicator(anMPS)
  overlap=(Multiplicator_Left_Clean(Prod,1))**(Multiplicator_Right_Clean(Prod,2))
  print *,overlap

  call anMPS%Canonize()
  call Prod%Reset(LEFT)
  call Prod%Reset(RIGHT)
  matrix=(Multiplicator_Left_Clean(Prod,2))*(Multiplicator_Right_Clean(Prod,3))
  call matrix%PrintDimensions()
  overlap=(Multiplicator_Left_Clean(Prod,0))**(Multiplicator_Right_Clean(Prod,1))
  print *,overlap
  assert_equal_within(abs(overlap),1.0d0,1.0e-8)

  call Prod%delete()
end test

end test_suite
