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
  type(MPS) :: anMPS
  type(Multiplicator) :: Prod
  integer :: length=10,spin=2,bond=20
  complex(8) :: overlap
  type(Tensor2) :: matrix

  anMPS=new_MPS(length,spin,bond)
  Prod=new_Multiplicator(anMPS)
  overlap=(LeftAtSite(Prod,1)).xx.(Multiplicator_Right_Clean(Prod,2))
  print *,overlap
  matrix=(Multiplicator_Left_Clean(Prod,2)).x.(Multiplicator_Right_Clean(Prod,3))
  call matrix%PrintDimensions()

  call anMPS%Canonize()
  call Prod%Reset(LEFT)
  call Prod%Reset(RIGHT)
  matrix=(Multiplicator_Left_Clean(Prod,2)).x.(Multiplicator_Right_Clean(Prod,3))
  call matrix%PrintDimensions()
  overlap=(Multiplicator_Left_Clean(Prod,0)).xx.(Multiplicator_Right_Clean(Prod,1))
  print *,overlap
  assert_equal_within(abs(overlap),1.0d0,1.0e-8)

  call Prod%delete()
end test

test twoMPSsizes
  type(MPS) :: anMPS,anotherMPS
  type(Multiplicator) :: Prod
  integer :: length=10,spin=2,bondU=20,bondD=10
  complex(8) :: overlap
  type(Tensor2) :: matrix
  integer :: dimensions(2)

  anMPS=new_MPS(length,spin,bondU)
  anotherMPS=new_MPS(length,spin,bondD)

  call anMPS%Canonize()
  call anotherMPS%Canonize()

  Prod=new_Multiplicator(anMPS,anotherMPS)
  matrix=(LeftAtSite(Prod,5))
  dimensions=matrix%GetDimensions()
  assert_equal(dimensions(1), 10 )
  assert_equal(dimensions(2), 20 )
  matrix=(RightAtSite(Prod,6))
  dimensions=matrix%GetDimensions()
  assert_equal(dimensions(2), 10 )
  assert_equal(dimensions(1), 20 )

  call Prod%delete()
end test

end test_suite
