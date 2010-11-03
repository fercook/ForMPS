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

test_suite PEPSTensor_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion
  type(PEPSTensor) :: A,B
  integer :: spinT=2,BondL=3,BondR=4,bondU=2,bondD=3

  A=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
  assert_false(WasThereError())
  assert_equal(A%delete(),Normal)
  call LowerFlag()
  B=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
  call B%Print('TESTING PRINTING')
  assert_equal(B%delete(),Normal)
  assert_false(WasThereError())
end test

!

test accesors_work
  type(PEPSTensor) :: A
  integer :: spinT=2,BondL=3,BondR=4,bondU=2,bondD=3

  A=new_PEPSTensor(SpinT,BondL,BondR,bondU,bondD)
  assert_equal(A%getspin(),SpinT)
  assert_equal(A%getDLeft(),BondL)
  assert_equal(A%getDRight(),BondR)
  assert_equal(A%getDDown(),BondD)
  assert_equal(A%getDUp(),BondU)
  assert_equal(A%delete(),Normal)
  assert_false(WasThereError())
end test

test CreateATransposedTensor
    type(PEPSTensor) :: A
    type(Tensor5) :: aTensor
    integer :: spin=2, bondL=3, bondR=4,bondU=2,bondD=3
    integer :: dims(5)

    aTensor=new_Tensor(bondR,spin,bondL,bondD,bondU)
    A=new_PEPSTensor(aTensor,SecondDimension,ThirdDimension,FirstDimension,FifthDimension,FourthDimension)
    dims=A%GetDimensions()
    assert_true( dims .equalvector. [bondL, bondR, bondU, bondD, spin] )
  assert_false(WasThereError())
end test


end test_suite

