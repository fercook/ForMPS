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

test_suite PEPO_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test PEPO_creation_deletion
   type(PEPO) :: aPEPO
   type(PEPOTensor) :: aTensor
   integer :: dims(6),RIGHTdims(6)

    aPEPO=new_PEPO(4,4,2,4)
    aTensor=GetPEPOTensorAtSite(aPEPO,2,4)
    dims=aTensor%GetDimensions()
    RIGHTdims=[4,4,1,4,2,2]
    assert_true(dims.equalvector.RIGHTdims)
    assert_false(WasThereError())
    Print *,'PEPO_creation_deletion'
    assert_equal(aPEPO%delete(),Normal)
end test

test PEPO_applied_to_PEPS
    type(PEPO) :: aPEPO
    type(PEPS) :: aPEPS,anotherPEPS
    type(PEPSTensor) :: aTensor
    integer :: bondPEPS=4,bondPEPO=3, Length=4, spin=2
    integer :: dims(5),RIGHTdims(5)

   aPEPO=new_PEPO(Length,Length,spin,bondPEPO)
   aPEPS=new_PEPS(Length,Length,spin,bondPEPS)
   anotherPEPS=aPEPO.ApplyPEPOTo.aPEPS
   aTensor=anotherPEPS%GetTensorAt(Length/2,Length/2)
    dims=aTensor%GetDimensions()
    RIGHTdims=[bondPEPO*bondPEPS,bondPEPO*bondPEPS,bondPEPO*bondPEPS,bondPEPO*bondPEPS,spin]
    assert_true(dims.equalvector.RIGHTdims)
   assert_false(WasThereError())
!   assert_equal(aPEPO%delete(),Normal)
!   assert_equal(aPEPS%delete(),Normal)
!   assert_equal(anotherPEPS%delete(),Normal)

end test


end test_suite

