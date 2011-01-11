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

test_suite PEPSAlgorithms_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test OverlapAlgorithm
  type(PEPS) :: aPEPS
  integer :: length=4,width=4,spin=2,bond=2, error
  complex(8) :: overlap12

  aPEPS=new_PEPS(length,width,spin,bond)
  overlap12 = Overlap_PEPS(aPEPS)
  print *,overlap12
  assert_false(abs(overlap12)**2.eq.1.0d0)

  call Normalize(aPEPS)
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  assert_equal_within(abs(overlap12)**2,1.0d0,1.0d-8)

  assert_false(WasThereError())

  error= aPEPS%Delete()

end test

test Progressive_truncation
  type(PEPS) :: aPEPS,smallPEPS
  integer :: length=4,width=4,spin=2,bond=6, error
  complex(8) :: overlap12
  integer :: smallbond

  aPEPS=new_PEPS(length,width,spin,bond)
  call Normalize(aPEPS)

  do smallbond=6,1,-1
      smallPEPS=ReduceMAXPEPSBond(aPEPS,smallbond)
      call Normalize(smallPEPS)
      overlap12 = Overlap_PEPS(aPEPS,smallPEPS)
      print *,'Bond: ',smallbond,', overlap: ',abs(overlap12)**2
  enddo

  assert_false(WasThereError())

end test

test  ApproximationAlgorithm
  type(PEPS) :: aPEPS, aBigPEPS
  integer :: length=4,width=4,spin=2,bondSmall=2,bondBig=4, error
  real(8) :: overlap12

  aBigPEPS=new_PEPS(width,length,spin,bondBig)
  call Normalize(aBigPEPS)

  print *,'About to approximate ---------'
  aPEPS=Approximate(aBigPEPS,bondSmall,overlap12)
  print *, 'Approximation overlaps to :', overlap12
  assert_false(abs(overlap12)**2.eq.1.0d0)

  error= aPEPS%Delete()
  error= aBigPEPS%Delete()
  assert_false(WasThereError())

end test

!test ApproximationAlgorithm
!  type(MPS) :: smallMPS,bigMPS
!  integer :: length=20,spin=2,bondBig=40,bondSmall=2,site
!  real(8) :: overlap12
!  type(MPSTensor) :: localTensor
!
!  bigMPS=new_MPS(length,spin,bondBig)
!  call bigMPS%Canonize()
!  call bigMPS%SetNorm(ONE)
!
!  smallMPS=Approximate(bigMPS,bondSmall,overlap12)
!  assert_equal_within(overlap12,1.0d0,1.0D-5)
!
!  print *,'Approximated overlap :',overlap12
!  print *,'Big bond:',bondBig
!  print *,'Small bond: ',smallMPS%GetMaxBond()
!
!  do site=1,smallMPS%GetSize()
!    localTensor=smallMPS.TensorAt.site
!    call localTensor%PrintDimensions('Dimensions of tensor')
!    print *,site
!  enddo
!end test

end test_suite
