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

!test MPS_MemoryDrawn
!
!  type(Tensor3) :: atensor
!type(Tensor1) :: avector
!type(Tensor2) :: amatrix
!
!  print *,'Before requesting memory'
!  aVector=New_Tensor(16*256*256)
!  print *,'V'
!  aMatrix=New_Tensor(16*256,256)
!  print *,'M'
!  aTensor=New_Tensor(256,256,16)
!  print *,'After requesting memory'
!
!end test

test OverlapAlgorithm
  type(PEPS) :: aPEPS
  integer :: length=4,width=4,spin=2,bond=2, error
  complex(8) :: overlap12

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(4.0d0)**(1.0d0/2.0d0))
  print *,'ABOUT TO COMPUTE OVERLAP'
  overlap12 = Overlap_PEPS(aPEPS)
  print *,'PREVIOUS OVERLAP',overlap12
  assert_false(abs(overlap12)**2.eq.1.0d0)

  print *,'About to NORMALIZE ----------'
  call Normalize(aPEPS)
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  print *,'OVERLAP AFTER NORMALIZATION',overlap12
!!  assert_equal_within(abs(overlap12)**2,1.0d0,1.0d-8)

  assert_false(WasThereError())

  error= aPEPS%Delete()

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(2.0d0)**(1.0d0/2.0d0))
  overlap12 = Overlap_PEPS(aPEPS)
  print *,'SECOND OVERLAP',overlap12

  print *,'About to NORMALIZE ----------'
  call Normalize(aPEPS)
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  print *,'OVERLAP AFTER NORMALIZATION',overlap12

end test

test Progressive_truncation
  type(PEPS) :: aPEPS,smallPEPS
  integer :: length=4,width=4,spin=2,bond=3, error
  complex(8) :: overlap12
  integer :: smallbond

  print *,' XXXXXXX HELLO TRUNCATION X----------'
  aPEPS=new_PEPS(length,width,spin,bond)
  print *,' XXXXXXX HELLO TRUNCATION 1----------'
  call aPEPS%ScaleBy(ONE/(4.0d0)**(1.0d0/2.0d0))
  print *,' XXXXXXX HELLO TRUNCATION 2----------'
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  call Normalize(aPEPS)
  print *,' XXXXXXX HELLO TRUNCATION 3----------'
  do smallbond=bond,1,-1
      smallPEPS=ReduceMAXPEPSBond(aPEPS,smallbond)
      call Normalize(smallPEPS)
      overlap12 = Overlap_PEPS(aPEPS,smallPEPS)
      print *,'Bond: ',smallbond,', overlap: ',abs(overlap12)**2
  enddo

  assert_false(WasThereError())

end test
!
!test  ApproximationAlgorithm
!  type(PEPS) :: aPEPS, aBigPEPS
!  integer :: length=4,width=4,spin=2,bondSmall=2,bondBig=4, error
!  real(8) :: overlap12
!
!  aBigPEPS=new_PEPS(width,length,spin,bondBig)
!  call Normalize(aBigPEPS)
!
!  print *,'About to approximate ---------'
!  aPEPS=Approximate(aBigPEPS,bondSmall,overlap12)
!  print *, 'Approximation overlaps to :', overlap12
!  assert_false(abs(overlap12)**2.eq.1.0d0)
!
!  error= aPEPS%Delete()
!  error= aBigPEPS%Delete()
!  assert_false(WasThereError())
!
!end test

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
