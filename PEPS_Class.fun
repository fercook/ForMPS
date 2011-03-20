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
test_suite PEPS_Class

setup

	  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown



test MPS_Core_Extraction

	use MPSAlgorithms_Class
	use Multiplicator_Class

    type(PEPS) :: aPEPS
    type(MPS) :: anMPS
    integer :: dims(2),x,y
    complex(8) :: theOverlap
    type(Multiplicator) :: aMultiplicator
    type(Tensor2) :: aMatrix,theId,LeftMatrices,RightMatrices
    type(PEPSTensor) :: aPTensor
    type(MPSTensor) :: aMtensor,bMTensor
    integer,parameter :: XLength=6,YLength=6,XcanonPos=4,YCanonPos=3

    aPEPS=new_PEPS(XLength,YLength,2,5)
    call aPEPS%CanonizeAt(XcanonPos,YCanonPos,HORIZONTAL,5,6)
    call aPEPS%PrintBondDimensions('PEPS map')
	anMPS=GetPEPSColAsMPS(aPEPS,XcanonPos,[1,YLength])
	call anMPS%PrintBondDimensions('Core MPS')

	aPTensor=aPEPS%GetTensorAt(XcanonPos,YcanonPos+1)
	aMatrix=aPTensor%CollapseAllIndicesBut(DOWN)
   dims=aMatrix%GetDimensions()
   TheId=Identity(dims(1))
   assert_equal_within(Norm(aMatrix-TheId),0.0d0,1.0d-12)
   aMtensor=aPTensor%AsMPSTensor(DOWN,UP)
   call aMTensor%PrintDimensions()
   print *,Norm(aMTensor),'NORM of MPS Tensor'
   bMtensor=anMPS%GetTensorAt(YcanonPos+1)
   print *,'Diff between direct and core mps tensor',norm(bMTensor-aMtensor)

   aMultiplicator=new_Multiplicator(anMPS)
   do y=1,YCanonPos-1
	   aMatrix=LeftAtSite(aMultiplicator,y)
	   dims=aMatrix%GetDimensions()
	   TheId=Identity(dims(1))
	   print *,'Diff of LEft at ',y,' with identity',Norm(aMatrix-TheId)
	   assert_equal_within(Norm(aMatrix-TheId),0.0d0,1.0d-12)
   enddo

   do y=YLength,YCanonPos+1,-1
   	aMatrix=RightAtSite(aMultiplicator,y)
	  dims=aMatrix%GetDimensions()
   	TheId=Identity(dims(1))
   	print *,'Diff of Right at ',y,' with identity',Norm(aMatrix-TheId)
   	assert_equal_within(Norm(aMatrix-TheId),0.0d0,1.0d-12)
   enddo

	theOverlap=Overlap_MPS(anMPS,anMPS)
	print *,theOverlap
	assert_equal_within(abs(theOverlap)**2,1.0d0,1.0d-6)

end test


test PEPS_Reduce_Bound_Dim
    type(PEPS) :: aPEPS,smallPEPS

    aPEPS=new_PEPS(4,4,2,4)
    smallPEPS=ReduceMAXPEPSBond(aPEPS,3)
    assert_equal(smallPEPS%GetMaxBond(),3)
    assert_true(smallPEPS%IsPEPSWellFormed())
    assert_true(aPEPS%IsPEPSWellFormed())
end test



test PEPS_Canonization_Routines_Horizontal
    type(PEPS) :: aPEPS,smallPEPS
    type(PEPSTensor) :: aTensor
    integer :: dims(4)
    type(Tensor2) :: aMatrix,theId

    aPEPS=new_PEPS(4,4,2,5)
    smallPEPS=aPEPS
    call smallPEPS%CanonizeAt(3,2,HORIZONTAL,6,6)
    dims=[4,2,36,16]
    aTensor=smallPEPS%GetTensorAt(3,2)
    assert_true(aTensor%GetBonds().equalvector.dims)
    assert_true(smallPEPS%IsPEPSWellFormed())
    call smallPEPS%PrintBondDimensions()

    !now check canonization of core
    aTensor=smallPEPS%GetTensorAt(3,4)
    aMatrix=aTensor%CollapseAllIndicesBut(DOWN)
    dims=aTensor%GetBonds()
    theId=Identity(dims(DOWN))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

    aTensor=smallPEPS%GetTensorAt(3,3)
    aMatrix=aTensor%CollapseAllIndicesBut(DOWN)
    dims=aTensor%GetBonds()
    theId=Identity(dims(DOWN))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

    aTensor=smallPEPS%GetTensorAt(3,1)
    aMatrix=aTensor%CollapseAllIndicesBut(UP)
    dims=aTensor%GetBonds()
    theId=Identity(dims(UP))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)
end test



test PEPS_Canonization_Routines_Vertical
    type(PEPS) :: aPEPS,smallPEPS
    type(PEPSTensor) :: aTensor
    integer :: dims(4)
    type(Tensor2) :: aMatrix,theId

    aPEPS=new_PEPS(4,4,2,5)
    smallPEPS=aPEPS
    call smallPEPS%CanonizeAt(3,2,VERTICAL,6,6)
    dims=[36,16,4,2]
    aTensor=smallPEPS%GetTensorAt(3,2)
    assert_true(aTensor%GetBonds().equalvector.dims)
    assert_true(smallPEPS%IsPEPSWellFormed())
    call smallPEPS%PrintBondDimensions()

    !now check canonization of core
    aTensor=smallPEPS%GetTensorAt(1,2)
    aMatrix=aTensor%CollapseAllIndicesBut(RIGHT)
    dims=aTensor%GetBonds()
    theId=Identity(dims(RIGHT))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

    aTensor=smallPEPS%GetTensorAt(2,2)
    aMatrix=aTensor%CollapseAllIndicesBut(RIGHT)
    dims=aTensor%GetBonds()
    theId=Identity(dims(RIGHT))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

    aTensor=smallPEPS%GetTensorAt(4,2)
    aMatrix=aTensor%CollapseAllIndicesBut(LEFT)
    dims=aTensor%GetBonds()
    theId=Identity(dims(LEFT))
    assert_equal_within(aMatrix.absdiff.TheId,0.0d0,1.0d-8)

end test


test MPS_Row_Extraction

	use Multiplicator_Class

    type(PEPS) :: aPEPS
    type(MPS) :: anMPS
    type(Multiplicator) :: aMultiplicator
    integer :: dims(2),theRow
    complex(8) :: theOverlap
    type(Tensor2) :: aMatrix,theId

    aPEPS=new_PEPS(6,4,2,5)
    call aPEPS%CanonizeAt(4,2,HORIZONTAL,6,6)
    do theRow=1,4
	    anMPS=GetRowAsMPS(aPEPS,theRow,[1,3])
    	aMultiplicator=new_Multiplicator(anMPS)
		aMatrix=LeftAtSite(aMultiplicator,3)
		dims=aMatrix%GetDimensions()
		TheId=Identity(dims(1))
		print *,'Diff with identity',Norm(aMatrix-TheId)
		assert_equal_within(Norm(aMatrix-TheId),0.0d0,1.0d-12)
	enddo
    do theRow=1,4
	    anMPS=GetRowAsMPS(aPEPS,theRow,[5,6])
    	aMultiplicator=new_Multiplicator(anMPS)
		aMatrix=RightAtSite(aMultiplicator,1)
		dims=aMatrix%GetDimensions()
		TheId=Identity(dims(1))
		print *,'Diff with identity',Norm(aMatrix-TheId)
		assert_equal_within(Norm(aMatrix-TheId),0.0d0,1.0d-12)
	enddo

end test


test PEPS_creation_deletion
   type(PEPS) :: aPEPS
   type(PEPSTensor) :: aTensor

    Print *,'Inside PEPS'
   aPEPS=new_PEPS(4,4,2,2)
   aTensor=aPEPS%GetTensorAt(3,4)
   call aTensor%PrintDimensions('Dims of a PEPSTensor')
   print *, aTensor%GetDimensions()
   assert_true(aTensor%GetDimensions().equalvector.[2,2,1,2,2])
   assert_false(WasThereError())
   print *,'After bussiness'
   assert_equal(aPEPS%delete(),Normal)
end test


end test_suite
