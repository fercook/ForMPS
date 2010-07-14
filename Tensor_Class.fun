test_suite Tensor_Class

setup
  !use ErrorHandling
  !use Constants
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion

  type(tensor3) :: aTensor
  integer error
  aTensor=new_Tensor(2,20,20)
  error=aTensor%delete()
  assert_equal(error,Normal)
  assert_false(WasThereError())

end test

test assignments_of_tensor

  type(tensor3) :: mps1,mps2
  integer error
  mps1=new_Tensor(10,2,10)
  mps2=mps1
  assert_equal_within(mps1.absdiff.mps2, 0.0d0, 1.0e-10)
  assert_equal(mps1%delete(),Normal)
  assert_equal(mps2%delete(),Normal)
  assert_false(WasThereError())

end test

test tensor3_joinIndices_first

  type(tensor3) :: aMPS
  type(tensor2) :: aMatrix,correct
  complex(8) :: data(2,3,4),matrix(6,4)
  integer error,i,j,k

  !Initialization
  forall (i=1:2 ,j=1:3, k=1:4) data(i,j,k)=ONE*(i+(j-1)*3+(k-1)*4)
  aMPS=new_Tensor(data)

  matrix=one*Reshape( data, [2*3,4])
  correct=new_Tensor(matrix)

  aMatrix=JoinIndices(aMPS,FIRSTANDSECOND,THIRD)
  assert_equal_within(amatrix.absdiff.correct, 0.0d0, 1.0e-8)

  correct=new_Tensor(transpose(matrix)) !This will be 4,2*3
  aMatrix=JoinIndices(aMPS,THIRD,FIRSTANDSECOND)
  assert_equal_within(amatrix.absdiff.correct, 0.0d0, 1.0e-8)

!TODO: New tests with all possible combinations of index bonding
!TODO: Some other tests about creation, products, and so on

   assert_equal(aMPS%delete(),Normal)
   assert_equal(aMatrix%delete(),Normal)
   assert_equal(correct%delete(),Normal)
   assert_false(WasThereError())

end test

!test canonization_of_mps
!  type(MatrixProductState) :: mps
!  integer site
!  real(8) result
!  mps=new_MatrixProductState(10,2,20)
!  assert_false(mps%isCanonized())
!  result=mps%RCanonize()
!  result=mps%LCanonize()
!  assert_equal_within(result, 1.0d0, 1.0e-10)
!  assert_true(mps%isCanonized())
!  result=mps%CanonizeAtSite(5)
!  assert_true(mps%isCanonized())
!  assert_equal_within(result, 5.0d0, 1.0e-10)
!end test



end test_suite
