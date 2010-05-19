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
