test_suite MPSTensor_Class
  type(MPSTensor) :: A,B,C
  logical :: Verbose=.false.
  integer error
  integer :: BondL=3
  integer :: BondR=5
  integer :: SpinT=2

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown

test type_creation_deletion
  A=new_MPSTensor(SpinT,BondL,BondR)
  assert_false(WasThereError(verbose))
  assert_false(B%Print().eq.Normal)
  call LowerFlag(verbose)
  B=new_MPSTensor(SpinT,BondL,BondR)
  assert_equal(B%Print(),Normal)
  assert_equal(B%delete(),Normal)
end test

test accesors_work
  A=new_MPSTensor(SpinT,BondL,BondR)
  assert_equal(A%spin(),SpinT)
  assert_equal(A%DLeft(),BondL)
  assert_equal(A%DRight(),BondR)
end test

test Comparison_Of_Dimensions
  A=new_MPSTensor(SpinT,BondL,BondR)
  B=new_MPSTensor(2,4,4)
  assert_false(A.equaldims.B)
  error=B%delete()
  B=new_MPSTensor(SpinT,BondL,BondR)
  assert_true(A.equaldims.B)
end test

end test_suite
