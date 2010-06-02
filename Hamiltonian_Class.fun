test_suite Hamiltonian_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test HamiltonianCreation
  type(Hamiltonian) :: Ham
  integer error
  Ham=new_Hamiltonian()
  error=Ham%delete()
  assert_equal(error,Normal)
  assert_false(WasThereError())
end test

test AddingInteractions
  type(Hamiltonian) :: Ham
  type(Operator) :: 

end test




end test_suite
