program TestSuite

!  use ErrorHandling
  use Tensor
  use TestHelper

  implicit none

  integer error
  type(MPSTensor) A
  integer,parameter :: SpinTest=2,BondTestL=10,BondTestR=20

  !Try to catch .print. with an non-initialized tensor
  error=.print.A
  if(error.eq.Normal) then !This fails so error code is not Normal
     call ThrowException('Test Suite','Tensor printing did not catch uninitialized tensor',error,Warning)
     call PrintTestMsg('PrintMPSTensor::Uninitialized tensor',FAILED)
  else
     call PrintTestMsg('PrintMPSTensor::Uninitialized tensor',PASSED)
  endif

  !Initialize tensor
  A = new_MPSTensor(2,2,3)

  !Now try and print tensor
  error=.print.A
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','Tensor printing failed',error,Warning)
     call PrintTestMsg('.print.tensor ',FAILED)
  else
     call PrintTestMsg('.print.tensor ',PASSED)
  endif
  
  !Try and delete the tensor now
  error=delete_MPSTensor(A)
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','Could not delete tensor',error,Warning)
     call PrintTestMsg('Delete tensor',FAILED)
  else
     call PrintTestMsg('Delete tensor',PASSED) 
 endif

  !Initialize tensor again
  A = new_MPSTensor(SpinTest,BondTestL,BondTestR)

  !Test accessor methods to properties
  error=(Spin(A).ne.SpinTest)
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','Spin dimension is wrong',error,Warning)
     call PrintTestMsg('Spin Accessor',FAILED) 
  else
     call PrintTestMsg('Spin Accessor',PASSED) 
  endif
  error=(DLeft(A).ne.BondTestL)
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','DLeft dimension is wrong',error,Warning)
     call PrintTestMsg('DLeft Accessor',FAILED) 
  else
     call PrintTestMsg('DLeft Accessor',PASSED) 
  endif  
  error=(DRight(A).ne.BondTestR)
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','DRight dimension is wrong',error,Warning)
     call PrintTestMsg('DRight Accessor',FAILED) 
  else
     call PrintTestMsg('DRight Accessor',PASSED) 
  endif 
  


! ##### END TESTING, FINAL REPORT


  print *,'    ###########################   '
  print *,'####      SUMMARY OF TESTS     ####'
  print *,'    ###########################   '
  print *,'Tests PASSED :',testsOK,' out of ',testsOK+testsNOT



  end program TestSuite

