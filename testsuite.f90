program TestSuite

!  use ErrorHandling
  use Tensor

  implicit none

  integer error
  type(MPSTensor) A
  integer,parameter :: SpinTest=2,BondTestL=10,BondTestR=20
  integer :: testsOK=0,testsNOT=0

  !Try to catch .print. with an non-initialized tensor
  error=.print.A
  if(error.eq.Normal) then !This fails so error code is not Normal
     call ThrowException('Test Suite','Tensor printing did not catch uninitialized tensor',error,Warning)
     testsNOT=testsNOT+1
  else
     print *,'PrintMPSTensor::Uninitialized tensor = PASSED '
     testsOK=testsOK+1
  endif

  !Initialize tensor
  A = new_MPSTensor(2,2,3)

  !Now try and print tensor
  error=.print.A
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','Tensor printing failed',error,Warning)
     testsNOT=testsNOT+1
  else
     print *,'.print.tensor                       = PASSED '
     testsOK=testsOK+1
  endif
  
  !Try and delete the tensor now
  error=delete_MPSTensor(A)
  if(error.ne.Normal) then !No failure so error = 0
     call ThrowException('Test Suite','Could not delete tensor',error,Warning)
     testsNOT=testsNOT+1
  else
     print *,'Delete tensor                       = PASSED '
     testsOK=testsOK+1
  endif

  !Initialize tensor
  A = new_MPSTensor(2,2,3)
  


! ##### END TESTING, FINAL REPORT


  print *,'    ###########################   '
  print *,'####      SUMMARY OF TESTS     ####'
  print *,'    ###########################   '
  print *,'Tests PASSED :',testsOK,' out of ',testsOK+testsNOT

  end program TestSuite
