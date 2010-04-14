program TestSuite
  use ErrorHandling
  use Tensor

  implicit none

  integer error=0

  call test_tensor(error)

  if(error.ne.0) then
     call ThrowException('Test Suite','test_tensor failed',error,1)
  else
     print *,'test_tensor : PASSED '
  endif

  contains
    subroutine test_tensor(error)
      integer error
      type(MPSTensor) A

      A = new_MPSTensor(2,10,10)
      print *,A(1:2,1:10,1:10)

    end subroutine test_tensor

  end program TestSuite
