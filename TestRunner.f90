
! TestRunner.f90 - runs fUnit test suites
!
! funit generated this file on Tue Jul 13 00:56:41 +0200 2010.

program TestRunner

    use Tensor_Class_fun
  
  implicit none

  integer, dimension(1) :: numTests, numAsserts, numAssertsTested, numFailures

    write(*,*)
  write(*,*) "Tensor_Class test suite:"
  call test_Tensor_Class &
    ( numTests(1), numAsserts(1), numAssertsTested(1), numFailures(1) )
  write(*,1) numAssertsTested(1), numAsserts(1), &
    numTests(1)-numFailures(1), numTests(1)
1 format('Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')
  
  write(*,*)
  write(*,'(a)') "==========[ SUMMARY ]=========="
      write(*,'(a14)',advance="no") " Tensor_Class:"
  if ( numFailures(1) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
