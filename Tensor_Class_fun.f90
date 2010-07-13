! Tensor_Class_fun.f90 - a unit test suite for Tensor_Class.f90
!
! funit generated this file from Tensor_Class.fun

module Tensor_Class_fun

 use Tensor_Class

 implicit none

 logical :: noAssertFailed

 public :: test_Tensor_Class

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0



 contains




 subroutine type_creation_deletion


  type(tensor3) :: aTensor
  integer error
  aTensor=new_Tensor(2,20,20)
  error=aTensor%delete()
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(error==Normal)) then
      print *, " *Assert_Equal failed* in test type_creation_deletion &
              &[Tensor_Class.fun:21]"
      print *, "  ", "error (",error,") is not", Normal
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (WasThereError()) then
      print *, " *Assert_False failed* in test type_creation_deletion &
              &[Tensor_Class.fun:22]"
      print *, "  ", "WasThereError() is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine type_creation_deletion


 subroutine assignments_of_tensor


  type(tensor3) :: mps1,mps2
  integer error
  mps1=new_Tensor(10,2,10)
  mps2=mps1
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-10) &
     .ge. &
     (mps1.absdiff.mps2) &
             .and. &
     (0.0d0 &
     -1.0e-10) &
     .le. &
     (mps1.absdiff.mps2) )) then
      print *, " *Assert_Equal_Within failed* in test assignments_of_tensor &
              &[Tensor_Class.fun:32]"
      print *, "  ", "mps1.absdiff.mps2 (",mps1.absdiff.mps2,") is not", &
 0.0d0,"within",1.0e-10
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(mps1%delete()==Normal)) then
      print *, " *Assert_Equal failed* in test assignments_of_tensor &
              &[Tensor_Class.fun:33]"
      print *, "  ", "mps1%delete() (",mps1%delete(),") is not", Normal
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(mps2%delete()==Normal)) then
      print *, " *Assert_Equal failed* in test assignments_of_tensor &
              &[Tensor_Class.fun:34]"
      print *, "  ", "mps2%delete() (",mps2%delete(),") is not", Normal
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (WasThereError()) then
      print *, " *Assert_False failed* in test assignments_of_tensor &
              &[Tensor_Class.fun:35]"
      print *, "  ", "WasThereError() is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine assignments_of_tensor


 subroutine tensor3_joinIndices_first


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
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (amatrix.absdiff.correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (amatrix.absdiff.correct) )) then
      print *, " *Assert_Equal_Within failed* in test tensor3_joinIndices_first &
              &[Tensor_Class.fun:54]"
      print *, "  ", "amatrix.absdiff.correct (",amatrix.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  correct=new_Tensor(transpose(matrix)) !This will be 4,2*3
  aMatrix=JoinIndices(aMPS,THIRD,FIRSTANDSECOND)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (amatrix.absdiff.correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (amatrix.absdiff.correct) )) then
      print *, " *Assert_Equal_Within failed* in test tensor3_joinIndices_first &
              &[Tensor_Class.fun:58]"
      print *, "  ", "amatrix.absdiff.correct (",amatrix.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(aMPS%delete()==Normal)) then
      print *, " *Assert_Equal failed* in test tensor3_joinIndices_first &
              &[Tensor_Class.fun:61]"
      print *, "  ", "aMPS%delete() (",aMPS%delete(),") is not", Normal
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(aMatrix%delete()==Normal)) then
      print *, " *Assert_Equal failed* in test tensor3_joinIndices_first &
              &[Tensor_Class.fun:62]"
      print *, "  ", "aMatrix%delete() (",aMatrix%delete(),") is not", Normal
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(correct%delete()==Normal)) then
      print *, " *Assert_Equal failed* in test tensor3_joinIndices_first &
              &[Tensor_Class.fun:63]"
      print *, "  ", "correct%delete() (",correct%delete(),") is not", Normal
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (WasThereError()) then
      print *, " *Assert_False failed* in test tensor3_joinIndices_first &
              &[Tensor_Class.fun:64]"
      print *, "  ", "WasThereError() is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine tensor3_joinIndices_first


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

 subroutine funit_setup
  !use ErrorHandling
  !use Constants
  !Set testing mode
  MaxErrorAllowed=CriticalError
  noAssertFailed = .true.
 end subroutine funit_setup


 subroutine funit_teardown

 end subroutine funit_teardown


 subroutine test_Tensor_Class( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call funit_setup
  call type_creation_deletion
  call funit_teardown

  call funit_setup
  call assignments_of_tensor
  call funit_teardown

  call funit_setup
  call tensor3_joinIndices_first
  call funit_teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_Tensor_Class

end module Tensor_Class_fun
