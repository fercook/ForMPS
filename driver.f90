! Tensor_Class_fun.f90 - a unit test suite for Tensor_Class.f90
!
! funit generated this file from Tensor_Class.fun

module Tensor_Class_fun

 use Tensor_Class
 use CommonModules

!This is the NEWNEW commend  Tensor_Class

implicit none
		      
logical :: noAssertFailed

public :: test_Tensor_Class

private

integer :: numTests          = 0
integer :: numAsserts        = 0
integer :: numAssertsTested  = 0
integer :: numFailures       = 0



 ! And before the contains goes...
 contains

!TODO: New tests with all possible combinations of index bonding

! use ErrorHandling
! use Constants




 subroutine type_creation_deletion


  type(tensor3) :: aTensor
  integer error
  aTensor=new_Tensor(2,20,20)
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (WasThereError()) then
      print *, " *Assert_False failed* in test type_creation_deletion &
              &[Tensor_Class.fun:24]"
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
              &[Tensor_Class.fun:34]"
      print *, "  ", "mps1.absdiff.mps2 (",mps1.absdiff.mps2,") is not", &
 0.0d0,"within",1.0e-10
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


 subroutine tensor3_joinindices_first

  type(tensor3) :: aMPS
  type(tensor2) :: aMatrix,correct
  complex(8) :: data(2,3,4),matrix(6,4)
  integer error,i,j,k

  !initialization
  forall (i=1:2 ,j=1:3, k=1:4) data(i,j,k)=ONE*(i+(j-1)*3+(k-1)*4)
  aMPS=new_Tensor(data)

  matrix=one*Reshape( data, [2*3,4])
  correct=new_Tensor(matrix)

  aMatrix=aMPS%Joinindices(FiRSTANDSECOND,THiRD)
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
      print *, " *Assert_Equal_Within failed* in test tensor3_joinindices_first &
              &[Tensor_Class.fun:53]"
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
  aMatrix=JoinindicesOf(aMPS,THiRD,FiRSTANDSECOND)
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
      print *, " *Assert_Equal_Within failed* in test tensor3_joinindices_first &
              &[Tensor_Class.fun:57]"
      print *, "  ", "amatrix.absdiff.correct (",amatrix.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
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
      print *, " *Assert_False failed* in test tensor3_joinindices_first &
              &[Tensor_Class.fun:59]"
      print *, "  ", "WasThereError() is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine tensor3_joinindices_first


 subroutine tensor4_joinindices

  type(tensor4) :: aTensor
  type(tensor2) :: aMatrix,correct
  complex(8) :: data(2,3,4,5),matrix(6,20)
  integer error,i,j,k,l

  !initialization
  forall (i=1:2 ,j=1:3, k=1:4, l=1:5) data(i,j,k,l)=ONE*(i+(j-1)*2+(k-1)*3*2+(l-1)*2*3*4)
  aTensor=new_Tensor(data)

  matrix=one*Reshape( data, [2*3,4*5])
  correct=new_Tensor(matrix)

  aMatrix=aTensor%Joinindices(FiRSTANDSECOND,THiRDANDFOURTH)
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
      print *, " *Assert_Equal_Within failed* in test tensor4_joinindices &
              &[Tensor_Class.fun:76]"
      print *, "  ", "amatrix.absdiff.correct (",amatrix.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  correct=new_Tensor(transpose(matrix))

  aMatrix=aTensor%Joinindices(THiRDANDFOURTH,FiRSTANDSECOND)

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
      print *, " *Assert_Equal_Within failed* in test tensor4_joinindices &
              &[Tensor_Class.fun:82]"
      print *, "  ", "amatrix.absdiff.correct (",amatrix.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine tensor4_joinindices



 subroutine tensor2_Splitindex

  type(tensor3) :: aMPS,correct
  type(tensor2) :: aMatrix
  complex(8) :: matrixdata(6,4)
  integer error,i,j,k

  !initialization
  forall (i=1:6 ,j=1:4) matrixdata(i,j)=ONE*(i+(j-1)*6)
  aMatrix=new_Tensor(matrixdata)

  correct=new_Tensor(one*Reshape( matrixdata, [2,3,4]))

  aMPS=SplitindexOf(aMatrix,FIRST,2)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (aMPS.absdiff.correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (aMPS.absdiff.correct) )) then
      print *, " *Assert_Equal_Within failed* in test tensor2_Splitindex &
              &[Tensor_Class.fun:100]"
      print *, "  ", "aMPS.absdiff.correct (",aMPS.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  correct=new_Tensor(one*Reshape( matrixdata, [6,2,2]))

  aMPS=SplitindexOf(aMatrix,SECOND,2)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (aMPS.absdiff.correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (aMPS.absdiff.correct) )) then
      print *, " *Assert_Equal_Within failed* in test tensor2_Splitindex &
              &[Tensor_Class.fun:105]"
      print *, "  ", "aMPS.absdiff.correct (",aMPS.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
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
      print *, " *Assert_False failed* in test tensor2_Splitindex &
              &[Tensor_Class.fun:107]"
      print *, "  ", "WasThereError() is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine tensor2_Splitindex


 subroutine transpose_of_tensor

    integer,parameter :: d1=2,d2=3,d3=4
    type(Tensor3) :: aTensor,bTensor,CorrectTensor
    complex(8) :: data(d1,d2,d3)
    complex(8),allocatable :: correct(:,:,:)
    integer :: i,j,k,neworder(3)

    forall (i=1:d1, j=1:d2, k=1:d3) data(i,j,k)=one*(i*100+j*10+II*k)
    aTensor=new_Tensor(data)

    neworder=[3,1,2]
    bTensor=ConjugateTranspose(aTensor,neworder)
    allocate(correct( d2,d3,d1 ))
    forall (i=1:d1, j=1:d2, k=1:d3) correct(j,k,i)=one*(i*100+j*10+II*k)
    CorrectTensor=new_Tensor(dconjg(correct))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (bTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (bTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test transpose_of_tensor &
              &[Tensor_Class.fun:125]"
      print *, "  ", "bTensor.absdiff.correctTensor (",bTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(correct)

    neworder=[2,1,3]
    bTensor=ConjugateTranspose(aTensor,neworder)
    allocate(correct( d2,d1,d3 ))
    forall (i=1:d1, j=1:d2, k=1:d3) correct(j,i,k)=one*(i*100+j*10+II*k)
    CorrectTensor=new_Tensor(dconjg(correct))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (bTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (bTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test transpose_of_tensor &
              &[Tensor_Class.fun:133]"
      print *, "  ", "bTensor.absdiff.correctTensor (",bTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(correct)

    neworder=[1,3,2]
    bTensor=ConjugateTranspose(aTensor,neworder)
    allocate(correct( d1,d3,d2 ))
    forall (i=1:d1, j=1:d2, k=1:d3) correct(i,k,j)=one*(i*100+j*10+II*k)
    CorrectTensor=new_Tensor(dconjg(correct))
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (bTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (bTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test transpose_of_tensor &
              &[Tensor_Class.fun:141]"
      print *, "  ", "bTensor.absdiff.correctTensor (",bTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(correct)


  numTests = numTests + 1

 end subroutine transpose_of_tensor


 subroutine transpose_of_tensor4

    integer,parameter :: d1=2,d2=3,d3=4,d4=5
    type(Tensor4) :: aTensor,bTensor,CorrectTensor
    complex(8) :: data(d1,d2,d3,d4)
    complex(8),allocatable :: correct(:,:,:,:)
    integer :: i,j,k,l,neworder(4)

    forall (i=1:d1, j=1:d2, k=1:d3, l=1:d4) data(i,j,k,l)=one*(i*1000+j*100+10*II*k+l)
    aTensor=new_Tensor(data)

    neworder=[1,2,3,4]
    bTensor=TensorTranspose(aTensor,neworder)
    allocate(correct(d1,d2,d3,d4))
    forall (i=1:d1, j=1:d2, k=1:d3, l=1:d4) correct(i,j,k,l)=one*(i*1000+j*100+10*II*k+l)
    CorrectTensor=new_Tensor(correct)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (bTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (bTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test transpose_of_tensor4 &
              &[Tensor_Class.fun:161]"
      print *, "  ", "bTensor.absdiff.correctTensor (",bTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(correct)

    neworder=[4,1,2,3]
    bTensor=TensorTranspose(aTensor,neworder)
    allocate(correct(d2,d3,d4,d1))
    forall (i=1:d1, j=1:d2, k=1:d3, l=1:d4) correct(j,k,l,i)=one*(i*1000+j*100+10*II*k+l)
    CorrectTensor=new_Tensor(correct)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (bTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (bTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test transpose_of_tensor4 &
              &[Tensor_Class.fun:169]"
      print *, "  ", "bTensor.absdiff.correctTensor (",bTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(correct)

    neworder=[3,2,4,1]
    bTensor=TensorTranspose(aTensor,neworder)
    allocate(correct(d4,d2,d1,d3))
    forall (i=1:d1, j=1:d2, k=1:d3, l=1:d4) correct(l,j,i,k)=one*(i*1000+j*100+10*II*k+l)
    CorrectTensor=new_Tensor(correct)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (bTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (bTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test transpose_of_tensor4 &
              &[Tensor_Class.fun:177]"
      print *, "  ", "bTensor.absdiff.correctTensor (",bTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(correct)


  numTests = numTests + 1

 end subroutine transpose_of_tensor4



 subroutine unfoldings_of_tensor3

    integer,parameter :: d1=2,d2=3,d3=4
    type(Tensor3) :: aTensor
    type(Tensor2) :: unfoldedTensor,CorrectTensor
    complex(8) :: data(d1,d2,d3)
    complex(8),allocatable :: matrixform(:,:)
    integer :: i1,i2,i3

    forall (i1=1:d1, i2=1:d2, i3=1:d3) data(i1,i2,i3)=one*(i1*100+i2*10+II*i3)
    aTensor=new_Tensor(data)

    unfoldedTensor=aTensor.unfold.1
    allocate(matrixform( d1,d2*d3 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3) matrixform (i1,(i3-1)*d2+i2)=one*(i1*100+i2*10+II*i3)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor3 &
              &[Tensor_Class.fun:198]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)

    unfoldedTensor=aTensor.unfold.2
    allocate(matrixform( d2,d1*d3 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3) matrixform (i2,(i3-1)*d1+i1)=one*(i1*100+i2*10+II*i3)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor3 &
              &[Tensor_Class.fun:205]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)

    unfoldedTensor=aTensor.unfold.3
    allocate(matrixform( d3,d1*d2 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3) matrixform (i3,(i1-1)*d2+i2)=one*(i1*100+i2*10+II*i3)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor3 &
              &[Tensor_Class.fun:212]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)


  numTests = numTests + 1

 end subroutine unfoldings_of_tensor3


 subroutine unfoldings_of_tensor4

    integer,parameter :: d1=2,d2=3,d3=4,d4=3
    type(Tensor4) :: aTensor
    type(Tensor2) :: unfoldedTensor,CorrectTensor
    complex(8) :: data(d1,d2,d3,d4)
    complex(8),allocatable :: matrixform(:,:)
    integer :: i1,i2,i3,i4

    forall (i1=1:d1, i2=1:d2, i3=1:d3, i4=1:d4) data(i1,i2,i3,i4)=one*(i1*1000+i2*100+II*10*i3+i4)
    aTensor=new_Tensor(data)

    unfoldedTensor=aTensor.unfold.1
    allocate(matrixform( d1,d2*d3*d4 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3,i4=1:d4) matrixform (i1,(i4-1)*d2*d3+(i3-1)*d2+i2)=one*(i1*1000+i2*100+II*10*i3+i4)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor4 &
              &[Tensor_Class.fun:232]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)

    unfoldedTensor=aTensor.unfold.2
    allocate(matrixform( d2,d1*d3*d4 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3,i4=1:d4) matrixform (i2,(i4-1)*d1*d3+(i3-1)*d1+i1)=one*(i1*1000+i2*100+II*10*i3+i4)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor4 &
              &[Tensor_Class.fun:239]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)

    unfoldedTensor=aTensor.unfold.3
    allocate(matrixform( d3,d2*d1*d4 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3,i4=1:d4) matrixform (i3,(i4-1)*d2*d1+(i1-1)*d2+i2)=one*(i1*1000+i2*100+II*10*i3+i4)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor4 &
              &[Tensor_Class.fun:246]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)

    unfoldedTensor=aTensor.unfold.4
    allocate(matrixform( d4,d2*d3*d1 ))
    forall (i1=1:d1, i2=1:d2, i3=1:d3,i4=1:d4) matrixform (i4,(i1-1)*d2*d3+(i3-1)*d2+i2)=one*(i1*1000+i2*100+II*10*i3+i4)
    CorrectTensor=new_Tensor(matrixform)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (unfoldedTensor.absdiff.correctTensor) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (unfoldedTensor.absdiff.correctTensor) )) then
      print *, " *Assert_Equal_Within failed* in test unfoldings_of_tensor4 &
              &[Tensor_Class.fun:253]"
      print *, "  ", "unfoldedTensor.absdiff.correctTensor (",unfoldedTensor.absdiff.correctTensor,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
    deallocate(matrixform)


  numTests = numTests + 1

 end subroutine unfoldings_of_tensor4




 subroutine matrix_product

    integer,parameter :: DLeft=6,DRight=8
    type(Tensor2) :: mat,mat1, mat2,correct
    complex(8) :: data(DLeft,DLeft),data1(DLeft,DRight),data2(DRight,DLeft)
    integer i,j

    forall (i=1:Dleft ,j=1:Dright) &
      & data1(i,j)=one*(exp(II*i*Pi/DLeft)+(j-1)*Dleft)
    forall (j=1:Dleft ,i=1:Dright) &
      & data2(i,j)=one*(exp(II*i*Pi/100)+(j-1)*DRight*6)
    data=matmul(data1,  data2)

    correct=new_Tensor(data)
    mat1=new_Tensor(data1)
    mat2=new_Tensor(data2)
    mat=mat1*mat2

  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-10) &
     .ge. &
     (mat .absdiff. correct) &
             .and. &
     (0.0d0 &
     -1.0e-10) &
     .le. &
     (mat .absdiff. correct) )) then
      print *, " *Assert_Equal_Within failed* in test matrix_product &
              &[Tensor_Class.fun:323]"
      print *, "  ", "mat .absdiff. correct (",mat .absdiff. correct,") is not", &
 0.0d0,"within",1.0e-10
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine matrix_product



 subroutine Singular_Value_Decomposition

   integer,parameter :: LeftDimension=6,RightDimension=8
   type(Tensor2) :: aMatrix,theU,theVt
   type(Tensor2) :: theSigma
   complex(8) :: data(LeftDimension,RightDimension)
   integer :: i,j,k

   !Input value is somewhat regular, perhaps try with random data at some point
   forall (i=1:LeftDimension ,j=1:RightDimension) &
      & data(i,j)=one*(exp(II*i*Pi/LeftDimension)+(j-1)*LeftDimension)
      !Input value is somewhat regular, perhaps try with random data at some point

    aMatrix=new_Tensor(data)

   call aMatrix%SVD(theU,theSigma,theVt)
   !Now test if the three output matrices form the original one
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-10) &
     .ge. &
     (aMatrix.absdiff.(theU*(theSigma*theVt))) &
             .and. &
     (0.0d0 &
     -1.0e-10) &
     .le. &
     (aMatrix.absdiff.(theU*(theSigma*theVt))) )) then
      print *, " *Assert_Equal_Within failed* in test Singular_Value_Decomposition &
              &[Tensor_Class.fun:343]"
      print *, "  ", "aMatrix.absdiff.(theU*(theSigma*theVt)) (",aMatrix.absdiff.(theU*(theSigma*theVt)),") is not", &
 0.0d0,"within",1.0e-10
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
      print *, " *Assert_False failed* in test Singular_Value_Decomposition &
              &[Tensor_Class.fun:344]"
      print *, "  ", "WasThereError() is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine Singular_Value_Decomposition



 subroutine Right_Compactification


  type(Tensor3) :: T_Tensor
  type(Tensor2) :: T_Compacted,T_Matrix,T_Correct
  integer,parameter :: DleftT=4, DrightT=3,spinT=2
  complex(8) :: data(DleftT,DrightT,spinT)
  complex(8) :: matrix(DrightT,DrightT)
  complex(8) :: CorrectResult(DleftT,DleftT)
  integer n,i,j,k,s

  do i=1,DleftT
     do j=1,DrightT
        do s=1,spinT
           data(i,j,s)=one*(i+(j-1)*DleftT+(s-1)*DrightT)
        enddo
     enddo
  enddo
  do i=1,DrightT
     do j=1,DrightT
        matrix(i,j)=(ii**i+(j-1)*DrightT)
     enddo
  enddo

  CorrectResult=one*reshape([3072 - 312*ii, 3528 - 312*ii, 3984 - 312*ii, 4440 - 312*ii, 3384 - &
       & 360*ii, 3888 - 360*ii, 4392 - 360*ii, 4896 - 360*ii, 3696 - &
       & 408*ii, 4248 - 408*ii, 4800 - 408*ii, 5352 - 408*ii, 4008 - &
       & 456*ii, 4608 - 456*ii, 5208 - 456*ii, 5808 - 456*ii], [DleftT,DleftT] )
  T_Tensor=new_Tensor(data)
  T_Matrix=new_Tensor(matrix)
  T_Correct=new_Tensor(CorrectResult)
  T_Compacted=CompactRight(T_Matrix,T_Tensor,Conjugate(T_Tensor),THiRD)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (T_Compacted.absdiff.T_Correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (T_Compacted.absdiff.T_Correct) )) then
      print *, " *Assert_Equal_Within failed* in test Right_Compactification &
              &[Tensor_Class.fun:485]"
      print *, "  ", "T_Compacted.absdiff.T_Correct (",T_Compacted.absdiff.T_Correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine Right_Compactification


 subroutine Left_Compactification

!
!  Mathematica code: NOTiCE THE TRANSPOSE TO GET THE ORDER RiGHT
! With[{DL = 3, DR = 4, spin = 2},
!   At = Table[
!     DR*(s - 1) + (j - 1)*DL + i, {s, 1, 2}, {i, 1, DL}, {j, 1, DR}];
!   mat = Table[i^i + (j - 1)*DL, {i, 1, DL}, {j, 1, DL}]];
!  Flatten[Transpose[LProduct[At, At, mat]]]
!
  type(Tensor3) :: T_Tensor
  type(Tensor2) :: T_Compacted,T_Matrix,T_Correct
  integer,parameter :: DleftT=3, DrightT=4,spinT=2
  complex(8) :: matrix(DleftT,DleftT)
  complex(8) :: CorrectResult(DrightT,DrightT)
  complex(8) :: data(DleftT,DrightT,spinT)
  integer n,i,j,k,s

   do i=1,DleftT
     do j=1,DrightT
        do s=1,spinT
           data(i,j,s)=one*(i+(j-1)*DleftT+(s-1)*DrightT)
        enddo
     enddo
  enddo
  do i=1,DleftT
     do j=1,DleftT
        matrix(i,j)=(ii**i+(j-1)*DleftT)
     enddo
  enddo

  CorrectResult=one*reshape([1104 - 48*ii, 1788 - 48*ii, 2472 - 48*ii, 3156 - 48*ii, 1680 - &
       &  84*ii, 2796 - 84*ii, 3912 - 84*ii, 5028 - 84*ii, 2256 - 120*ii, 3804 - &
       &  120*ii, 5352 - 120*ii, 6900 - 120*ii, 2832 - 156*ii, 4812 - &
       &  156*ii, 6792 - 156*ii, 8772 - 156*ii], [DrightT,DrightT] )

  T_Tensor=new_Tensor(data)
  T_Matrix=new_Tensor(matrix)
  T_Correct=new_Tensor(CorrectResult)
  T_Compacted=CompactLeft(T_Matrix,T_Tensor,Conjugate(T_Tensor),THiRD)
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (T_Compacted.absdiff.T_Correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (T_Compacted.absdiff.T_Correct) )) then
      print *, " *Assert_Equal_Within failed* in test Left_Compactification &
              &[Tensor_Class.fun:528]"
      print *, "  ", "T_Compacted.absdiff.T_Correct (",T_Compacted.absdiff.T_Correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine Left_Compactification


 subroutine Compact_From_Below_T3_T4

!Mathematica Code:
!
!T3 = Table[a*100 + b*10 + 1.0 c, {a, 1, 2}, {b, 1, 3}, {c, 1, 4}];
!T4 = Table[a*1000 + b*100 + c 10.0 + d I, {a, 1, 2}, {b, 1, 2}, {c, 1, 2}, {d,1, 2}];
!math = Round[
!   Flatten[Transpose[T3, {3, 1, 2}].Transpose[
!      T4, {2, 3, 1, 4}], {{5}, {3, 1}, {4, 2}}]];
!Flatten[Transpose[math, {3, 2, 1}]]
  type(Tensor3) :: aT3,correct,result
  type(Tensor4) :: aT4
  complex(8) :: origArray(2,3,4), origTensor(2,2,2,2), correctArray(2,6,8)
  integer :: i,j,k,l

  forall (i=1:2, j=1:3, k=1:4) origArray(i,j,k)=100*i+10*j+k
  forall (i=1:2, j=1:2, k=1:2, l=1:2) origTensor(i,j,k,l)=1000*i+100*j+10*k+ii*l

  correctArray= ONE*reshape( [359530 + 322 *II, 359530 + 644 *II, 381830 + 342 *II, 381830 +  &
	&	 684 *II, 404130 + 362 *II, 404130 + 724 *II, 681530 + 322 *II, 681530 + &
	&	 644 *II, 723830 + 342 *II, 723830 + 684 *II, 766130 + 362 *II, 766130 + &
	&	 724 *II, 361760 + 324 *II, 361760 + 648 *II, 384060 + 344 *II, 384060 + &
	&	 688 *II, 406360 + 364 *II, 406360 + 728 *II, 685760 + 324 *II, 685760 + &
	&	 648 *II, 728060 + 344 *II, 728060 + 688 *II, 770360 + 364 *II, 770360 + &
	&	 728 *II, 363990 + 326 *II, 363990 + 652 *II, 386290 + 346 *II, 386290 + &
	&	 692 *II, 408590 + 366 *II, 408590 + 732 *II, 689990 + 326 *II, 689990 + &
	&	 652 *II, 732290 + 346 *II, 732290 + 692 *II, 774590 + 366 *II, 774590 + &
	&	 732 *II, 366220 + 328 *II, 366220 + 656 *II, 388520 + 348 *II, 388520 + &
	&	 696 *II, 410820 + 368 *II, 410820 + 736 *II, 694220 + 328 *II, 694220 + &
	&	 656 *II, 736520 + 348 *II, 736520 + 696 *II, 778820 + 368 *II, 778820 + &
	&	 736 *II, 391730 + 322 *II, 391730 + 644 *II, 416030 + 342 *II, 416030 + &
	&	 684 *II, 440330 + 362 *II, 440330 + 724 *II, 713730 + 322 *II, 713730 + &
	&	 644 *II, 758030 + 342 *II, 758030 + 684 *II, 802330 + 362 *II, 802330 + &
	&	 724 *II, 394160 + 324 *II, 394160 + 648 *II, 418460 + 344 *II, 418460 + &
	&	 688 *II, 442760 + 364 *II, 442760 + 728 *II, 718160 + 324 *II, 718160 + &
	&	 648 *II, 762460 + 344 *II, 762460 + 688 *II, 806760 + 364 *II, 806760 + &
	&	 728 *II, 396590 + 326 *II, 396590 + 652 *II, 420890 + 346 *II, 420890 + &
	&	 692 *II, 445190 + 366 *II, 445190 + 732 *II, 722590 + 326 *II, 722590 + &
	&	 652 *II, 766890 + 346 *II, 766890 + 692 *II, 811190 + 366 *II, 811190 + &
	&	 732 *II, 399020 + 328 *II, 399020 + 656 *II, 423320 + 348 *II, 423320 + &
	&	 696 *II, 447620 + 368 *II, 447620 + 736 *II, 727020 + 328 *II, 727020 + &
	&	 656 *II, 771320 + 348 *II, 771320 + 696 *II, 815620 + 368 *II, 815620 + 736 *II ], [2,6,8] )

    aT3=new_Tensor(origArray)
    aT4=new_Tensor(origTensor)
    correct=new_Tensor(correctArray)
    result=CompactBelow(aT3,FIRST,aT4,THIRD,FOURTH)

  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (result.absdiff.correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (result.absdiff.correct) )) then
      print *, " *Assert_Equal_Within failed* in test Compact_From_Below_T3_T4 &
              &[Tensor_Class.fun:579]"
      print *, "  ", "result.absdiff.correct (",result.absdiff.correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine Compact_From_Below_T3_T4





 subroutine Matrices_times_tensor3s


  type(Tensor3) :: T_result,T_Tensor,T_Correct
  type(Tensor2) :: T_Matrix
  integer,parameter :: left=2,center=3,right=4
  complex(8) :: matrixL(center,left),matrixR(right,left)
  complex(8) :: CorrectL(center,center,right),correctR(left,center,left)
  complex(8) :: tensorarray(left,center,right)
  integer n,i,j,k,s

  forall (i=1:left, j=1:center, k=1:right) tensorarray(i,j,k)=1000*i+100*j*II+10*k
  forall (i=1:center, j=1:left) matrixL(i,j)=127*i+13*j*II
  forall (i=1:right, j=1:left) matrixR(i,j)=41*i+7*j*II
  do k=1,right
    do j=1,center
        do i=1,center
            correctL(i,j,k)=sum( matrixL(i,:)*tensorarray(:,j,k) )!dot_product(dconjg(matrixL(i,:)),tensorarray(:,j,k))
        enddo
    enddo
  enddo
  T_Tensor=new_Tensor(tensorarray)
  T_matrix=new_Tensor(matrixL)
  T_result=T_matrix*T_Tensor
  T_correct=new_Tensor(correctL)

  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (T_result.absdiff.T_correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (T_result.absdiff.T_correct) )) then
      print *, " *Assert_Equal_Within failed* in test Matrices_times_tensor3s &
              &[Tensor_Class.fun:691]"
      print *, "  ", "T_result.absdiff.T_correct (",T_result.absdiff.T_correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  do k=1,left
    do j=1,center
        do i=1,left
            correctR(i,j,k)=sum( tensorarray(i,j,:) * matrixR(:,k) )
        enddo
    enddo
  enddo
  T_matrix=new_Tensor(matrixR)
  T_result=T_Tensor*T_matrix
  T_correct=new_Tensor(correctR)

  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-8) &
     .ge. &
     (T_result.absdiff.T_correct) &
             .and. &
     (0.0d0 &
     -1.0e-8) &
     .le. &
     (T_result.absdiff.T_correct) )) then
      print *, " *Assert_Equal_Within failed* in test Matrices_times_tensor3s &
              &[Tensor_Class.fun:704]"
      print *, "  ", "T_result.absdiff.T_correct (",T_result.absdiff.T_correct,") is not", &
 0.0d0,"within",1.0e-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine Matrices_times_tensor3s



 subroutine xplusWithTensor4


    type(tensor4) :: aT4,anotherT4,newT4,correctT4
    integer,parameter :: lft=4,rgt=5,ump=2,dwn=3
    real(8) :: data1R(lft,rgt),data1I(lft,rgt),data2R(rgt,lft),data2I(rgt,lft)
    complex(8) :: fulltensor1(lft,rgt,ump,dwn), fulltensor2(rgt,lft,ump,dwn),correctTensor(lft,lft,ump*ump,dwn*dwn)
    integer :: i,j,k,l

    do i=1,ump
      do j=1,dwn
        call random_number(data1R)
        call random_number(data1I)
        call random_number(data2R)
        call random_number(data2I)
        fulltensor1(:,:,i,j)=data1R+II*data1I
        fulltensor2(:,:,i,j)=data2R+II*data2I
      enddo
    enddo
    correctTensor=ZERO
    do i=1,ump
    do j=1,ump
      do k=1,dwn
      do l=1,dwn
        correctTensor(:,:,(j-1)*ump+i,(l-1)*dwn+k)=matmul(fulltensor1(:,:,i,k),fulltensor2(:,:,j,l))
      enddo
      enddo
    enddo
    enddo

    correctT4=new_Tensor(correctTensor)
    aT4=new_Tensor(fulltensor1)
    anotherT4=new_Tensor(fulltensor2)
    newT4=aT4.xplus.anotherT4

  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-10) &
     .ge. &
     (newT4.absdiff.correctT4) &
             .and. &
     (0.0d0 &
     -1.0e-10) &
     .le. &
     (newT4.absdiff.correctT4) )) then
      print *, " *Assert_Equal_Within failed* in test xplusWithTensor4 &
              &[Tensor_Class.fun:757]"
      print *, "  ", "newT4.absdiff.correctT4 (",newT4.absdiff.correctT4,") is not", &
 0.0d0,"within",1.0e-10
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine xplusWithTensor4


 subroutine TensorTraceWithTensor4


    type(tensor4) :: aT4
    type(tensor2) :: aT2,correctT2
    integer,parameter :: lft=5,rgt=5,ump=2,dwn=3
    real(8) :: dataR(lft,rgt),dataI(lft,rgt)
    complex(8) :: fulltensor(lft,rgt,ump,dwn), correctTensor(ump,dwn)
    integer :: i,j,k,l

    do i=1,ump
      do j=1,dwn
        call random_number(dataR)
        call random_number(dataI)
        fulltensor(:,:,i,j)=dataR+II*dataI
      enddo
    enddo
    correctTensor=ZERO
    do i=1,ump
      do j=1,dwn
        do k=1,lft
            correctTensor(i,j)=correctTensor(i,j)+ (fulltensor(k,k,i,j))
        enddo
      enddo
    enddo

    correctT2=new_Tensor(correctTensor)
    aT4=new_Tensor(fulltensor)
    aT2=TensorTrace(aT4,THIRDANDFOURTH)
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.(aT2%GetDimensions().equalvector.[ump,dwn])) then
      print *, " *Assert_True failed* in test TensorTraceWithTensor4 &
              &[Tensor_Class.fun:789]"
      print *, "  ", "aT2%GetDimensions().equalvector.[ump,dwn] is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0e-10) &
     .ge. &
     (aT2.absdiff.correctT2) &
             .and. &
     (0.0d0 &
     -1.0e-10) &
     .le. &
     (aT2.absdiff.correctT2) )) then
      print *, " *Assert_Equal_Within failed* in test TensorTraceWithTensor4 &
              &[Tensor_Class.fun:790]"
      print *, "  ", "aT2.absdiff.correctT2 (",aT2.absdiff.correctT2,") is not", &
 0.0d0,"within",1.0e-10
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif


  numTests = numTests + 1

 end subroutine TensorTraceWithTensor4


 subroutine LinearSolver

    integer,parameter :: lft=5,rgt=2
    real(8) :: tempR(lft,lft),tempI(lft,lft),vecR(lft,rgt),vecI(lft,rgt)
    complex(8) :: thedata(lft,lft),vectordata(lft,rgt)
    type(tensor2) :: matrix,vector,solution,originalMatrix
    integer :: errorpoint

    call random_number(tempR)
    call random_number(tempI)
    thedata=tempR+II*tempI
    call random_number(vecR)
    call random_number(vecI)
    vectordata=vecR+II*vecI

    vector=new_Tensor(vectorData)
    matrix=new_Tensor(thedata)
    originalMatrix=new_Tensor(matrix)
    solution=SolveLinearProblem(matrix, vector, 1.0d-8 )
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((0.0d0 &
     +1.0d-8) &
     .ge. &
     ( matrix * solution.absdiff.vector) &
             .and. &
     (0.0d0 &
     -1.0d-8) &
     .le. &
     ( matrix * solution.absdiff.vector) )) then
      print *, " *Assert_Equal_Within failed* in test LinearSolver &
              &[Tensor_Class.fun:812]"
      print *, "  ", " matrix * solution.absdiff.vector (", matrix * solution.absdiff.vector,") is not", &
 0.0d0,"within",1.0d-8
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

    !call ZGELSD( thedata, vectordata, 1.0d-8 )


  numTests = numTests + 1

 end subroutine LinearSolver


 subroutine funit_setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
  CALL random_seed()
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
  call tensor3_joinindices_first
  call funit_teardown

  call funit_setup
  call tensor4_joinindices
  call funit_teardown

  call funit_setup
  call tensor2_Splitindex
  call funit_teardown

  call funit_setup
  call transpose_of_tensor
  call funit_teardown

  call funit_setup
  call transpose_of_tensor4
  call funit_teardown

  call funit_setup
  call unfoldings_of_tensor3
  call funit_teardown

  call funit_setup
  call unfoldings_of_tensor4
  call funit_teardown

  call funit_setup
  call matrix_product
  call funit_teardown

  call funit_setup
  call Singular_Value_Decomposition
  call funit_teardown

  call funit_setup
  call Right_Compactification
  call funit_teardown

  call funit_setup
  call Left_Compactification
  call funit_teardown

  call funit_setup
  call Compact_From_Below_T3_T4
  call funit_teardown


  call funit_setup
  call Matrices_times_tensor3s
  call funit_teardown

  call funit_setup
  call xplusWithTensor4
  call funit_teardown

  call funit_setup
  call TensorTraceWithTensor4
  call funit_teardown

  call funit_setup
  call LinearSolver
  call funit_teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_Tensor_Class

end module Tensor_Class_fun
