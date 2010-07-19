test_suite Tensor_Class

!TODO: New tests with all possible combinations of index bonding
!TODO: Some other tests about creation, products, and so on

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

  aMatrix=aMPS%JoinIndices(FIRSTANDSECOND,THIRD)
  assert_equal_within(amatrix.absdiff.correct, 0.0d0, 1.0e-8)

  correct=new_Tensor(transpose(matrix)) !This will be 4,2*3
  aMatrix=JoinIndicesOf(aMPS,THIRD,FIRSTANDSECOND)
  assert_equal_within(amatrix.absdiff.correct, 0.0d0, 1.0e-8)

   assert_equal(aMPS%delete(),Normal)
   assert_equal(aMatrix%delete(),Normal)
   assert_equal(correct%delete(),Normal)
   assert_false(WasThereError())
end test

test tensor4_joinIndices
  type(tensor4) :: aTensor
  type(tensor2) :: aMatrix,correct
  complex(8) :: data(2,3,4,5),matrix(6,20)
  integer error,i,j,k,l

  !Initialization
  forall (i=1:2 ,j=1:3, k=1:4, l=1:5) data(i,j,k,l)=ONE*(i+(j-1)*2+(k-1)*3*2+(l-1)*2*3*4)
  aTensor=new_Tensor(data)

  matrix=one*Reshape( data, [2*3,4*5])
  correct=new_Tensor(matrix)

  aMatrix=aTensor%JoinIndices(FIRSTANDSECOND,THIRDANDFOURTH)
  assert_equal_within(amatrix.absdiff.correct, 0.0d0, 1.0e-8)

  error=correct%Delete()
  correct=new_Tensor(transpose(matrix))

  aMatrix=aTensor%JoinIndices(THIRDANDFOURTH,FIRSTANDSECOND)
  print *,'returning'
  assert_equal_within(amatrix.absdiff.correct, 0.0d0, 1.0e-8)

end test



test tensor2_SplitIndex
  type(tensor3) :: aMPS,correct
  type(tensor2) :: aMatrix
  complex(8) :: matrixdata(6,4)
  integer error,i,j,k

  !Initialization
  forall (i=1:6 ,j=1:4) matrixdata(i,j)=ONE*(i+(j-1)*6)
  aMatrix=new_Tensor(matrixdata)

  correct=new_Tensor(one*Reshape( matrixdata, [2,3,4]))

  aMPS=SplitIndexOf(aMatrix,FIRST,2)
  assert_equal_within(aMPS.absdiff.correct, 0.0d0, 1.0e-8)

  correct=new_Tensor(one*Reshape( matrixdata, [6,2,2]))

  aMPS=SplitIndexOf(aMatrix,SECOND,2)
  assert_equal_within(aMPS.absdiff.correct, 0.0d0, 1.0e-8)

  assert_false(WasThereError())
end test

test transpose_of_tensor
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
    assert_equal_within(bTensor.absdiff.correctTensor, 0.0d0, 1.0e-8)
    deallocate(correct)

    neworder=[2,1,3]
    bTensor=ConjugateTranspose(aTensor,neworder)
    allocate(correct( d2,d1,d3 ))
    forall (i=1:d1, j=1:d2, k=1:d3) correct(j,i,k)=one*(i*100+j*10+II*k)
    CorrectTensor=new_Tensor(dconjg(correct))
    assert_equal_within(bTensor.absdiff.correctTensor, 0.0d0, 1.0e-8)
    deallocate(correct)

    neworder=[1,3,2]
    bTensor=ConjugateTranspose(aTensor,neworder)
    allocate(correct( d1,d3,d2 ))
    forall (i=1:d1, j=1:d2, k=1:d3) correct(i,k,j)=one*(i*100+j*10+II*k)
    CorrectTensor=new_Tensor(dconjg(correct))
    assert_equal_within(bTensor.absdiff.correctTensor, 0.0d0, 1.0e-8)
    deallocate(correct)

end test

test matrix_product
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

    assert_equal_within(mat .absdiff. correct, 0.0d0, 1.0e-10)
end test


test Singular_Value_Decomposition
   integer,parameter :: DLeft=6,DRight=8
   type(Tensor2) :: aMatrix,theU,theVt
   type(Tensor2) :: theSigma
   complex(8) :: data(DLeft,DRight)
   integer :: i,j,k

   forall (i=1:Dleft ,j=1:Dright) &
      & data(i,j)=one*(exp(II*i*Pi/DLeft)+(j-1)*Dleft)
    aMatrix=new_Tensor(data)

   call aMatrix%SVD(theU,theSigma,theVt)

   assert_equal_within(aMatrix.absdiff.(theU*(theSigma*theVt)), 0.0d0, 1.0e-10)
   assert_false(WasThereError())

end test



end test_suite
