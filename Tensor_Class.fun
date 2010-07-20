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

test Right_Compactification

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
        matrix(i,j)=(II**i+(j-1)*DrightT)
     enddo
  enddo

  CorrectResult=one*reshape([3072 - 312*II, 3528 - 312*II, 3984 - 312*II, 4440 - 312*II, 3384 - &
       & 360*II, 3888 - 360*II, 4392 - 360*II, 4896 - 360*II, 3696 - &
       & 408*II, 4248 - 408*II, 4800 - 408*II, 5352 - 408*II, 4008 - &
       & 456*II, 4608 - 456*II, 5208 - 456*II, 5808 - 456*II], [DleftT,DleftT] )
  T_Tensor=new_Tensor(data)
  T_Matrix=new_Tensor(matrix)
  T_Correct=new_Tensor(CorrectResult)
  T_Compacted=CompactRight(T_Matrix,T_Tensor,T_Tensor,THIRD)
  assert_equal_within(T_Compacted.absdiff.T_Correct, 0.0d0, 1.0e-8)

end test

test Left_Compactification
!
!  Mathematica code: NOTICE THE TRANSPOSE TO GET THE ORDER RIGHT
! With[{DL = 3, DR = 4, spin = 2},
!   At = Table[
!     DR*(s - 1) + (j - 1)*DL + i, {s, 1, 2}, {i, 1, DL}, {j, 1, DR}];
!   mat = Table[I^i + (j - 1)*DL, {i, 1, DL}, {j, 1, DL}]];
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
        matrix(i,j)=(II**i+(j-1)*DleftT)
     enddo
  enddo

  CorrectResult=one*reshape([1104 - 48*II, 1788 - 48*II, 2472 - 48*II, 3156 - 48*II, 1680 - &
       &  84*II, 2796 - 84*II, 3912 - 84*II, 5028 - 84*II, 2256 - 120*II, 3804 - &
       &  120*II, 5352 - 120*II, 6900 - 120*II, 2832 - 156*II, 4812 - &
       &  156*II, 6792 - 156*II, 8772 - 156*II], [DrightT,DrightT] )

  T_Tensor=new_Tensor(data)
  T_Matrix=new_Tensor(matrix)
  T_Correct=new_Tensor(CorrectResult)
  T_Compacted=CompactLeft(T_Matrix,T_Tensor,T_Tensor,THIRD)
  assert_equal_within(T_Compacted.absdiff.T_Correct, 0.0d0, 1.0e-8)

end test

end test_suite
