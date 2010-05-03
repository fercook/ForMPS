test_suite MPSTensor_Class
  type(MPSTensor) :: A,B,C,D,E,F
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
  assert_false(WasThereError())
  assert_false(B%Print().eq.Normal)
  call LowerFlag()
  B=new_MPSTensor(SpinT,BondL,BondR)
  assert_equal(B%Print(),Normal)
  assert_equal(B%delete(),Normal)
  assert_false(WasThereError())
end test

test accesors_work
  A=new_MPSTensor(SpinT,BondL,BondR)
  assert_equal(A%spin(),SpinT)
  assert_equal(A%DLeft(),BondL)
  assert_equal(A%DRight(),BondR)
  assert_false(WasThereError())
end test

test Comparison_Of_Dimensions
  A=new_MPSTensor(SpinT,BondL,BondR)
  B=new_MPSTensor(2,4,4)
  assert_false(A.equaldims.B)
  error=B%delete()
  B=new_MPSTensor(SpinT,BondL,BondR)
  assert_true(A.equaldims.B)
  assert_false(WasThereError())
end test

test Multiplication_By_Number
  integer :: I = 2
  real :: R = 2.0
  real(8) :: R8 = 2.0d0
  complex :: C = 2.0
  complex(8) :: C8 = 2.0d0
  A=new_MPSTensor(SpinT,BondL,BondR)
  B=new_MPSTensor(A)
  assert_equal_within((I*A).diff.(R*B),0.0d0, 1.0e-8)
  assert_equal_within((R8*A).diff.(C*B),0.0d0, 1.0e-8)
  assert_equal_within((I*A).diff.(C8*B),0.0d0, 1.0e-8)
  assert_false(WasThereError())
end test

test L_product_single
  integer,parameter :: DleftT=1, DrightT=4
  complex(8) :: data(DleftT,DrightT,spinT)
  complex(8) :: CorrectResult(DrightT,DrightT,1)
  integer :: shape(3)
  integer n,i,j,k,s
  do i=1,DleftT
     do j=1,DrightT
        do s=1,spinT
           data(i,j,s)=one*(i+(j-1)*DleftT+(s-1)*DrightT)
        enddo
     enddo
  enddo
  shape = [DrightT,DrightT,1]
  CorrectResult=reshape([26, 32, 38, 44, 32, 40, 48, 56, 38, 48, 58, 68, 44, 56, 68, 80], shape)
  A=new_MPSTensor(SpinT,DLeftT,DRightT,data)
  C=new_MPSTensor(1,DrightT,DrightT,CorrectResult)
  D=MPS_Left_Product(A,A)
  assert_equal_within(D.diff.C, 0.0d0, 1.0e-8)
  assert_false(WasThereError())
end test


test L_product_withMatrix
!
!  Mathematica code: NOTICE THE TRANSPOSE TO GET THE ORDER RIGHT
! With[{DL = 3, DR = 4, spin = 2}, 
!   At = Table[
!     DR*(s - 1) + (j - 1)*DL + i, {s, 1, 2}, {i, 1, DL}, {j, 1, DR}];
!   mat = Table[I^i + (j - 1)*DL, {i, 1, DL}, {j, 1, DL}]];
!  Flatten[Transpose[LProduct[At, At, mat]]]
!
  integer,parameter :: DleftT=3, DrightT=4
  complex(8) :: data(DleftT,DrightT,spinT)
  complex(8) :: matrix(DleftT,DleftT,1)
  complex(8) :: CorrectResult(DrightT,DrightT,1)
  integer :: shape(3)
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
        matrix(i,j,1)=(II**i+(j-1)*DleftT)
     enddo
  enddo
  shape = [DrightT,DrightT,1]
  CorrectResult=reshape([1104 - 48*II, 1788 - 48*II, 2472 - 48*II, 3156 - 48*II, 1680 - &
       &  84*II, 2796 - 84*II, 3912 - 84*II, 5028 - 84*II, 2256 - 120*II, 3804 - &
       &  120*II, 5352 - 120*II, 6900 - 120*II, 2832 - 156*II, 4812 - &
       &  156*II, 6792 - 156*II, 8772 - 156*II], shape)
  A=new_MPSTensor(SpinT,DLeftT,DRightT,data)
  B=new_MPSTensor(1,DleftT,DleftT,matrix)
  C=new_MPSTensor(1,DrightT,DrightT,CorrectResult)
  D=MPS_Left_Product(A,A,B)
  assert_equal_within(D.diff.C, 0.0d0, 1.0e-8)
  assert_false(WasThereError())
end test

test R_product_single
!!  Mathematica test code:
!!  With[{DL = 4, DR = 1, spin = 2}, 
!!  At = Table[
!!    DR*(s - 1) + (j - 1)*DL + i, {s, 1, 2}, {i, 1, DL}, {j, 1, DR}]];
!!  Flatten[Transpose[RProduct[At, At]]]
!!
  integer,parameter :: DleftT=4, DrightT=1
  complex(8) :: data(DleftT,DrightT,spinT)
  complex(8) :: CorrectResult(DleftT,DleftT,1)
  integer :: shape(3)
  integer n,i,j,k,s
  do i=1,DleftT
     do j=1,DrightT
        do s=1,spinT
           data(i,j,s)=one*(i+(j-1)*DleftT+(s-1)*DrightT)
        enddo
     enddo
  enddo
  shape = [DleftT,DleftT,1]
  CorrectResult=reshape([5, 8, 11, 14, 8, 13, 18, 23, 11, 18, 25, 32, 14, 23, 32, 41], shape)
  A=new_MPSTensor(SpinT,DLeftT,DRightT,data)
  C=new_MPSTensor(1,DLeftT,DLeftT,CorrectResult)
  D=MPS_Right_Product(A,A)
  assert_equal_within(D.diff.C, 0.0d0, 1.0e-8)
  assert_false(WasThereError())
end test


test R_product_withMatrix
!! Mathematica code: (NOTICE TRANSPOSE TO GET ORDER RIGHT)
!!
!! With[{DL = 4, DR = 3, spin = 2}, 
!!  At = Table[
!!    DR*(s - 1) + (j - 1)*DL + i, {s, 1, 2}, {i, 1, DL}, {j, 1, DR}];
!!  mat = Table[I^i + (j - 1)*DR, {i, 1, DR}, {j, 1, DR}]];
!! Flatten[Transpose[RProduct[At, At, mat]]]
!!
  integer,parameter :: DleftT=4, DrightT=3
  complex(8) :: data(DleftT,DrightT,spinT)
  complex(8) :: matrix(DrightT,DrightT,1)
  complex(8) :: CorrectResult(DleftT,DleftT,1)
  integer :: shape(3)
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
        matrix(i,j,1)=(II**i+(j-1)*DrightT)
     enddo
  enddo
  shape = [DleftT,DleftT,1]
  CorrectResult=reshape([3072 - 312*II, 3528 - 312*II, 3984 - 312*II, 4440 - 312*II, 3384 - &
       & 360*II, 3888 - 360*II, 4392 - 360*II, 4896 - 360*II, 3696 - & 
       & 408*II, 4248 - 408*II, 4800 - 408*II, 5352 - 408*II, 4008 - &
       & 456*II, 4608 - 456*II, 5208 - 456*II, 5808 - 456*II], shape)
  A=new_MPSTensor(SpinT,DLeftT,DRightT,data)
  B=new_MPSTensor(1,DrightT,DrightT,matrix)
  C=new_MPSTensor(1,DleftT,DleftT,CorrectResult)
  D=MPS_Right_Product(A,A,B)
  assert_equal_within(D.diff.C, 0.0d0, 1.0e-8)
  assert_false(WasThereError())
end test

test Mult_by_Matrix
  integer,parameter :: DleftT=4, DrightT=3
  complex(8) :: data(DleftT,DrightT,spinT)
  complex(8) :: matrix(DrightT,DLeftT,1)
  complex(8) :: CorrectByRight(DleftT,DLeftT,spinT)
  complex(8) :: CorrectByLeft(DrightT,DrightT,spinT)

  integer :: i,j,k
  !Initialization
  forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT)
  forall (j=1:DleftT ,i=1:DrightT ) matrix(i,j,1)=(II**i+(j-1)*DrightT)

  CorrectByRight(:,:,1)=Reshape( [-5 - 8*II, -6 - 8*II, -7 - 8*II, -8 - 8*II, 40 - 8*II, 48 - 8*II, 56 - &
      &  8*II, 64 - 8*II, 85 - 8*II, 102 - 8*II, 119 - 8*II, 136 - 8*II, 130 - &
      &  8*II, 156 - 8*II, 182 - 8*II, 208 - 8*II], [DleftT,DleftT])
  CorrectByRight(:,:,2)=Reshape( [-8 - 8*II, -9 - 8*II, -10 - 8*II, -11 - 8*II, 64 - 8*II, 72 - 8*II, 80 - &
      & 8*II, 88 - 8*II, 136 - 8*II, 153 - 8*II, 170 - 8*II, 187 - 8*II, 208 - &
      & 8*II, 234 - 8*II, 260 - 8*II, 286 - 8*II],[DLeftT,DLeftT])
  CorrectByLeft(:,:,1)=Reshape( [60 + 10*II, 50, 60 - 10*II, 132 + 26*II, 106, 132 - 26*II, 204 + &
      & 42*II, 162, 204 - 42*II],[DRightT,DRightT])
  CorrectByLeft(:,:,2)=Reshape([114 + 22*II, 92, 114 - 22*II, 186 + 38*II, 148, 186 - 38*II, 258 + &
      & 54*II, 204, 258 - 54*II],[DRightT,DRightT])

   A=new_MPSTensor(SpinT,DLeftT,DrightT,data)
   B=new_MPSTensor(1,DrightT,DLeftT,matrix)

   C=new_MPSTensor(spinT,DrightT,DrightT,CorrectByLeft)
   D=new_MPSTensor(2,DrightT,DrightT)
   D=B*A !Matrix_times_MPSTensor(B,A)
   assert_equal_within(D.diff.C, 0.0d0, 1.0e-8)

   i=C%delete()
   i=D%delete()
   C=new_MPSTensor(spinT,DLeftT,DLeftT,CorrectByRight)
   D=new_MPSTensor(2,DLeftT,DLeftT)
   D=A*B !Matrix_times_MPSTensor(A,B)
   assert_equal_within(D.diff.C, 0.0d0, 1.0e-8)
   assert_false(WasThereError())	 				     
end test


test Collapse1
   integer,parameter :: DleftT=4, DrightT=3
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: matrix(SpinT*DleftT,DrightT)
   complex(8) :: Correct(SpinT*DleftT,DrightT)
   integer :: i,j,k

   !Initialization
   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT) 
   Correct=Reshape( [1., 2., 3., 4., 4., 5., 6., 7., 5., 6., 7., 8., 8., 9., 10., 11., 9., &
                   & 10., 11., 12., 12., 13., 14., 15.], [SpinT*DleftT,DrightT])
     A=new_MPSTensor(SpinT,DleftT,DrightT,data)
     call CollapseSpinWithBond(A,matrix,FirstDimension)
     assert_equal_within(Difference_btw_Matrices(matrix,correct), 0.0d0, 1.0e-8)
     B=new_MPSTensor(SpinT,DleftT,DrightT,zero)
     call SplitSpinFromBond(matrix,B,FirstDimension,DleftT,DrightT)
     assert_equal_within(A.diff.B, 0.0d0, 1.0e-8)
     assert_false(WasThereError())

end test
     
test Collapse2
   integer,parameter :: DleftT=4, DrightT=3
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: matrix(DleftT,SpinT*DrightT)
   complex(8) :: Correct(DleftT,SpinT*DrightT)
   integer :: i,j,k

   !Initialization
   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT) 
   Correct=Reshape( [1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 4., 5., 6., 7., &
   		    & 8., 9., 10., 11., 12., 13., 14., 15.], [DleftT,SpinT*DrightT])
     A=new_MPSTensor(SpinT,DleftT,DrightT,data)
     call CollapseSpinWithBond(A,matrix,SecondDimension)
     assert_equal_within(Difference_btw_Matrices(matrix,correct), 0.0d0, 1.0e-8)
     B=new_MPSTensor(SpinT,DleftT,DrightT,zero)
     call SplitSpinFromBond(matrix,B,SecondDimension,DleftT,DrightT)
     assert_equal_within(A.diff.B, 0.0d0, 1.0e-8)
     assert_false(WasThereError())
end test

test SingularVD
   integer,parameter :: DleftT=4, DrightT=3
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: matrix(SpinT*DleftT,DrightT)
   complex(8) :: U(SpinT*DleftT,SpinT*DleftT)
   complex(8) :: V(DrightT,DrightT)
   real(8) :: sigma(DrightT)
   complex(8) :: diagonal(SpinT*DleftT,DrightT)
   complex(8) :: Correct(SpinT*DleftT,DrightT)
   integer :: i,j,k

   !Initialization
   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT) 
   A=new_MPSTensor(SpinT,DleftT,DrightT,data)
   call CollapseSpinWithBond(A,matrix,FirstDimension)
   Correct=matrix
   i=SingularValueDecomposition(matrix,U,Sigma,V)
   assert_equal(i,0)
   diagonal=zero
   do i=1,DRightT
      diagonal(i,i)=Sigma(i)
   enddo
   assert_equal_within(Difference_btw_Matrices(Correct,matmul(Matmul(U,diagonal),V)), 0.0d0, 1.0e-10)
   assert_false(WasThereError())
end test

test LeftCanon
   integer,parameter :: DleftT=4, DrightT=3
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: matrix(DleftT,DLeftT,1)
   complex(8) :: Correct(DleftT,DrightT,spinT)
   integer :: i,j,k

   !Initialization
   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT) 
   Correct(:,:,1)=Reshape( [-0.120802, -0.736476, 0.639716, -0.0536866, -0.304505, -0.320898, &
-0.325939, -0.200223, -0.488209, 0.0946794, 0.17894, 0.76573], [DleftT,DrightT])
   Correct(:,:,1)=Reshape( [-0.258579, -0.424792, -0.512473, -0.130085, -0.442283, -0.00921493, &
-0.298086, 0.183901, -0.625986, 0.406363, 0.317842, -0.565636], [DleftT,DrightT])
   matrix(:,:,1)=Reshape ( [-18.1216, -20.362, -22.6023, -24.8427, 1.61461, 0.624272, -0.366067, &
   			 & -1.35641, 0., 0., 0., 0., 0., 0., 0., 0.], [DLeftT,DLeftT] )
     A=new_MPSTensor(SpinT,DleftT,DrightT,data)
     B=new_MPSTensor(SpinT,DleftT,DrightT,Correct)
     C=A%LCanonize()
     D=new_MPSTensor(SpinT,DleftT,DleftT,matrix)
     assert_equal_within(A.diff.B, 0.0d0, 1.0e-8)
     assert_equal_within(C.diff.D, 0.0d0, 1.0e-8)
     assert_false(WasThereError())

end test
     


end test_suite



!  print *,shape(data(:,:,:))
! print *,'And now'
!  print *,shape(transpose(data(:,:,:)))
!!  print *,'And finally'
!  print *,shape(reshape(data,[DLeftT,spinT*DrightT]))
!  assert_true(.true.)
