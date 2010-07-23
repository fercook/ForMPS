!!   Copyright 2010 Fernando M. Cucchietti
!
!    This file is part of FortranMPS
!
!    FortranMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    FortranMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with FortranMPS.  If not, see <http://www.gnu.org/licenses/>.


test_suite MPSTensor_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown


test type_creation_deletion
  type(MPSTensor) :: A,B
  integer :: spinT=2,BondL=3,BondR=4

  A=new_MPSTensor(SpinT,BondL,BondR)
  assert_false(WasThereError())
  assert_equal(A%delete(),Normal)
  call LowerFlag()
  B=new_MPSTensor(SpinT,BondL,BondR)
  call B%Print()
  assert_equal(B%delete(),Normal)
  assert_false(WasThereError())
end test

!

test accesors_work
  type(MPSTensor) :: A
  integer :: spinT=2,BondL=3,BondR=4

  A=new_MPSTensor(SpinT,BondL,BondR)
  assert_equal(A%getspin(),SpinT)
  assert_equal(A%getDLeft(),BondL)
  assert_equal(A%getDRight(),BondR)
  assert_equal(A%delete(),Normal)
  assert_false(WasThereError())
end test



test LeftCanonTensor
    !For some reason the output matrix U from SVD is different than the one from Mathematica
    ! but only in the non-important vectors (where Sigma=0)
    !Then, this test is changed to an automatic test that sees if the original tensor can be
    !reconstructed by multiplying the canonized tensor with the output matrix
   type(MPSTensor) :: TheTensor,CorrectTensor,origTensor
   type(Tensor2) :: outputMatrix,correctMatrix
   integer,parameter :: spinT=2,DleftT=4, DrightT=3
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: matrix(DrightT,DrightT)
   complex(8) :: Correct(DleftT,DrightT,spinT)
   integer :: i,j,k

   !Initialization
   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT)
   Correct(:,:,1)=one*Reshape( [-0.234451, -0.271825, -0.309198, -0.346571, -0.625669, -0.440419, &
    & -0.255169, -0.0699187, -0.295271, 0.429314, 0.0856099, -0.126573], [DleftT,DrightT])
   Correct(:,:,2)=Reshape( [-0.346571, -0.383944, -0.421317, -0.45869, -0.0699187, 0.115331, &
    & 0.300581, 0.485831, -0.126573, 0.354249, -0.66644, 0.345686], [DleftT,DrightT])
   matrix(:,:)=Reshape( [-12.1367, 2.9496, 0., -23.227, 0.712199, 0., &
           & -34.3173, -1.5252, 0.], [DrightT,DrightT] )
   TheTensor=new_MPSTensor(data(:,:,1),data(:,:,2))
   origTensor=TheTensor
   !CorrectTensor=new_MPSTensor(Correct(:,:,1),Correct(:,:,2))
   outputMatrix=LeftCanonize(TheTensor) !  C=A%LCanonize()
   correctMatrix=new_Tensor(matrix)
   call outputMatrix%PrintDimensions()
   call correctMatrix%PrintDimensions()
   !call theTensor%Print()
   !call correctTensor%Print()
   CorrectTensor=MPSTensor_times_matrix(TheTensor,outputMatrix)
   assert_equal_within(outputMatrix.absdiff.correctMatrix, 0.0d0, 1.0e-4)
   assert_equal_within(OrigTensor.absdiff.CorrectTensor, 0.0d0, 1.0e-5)

end test


test RightCanonTensor
    !For some reason the output matrix U from SVD is different than the one from Mathematica
    ! but only in the non-important vectors (where Sigma=0)
    !Then, this test is changed to an automatic test that sees if the original tensor can be
    !reconstructed by multiplying the canonized tensor with the output matrix
   type(MPSTensor) :: TheTensor,CorrectTensor,origTensor
   type(Tensor2) :: outputMatrix
   integer,parameter :: spinT=2,DleftT=7, DrightT=5
   complex(8) :: data(DleftT,DrightT,spinT)
   integer :: i,j,k

   !Initialization
   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i*II+(j-1)*DleftT+(k-1)*DrightT)
   TheTensor=new_MPSTensor(data(:,:,1),data(:,:,2))
   origTensor=TheTensor
   outputMatrix=RightCanonize(TheTensor) !  C=A%LCanonize()

   CorrectTensor=matrix_times_MPSTensor(outputMatrix,TheTensor)

   assert_equal_within(OrigTensor.absdiff.CorrectTensor, 0.0d0, 1.0e-5)

end test
!
!
!
!
!test RightCanonTensor
!   type(MPSTensor) :: A,B,C,D
!   integer,parameter :: DleftT=4, DrightT=3
!   complex(8) :: data(DleftT,DrightT,spinT)
!   complex(8) :: matrix(DleftT,DleftT,1)
!   complex(8) :: Correct(DleftT,DrightT,spinT)
!   integer :: i,j,k
!
!   !Initialization
!   forall (i=1:DleftT ,j=1:DrightT, k=1:SpinT) data(i,j,k)=one*(i+(j-1)*DleftT+(k-1)*DrightT)
!   Correct(:,:,1)=one*Reshape( [-0.120802, -0.736476, 0.639716, -0.0536866, -0.304505, -0.320898, &
!   			   & -0.325939, -0.200223, -0.488209, 0.0946794, 0.17894, 0.76573], [DleftT,DrightT])
!   Correct(:,:,2)=one*Reshape( [-0.258579, -0.424792, -0.512473, -0.130085, -0.442283, -0.00921493, &
!   			   & -0.298086, 0.183901, -0.625986, 0.406363, 0.317842, -0.565636], [DleftT,DrightT])
!   matrix(:,:,1)=one*Reshape ( [-18.1216, -20.362, -22.6023, -24.8427, 1.61461, 0.624272, -0.366067, &
!   			 & -1.35641, 0., 0., 0., 0., 0., 0., 0., 0.], [DleftT,DleftT] )
!     A=new_MPSTensor(SpinT,DleftT,DrightT,data)
!     B=new_MPSTensor(SpinT,DleftT,DrightT,Correct)
!     C=Right_Canonize_MPSTensor(A)!C=A%RCanonize()
!     D=new_MPSTensor(MatrixSpin,DleftT,DleftT,matrix)
!     assert_equal_within(A.absdiff.B, 0.0d0, 1.0e-5)
!     assert_equal_within(C.absdiff.D, 0.0d0, 1.0e-4)
!     assert_equal(A%delete(),Normal)
!     assert_equal(B%delete(),Normal)
!     assert_equal(C%delete(),Normal)
!     assert_equal(D%delete(),Normal)
!     assert_false(WasThereError())
!end test
!
!

end test_suite



!  print *,shape(data(:,:,:))
! print *,'And now'
!  print *,shape(transpose(data(:,:,:)))
!!  print *,'And finally'
!  print *,shape(reshape(data,[DLeftT,spinT*DrightT]))
!  assert_true(.true.)

