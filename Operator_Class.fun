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

test_suite Operator_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown

test ApplyOperator
   type(MPSTensor) :: A,B,C
   integer,parameter :: DleftT=4, DrightT=3,spinT=2
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: Correct(DleftT,DrightT,spinT)
   integer error

   !Initialization
   data(:,:,1)=one*Reshape( [-0.120802, -0.736476, 0.639716, -0.0536866, -0.304505, -0.320898, &
   			   & -0.325939, -0.200223, -0.488209, 0.0946794, 0.17894, 0.76573], [DleftT,DrightT])
   data(:,:,2)=one*Reshape( [-0.258579, -0.424792, -0.512473, -0.130085, -0.442283, -0.00921493, &
   			   & -0.298086, 0.183901, -0.625986, 0.406363, 0.317842, -0.565636], [DleftT,DrightT])
   Correct(:,:,1)=data(:,:,1)
   Correct(:,:,2)=(-1.0d0)*data(:,:,2)

     A=new_MPSTensor(data(:,:,1),data(:,:,2))
     B=new_MPSTensor(Correct(:,:,1),Correct(:,:,2))
     C=A.Apply.PauliSigma(3)
     assert_equal_within(B.absdiff.C, 0.0d0, 1.0e-5)
     assert_false(WasThereError())
end test

end test_suite

