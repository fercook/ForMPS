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
!    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

module Constants
  
  complex(8), parameter :: II=(0.0d0,1.0d0)
  complex(8), parameter :: one=(1.0d0,0.0d0)
  complex(8), parameter :: zero=(0.0d0,0.0d0)
  real(8),parameter :: Pi=3.141592653589793

  integer,parameter :: FirstDimension=1
  integer,parameter :: SecondDimension=2
  integer,parameter :: ThirdDimension=3
  integer,parameter :: FourthDimension=4

  integer,parameter,dimension(1) :: First  = [1]
  integer,parameter,dimension(1) :: Second = [2]
  integer,parameter,dimension(1) :: Third = [3]
  integer,parameter,dimension(1) :: Fourth = [4]
  integer,parameter,dimension(1) :: Fifth = [5]

  integer,parameter,dimension(2) :: FirstAndSecond = [1,2]
  integer,parameter,dimension(2) :: FirstAndThird = [1,3]
  integer,parameter,dimension(2) :: SecondAndFirst = [2,1]
  integer,parameter,dimension(2) :: SecondAndThird = [2,3]
  integer,parameter,dimension(2) :: ThirdAndSecond = [3,2]
  integer,parameter,dimension(2) :: ThirdAndFirst = [3,1]

  integer,parameter,dimension(2) :: FirstAndFourth = [1,4]
  integer,parameter,dimension(2) :: SecondAndFourth = [2,4]
  integer,parameter,dimension(2) :: ThirdAndFourth = [3,4]
  integer,parameter,dimension(2) :: FourthAndFirst = [4,1]
  integer,parameter,dimension(2) :: FourthAndSecond = [4,2]
  integer,parameter,dimension(2) :: FourthAndThird = [4,3]

  integer,parameter :: Right = 1
  integer,parameter :: Left = -1
  integer,parameter :: No = 0
  
  interface operator (.equalvector.)
     module procedure AreIntVectorsEqual
  end interface

  contains

  logical function AreIntVectorsEqual(vector1,vector2) result(AreEqual)
     integer,intent(IN) :: vector1(:),vector2(:)
     integer n

     if(size(vector1).eq.size(vector2)) then
        AreEqual=.true.
        do n=1,size(vector1)
           AreEqual=AreEqual.and.(vector1(n).eq.vector2(n))
        enddo
     else
        AreEqual=.false.
     endif
   end function

end module Constants
