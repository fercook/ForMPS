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
  integer, parameter :: integerOne=1
  integer, parameter :: integerZero=0

  integer,parameter :: SpinDimension=0
  integer,parameter :: FirstDimension=1
  integer,parameter :: SecondDimension=2
  integer,parameter :: ThirdDimension=3
  integer,parameter :: FourthDimension=4
  integer,parameter :: FifthDimension=5
  integer,parameter :: SixthDimension=6

  integer,parameter,dimension(1) :: First  = [1]
  integer,parameter,dimension(1) :: Second = [2]
  integer,parameter,dimension(1) :: Third = [3]
  integer,parameter,dimension(1) :: Fourth = [4]
  integer,parameter,dimension(1) :: Fifth = [5]
  integer,parameter,dimension(1) :: Sixth = [6]

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

  enum,bind(C)
  enumerator :: NOWHERE,LEFT, RIGHT, UP, DOWN
  end enum

  enum,bind(C)
  enumerator :: HORIZONTAL,VERTICAL
  end enum

  integer,parameter :: No = 0
  integer,parameter :: Yes = 1
  integer,parameter :: DONOTCONJUGATE=NO
  integer,parameter :: UNDEFINED = -999
  integer,parameter :: ALLTENSORS = 1313789
  integer,parameter :: NOLIMIT = 9999999

  logical,parameter :: Verbose = .true.

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

   integer function DirectionOppositeTo(aDirection) result(theOpposite)
      integer, intent(IN) :: aDirection

      select case (aDirection)
         case (UP)
            theOpposite=DOWN
         case (DOWN)
            theOpposite=UP
         case (LEFT)
            theOpposite=RIGHT
         case (RIGHT)
            theOpposite=LEFT
      end select
   end function DirectionOppositeTo


end module Constants
