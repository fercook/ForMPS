module Constants
  
  complex(8), parameter :: II=(0.0d0,1.0d0)
  complex(8), parameter :: one=(1.0d0,0.0d0)
  complex(8), parameter :: zero=(0.0d0,0.0d0)
  real(8),parameter :: Pi=3.141592653589793

  integer,parameter :: MatrixSpin =1
  integer,parameter :: FirstDimension=1
  integer,parameter :: SecondDimension=2
  integer,parameter :: ThirdDimension=3
  integer,parameter :: FourthDimension=4

  integer,parameter :: First  = 1
  integer,parameter :: Second = 2
  integer,parameter :: Third = 3
  integer,parameter :: FirstAndSecond = 12
  integer,parameter :: FirstAndThird = 13
  integer,parameter :: SecondAndFirst = 21
  integer,parameter :: SecondAndThird = 23
  integer,parameter :: ThirdAndSecond = 32
  integer,parameter :: ThirdAndFirst = 31


  integer,parameter :: Right = 1
  integer,parameter :: Left =2
  
end module Constants
