!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010

!!TODO: Reimplement extendind the type from a Matrix class that performs some low level functions

module MPSProductStorer_Class

  use ErrorHandling
  use Constants
  use MPSTensor_Class
  use Operator_Class
  use MatrixProductState_Class

  implicit none

!###############################
!#####  The class main object
!###############################
  type MPSProductStorer
     private
     complex(8),allocatable,dimension(:) :: storage(:,:) 
     integer :: Length
     !Next is Right or Left, defined in Constants
     !It is ugly because gfortran has a bug with polymorphism
     !Will be changed in the future to a factory that provides right or left Storers
     integer :: Sense 
!   contains
!     procedure :: ApplyToTensor => ApplyToMPSTensor
  end type MPSProductStorer

!################# Methods

 contains

 end module MPSProductStorer_Class
