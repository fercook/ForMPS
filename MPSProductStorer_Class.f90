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
