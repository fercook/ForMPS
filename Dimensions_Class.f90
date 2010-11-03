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



! THIS CLASS IS NOT USED -- IT MIGHT BE USED WHEN CODE IS FURTHER
! REFACTORED, OR NOT.


module Dimensions_Class

  use ErrorHandling
  use Constants

  implicit none

  private

  type :: Dimensions
  	private
  	integer,allocatable :: dimensionList(:)
  	integer :: numberOfDimensions=0
	integer :: IndexOfSpecialDimension=0
  contains
  	procedure :: IsInitialized => Is_Initialized
    procedure :: delete => delete_Dimensions
  	procedure :: print => print_Dimensions
    procedure :: getOrderedDimensions => get_OrderedDimensions
    procedure :: getDimensions => get_UnorderedDimensions
    procedure :: length => getNumberOfDimensions
	procedure ::
  end type Dimensions

	interface new_Dimensions
		new_Dimensions_fromList
	end interface
contains

	function new_Dimensions_fromList( dimList, specialIndex)
		type(Dimensions) :: this
		integer :: specialIndex
		integer :: dimList(:)
		integer :: Length

		Length=size(dimList)
		allocate(this%dimensionList(Length))
		this%dimensionList=dimList
		this%numberOfDimensions=Length
		this%IndexOfSpecialDimension=specialIndex
		return
	end function new_Dimensions_fromList


	integer function delete_Dimensions(this)
		type(Dimensions) :: this
		if(this%IsInitialized()) then
			deallocate(this%dimensionList(:))
			this%numberOfDimensions=0
			this%IndexOfSpecialDimension=0
		endif
		delete_Dimensions=Normal
	end function delete_Dimensions


	logical function Is_Initialized(this)
		type(Dimensions) :: this
		Is_Initialized=allocated(this%dimensionList)
	end	function Is_Initialized


end module Dimensions_Class
