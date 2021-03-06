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

!!  This module contains a class that describes Hamiltonians, in a form suitable
!!  for MPS algorithms.

!!TODO: Reimplement as extension of tensor type
!!TODO: Reimplement with a Visitor pattern that travels around the interactions

module Hamiltonian_Class

  use ErrorHandling 
  use Constants
  use Matrix_Helper
  use Operator_Class

  implicit none

  integer,parameter :: MAX_spin = 2

!################ Operators con know to which site they are applied, and with which amplitude
  type :: SiteTerm
     private
     type(SpinOperator) :: Op
     integer :: site
     complex(8) :: amplitude
  end type SiteTerm

  type :: BondTerm
     private
     type(SpinOperator) :: OpA,OpB
     integer :: siteA, siteB
     complex(8) :: amplitude
  end type BondTerm

!###############################
!#####  The class main object
!###############################
  type Hamiltonian
     private
     integer :: Length
     integer :: InteractionRange

!   contains
!     procedure,public :: addTerm => addTermToHamiltonian
!     procedure,public :: FieldsAtSite =>
!     procedure,public :: InteractionsAtSite =>
  end type Hamiltonian


!THIS class must provide "Iterators", or lists of terms that involve a given site.
!Example:
!FieldTermsOfSite(i)%next() and so on
!BondTermsOfSite(i)%ExistsNext()
!
!Look at common terminology for Iterators
! Also this class must contain a function like
! AddBondTerm(Op,i,j)
 contains


   function new_SiteTerm(Op,site,amplitude) result(this)
     type(SiteTerm) :: this
     type(SpinOperator),target :: Op
     integer :: site
     complex(8) :: amplitude

     this%Op => Op
     this%site = site
     this%amplitude = amplitude
   end function new_SiteTerm



   function new_SiteTermList(Op,site,amplitude) result(this)
     type(SiteTermList) :: this
     type(SpinOperator),target :: Op
     integer :: site
     complex(8) :: amplitude
     this%Op => Op
     this%site = site
     this%amplitude = amplitude
     this%next => null()
   end function new_SiteTermList



end module Hamiltonian_Class
