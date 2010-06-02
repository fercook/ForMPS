!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010

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


!###############################
!#####  The class main object
!###############################
  type Hamiltonian
     private
     integer :: Length
!   contains
!     procedure :: delete => delete_MPSTensor
  end type Hamiltonian


 contains

end module Hamiltonian_Class
