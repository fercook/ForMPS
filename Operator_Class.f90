!! Matrix Product States algorithms
!! Author: Fernando M. Cucchietti 2010

!!TODO: Reimplement extendind the type from a Matrix class that performs some low level functions

module Operator_Class

  use ErrorHandling
  use Constants
  use MPSTensor_Class

  implicit none

  integer,parameter :: OPERATOR_BASIS_SIZE = 3

!###############################
!#####  The class main object
!###############################
  type Operator
     private
     complex(8) :: data(MAX_spin, MAX_spin) !!$TODO: Change to extend a matrix (maybe abstract) type
   contains
     procedure :: ApplyTo => ApplyToMPSTensor
  end type Operator

  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaX = Reshape([zero,one,one,zero], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaY = Reshape([zero,II,-II,zero], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaZ = Reshape([one,zero,zero,zero,-one], [MAX_spin,MAX_spin]) 

  type(Operator),dimension(OPERATOR_BASIS_SIZE) :: PauliSigma 
  data PauliSigma/  Operator(PauliSigmaX) ,  Operator(PauliSigmaY) ,  Operator(PauliSigmaZ) /

 contains

   function ApplyToMPSTensor (this,aMPSTensor) result(newMPSTensor)
     TYPEORCLASS(Operator),intent(IN) :: this   !!<<TYPE>>!!
     type(MPSTensor),intent(IN) :: aMPSTensor
     type(MPSTensor) :: newMPSTensor
     
     newMPSTensor = aMPSTensor%ApplyOperator(this)

   end function ApplyToMPSTensor

 end module Operator_Class
