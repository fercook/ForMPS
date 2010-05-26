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

  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaXMatrix = Reshape([zero,one,one,zero], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaYMatrix = Reshape([zero,II,-II,zero], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaZMatrix = Reshape([one,zero,zero,zero,-one], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: IdentityOperatorMatrix = Reshape([one,zero,zero,zero,one], [MAX_spin,MAX_spin]) 

  type(Operator),dimension(0:OPERATOR_BASIS_SIZE) :: PauliSigma 
  data PauliSigma/ Operator(IdentityOperatorMatrix), Operator(PauliSigmaXMatrix) ,  Operator(PauliSigmaYMatrix) ,  Operator(PauliSigmaZMatrix) /

 contains

   function ApplyToMPSTensor (this,aMPSTensor) result(newMPSTensor)
     TYPEORCLASS(Operator),intent(IN) :: this   !!<<TYPE>>!!
     type(MPSTensor),intent(IN) :: aMPSTensor
     type(MPSTensor) :: newMPSTensor
     
     newMPSTensor = aMPSTensor%ApplyOperator(this%data)

   end function ApplyToMPSTensor

 end module Operator_Class
