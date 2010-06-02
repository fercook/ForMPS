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
  type SpinOperator
     private
     complex(8) :: data(MAX_spin, MAX_spin) !!$TODO: Change to extend a matrix (maybe abstract) type
   contains
     procedure :: ApplyToTensor => ApplyToMPSTensor
  end type SpinOperator

  interface operator (.applyto.) 
     module procedure ApplyToMPSTensor
  end interface

  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaXMatrix = Reshape([zero,one,one,zero], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaYMatrix = Reshape([zero,II,-II,zero], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: PauliSigmaZMatrix = Reshape([one,zero,zero,-one], [MAX_spin,MAX_spin]) 
  complex(8),parameter,dimension(MAX_spin,MAX_spin):: IdentityOperatorMatrix = Reshape([one,zero,zero,one], [MAX_spin,MAX_spin]) 

  type(SpinOperator),dimension(0:OPERATOR_BASIS_SIZE),target :: PauliSigma 
  data PauliSigma/ SpinOperator(IdentityOperatorMatrix), SpinOperator(PauliSigmaXMatrix) ,  SpinOperator(PauliSigmaYMatrix) ,  SpinOperator(PauliSigmaZMatrix) /

!################ Helper classes contained inside Hamiltonians
  type :: SiteTerm
     private
     type(SpinOperator),pointer :: Op
     integer :: site
     complex(8) :: amplitude
  end type SiteTerm

  type :: BondTerm
     private
     type(SpinOperator),pointer :: OpL,OpR
     integer :: siteL, siteR
     complex(8) :: amplitude
  end type BondTerm


  type, extends(SiteTerm) :: SiteTermList
     private
     type(SiteTermList),pointer :: next
  end type SiteTermList


!################# Methods

 contains


   function ApplyToMPSTensor (this,aMPSTensor) result(newMPSTensor)
     TYPEORCLASS(SpinOperator),intent(IN) :: this   !!<<TYPE>>!!
     type(MPSTensor),intent(IN) :: aMPSTensor
     type(MPSTensor) :: newMPSTensor
     
     newMPSTensor = aMPSTensor%ApplyOperator(this%data)

   end function ApplyToMPSTensor



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


 end module Operator_Class
