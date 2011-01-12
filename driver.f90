module Tensor_Class_fun

 use Tensor_Class

implicit none

private

 contains

 subroutine type_creation_deletion

  type(tensor3) :: aTensor
  integer :: dims(3)

  CALL random_seed()
  aTensor=new_Tensor(2,20,20)
  dims=aTensor%GetDimensions()

 end subroutine type_creation_deletion

end module Tensor_Class_fun
