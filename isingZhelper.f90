!!   Copyright 2011 Fernando M. Cucchietti
!
!    This file is part of ForMPS
!
!    ForMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ForMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ForMPS.  If not, see <http://www.gnu.org/licenses/>.

module isingZhelper

   use PEPS_Class
   use PEPOTensor_Class
   use PEPSTensor_Class
   use Constants
   use ErrorHandling

    implicit none

contains

    function DirectionalPEPOIsing(direction,beta) result(aPepo)
      integer,intent(IN) :: direction
      real(8),intent(IN) :: beta
      type(PEPOTensor) :: aPEPO
      integer,parameter :: bondDim = 2, spinDim=2
      real(8) :: sqSinh,sqCosh
      complex(8) :: matrix(2,2,2)
      complex(8),allocatable :: tensor(:,:,:,:,:,:)

      sqSinh=Sqrt(Sinh(beta))
      sqCosh=sqrt(cosh(beta))
      matrix=ZERO
      matrix(1,1,1)= ONE*sqCosh  !bond index is last, spin indices are first
      matrix(2,2,1)= ONE*sqCosh
      matrix(1,1,2)= II*sqSinh
      matrix(2,2,2)= -II*sqSinh

      select case (direction)
         case (LEFT)
            allocate(tensor(spinDim,spinDim,bondDim,integerONE,integerONE,integerONE))
            tensor(:,:,:,1,1,1)=matrix
         case (RIGHT)
            allocate(tensor(spinDim,spinDim,integerONE,bondDim,integerONE,integerONE))
            tensor(:,:,1,:,1,1)=matrix
         case (UP)
            allocate(tensor(spinDim,spinDim,integerONE,integerONE,bondDim,integerONE))
            tensor(:,:,1,1,:,1)=matrix
         case (DOWN)
            allocate(tensor(spinDim,spinDim,integerONE,integerONE,integerONE,bondDim))
            tensor(:,:,1,1,1,:)=matrix
         case default
            call ThrowException('DirectionalPEPOIsing','Unknown direction',direction,CriticalError)
      end select

      aPEPO=new_PEPOTensor(tensor)

    end function DirectionalPEPOIsing



   subroutine preparePEPS(thePEPS,Xlength,Ylength,beta)
      integer,intent(IN) :: Xlength,Ylength
      real(8) :: beta
      type(PEPOTensor) :: PEPOinDir(LEFT:DOWN)
      type(PEPS) :: thePEPS
      type(PEPSTensor) :: tempPEPS
      integer :: x,y,dir
      integer,parameter :: spinDim=2

      do dir=LEFT,DOWN
         PEPOinDir(dir)=DirectionalPEPOIsing(dir,beta/2.0d0)
      enddo

      thePEPS=new_PEPS(Xlength,Ylength,spinDim,integerONE)

      do y=1,Ylength
         do x=1,XLength
            tempPEPS=new_Tensor( integerONE,integerONE,integerONE,integerONE,2,ONE/sqrt(2.0d0) )
            if (x.gt.1) tempPEPS=PEPOinDir(LEFT).applyTo.tempPEPS
            if (x.lt.Xlength) tempPEPS=PEPOinDir(RIGHT).applyTo.tempPEPS
            if (y.gt.1) tempPEPS=PEPOinDir(DOWN).applyTo.tempPEPS
            if (y.lt.Ylength) tempPEPS=PEPOinDir(UP).applyTo.tempPEPS
            call thePEPS%SetTensorAt(x,y,tempPEPS)
         enddo
      enddo

   end subroutine

end module isingZhelper
