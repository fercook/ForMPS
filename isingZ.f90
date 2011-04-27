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

program isingZ


   use PEPS_Class
   use PEPOTensor_Class
   use PEPSTensor_Class
   use isingZhelper
   use Constants
   use ErrorHandling

   implicit none

   integer,parameter :: steps=10, Xlength=6,Ylength=6
   integer,parameter :: canonX=3,canonY=3
   integer,parameter :: LongBond=10, TransBond=6
   real(8),parameter :: betaMax=2.0d0, betaMin=0.2D0
   type(PEPS) :: thePEPS
   integer :: nbeta
   real(8) :: beta, theNorm

   do nbeta=0,steps
      beta=(betaMax-betaMin)*nbeta/steps + betaMin
      call preparePEPS(thePEPS,Xlength,Ylength,beta)
      call thePEPS%CanonizeAt(canonX,canonY, HORIZONTAL, LongBond, TransBond, theNorm)
      print *, beta, theNorm
   enddo


end program isingZ
