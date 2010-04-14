module ErrorHandling
  
  integer,parameter :: MaxErrorAllowed=0
  integer,parameter :: Warning = -1
  integer,parameter :: Normal = 0
  integer,parameter :: MinorError = 1
  integer,parameter :: CriticalError = 2
  integer,parameter :: NoErrorCode = 0
  
contains
  
  subroutine ThrowException(caller,message,errorcode,errordegree)
    integer,intent(IN) :: errorcode,errordegree
    character*(*),intent(IN) :: caller,message

    print *,'####### WARNING ####### '
    print *,'Error from :',caller
    print *,'Says       :',message
    print *,'Sent error :',errorcode
    if (errordegree.gt.MaxErrorAllowed) then
       print *,'####### ABORTING ####### '
       stop
    else
       print *,'Continuing...'
       print *,'####### WARNING ####### '
    endif
  end subroutine ThrowException

end module ErrorHandling
