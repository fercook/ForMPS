! Module for standard error handling
! it exports constants that describe normal behavior or slight or
!  critical errors. Main routine is ThrowException, that is called
!   when an error is found, it informs the user, and then evaluates
!    if aborting the program is necessary.

module ErrorHandling
  
  integer :: MaxErrorAllowed=0

  integer,parameter :: Warning = -1
  integer,parameter :: Normal = 0
  integer,parameter :: MinorError = 1
  integer,parameter :: CriticalError = 2
  integer,parameter :: NoErrorCode = 0

  logical :: VerboseLevel = .false.

  logical :: ErrorFlagged = .false.
  integer :: CurrentErrorLevel = Normal

contains
  
  subroutine ThrowException(caller,message,errorcode,errordegree)
    integer,intent(IN) :: errorcode,errordegree
    character*(*),intent(IN) :: caller,message

    print *,'####### WARNING ####### '
    print *,'Error from :',caller
    print *,'Says       :',message
    print *,'Sent error :',errorcode
    print *,'####  END WARNING  #### '
    CurrentErrorLevel = errordegree
    ErrorFlagged=.true.

  end subroutine ThrowException

  subroutine LowerFlag()
    if(VerboseLevel) print *,'Lowering error flag'
    ErrorFlagged=.false.
  end subroutine LowerFlag

  subroutine  ProcessException(caller)
    character*(*),intent(IN) :: caller
    if (CurrentErrorLevel.gt.MaxErrorAllowed) then
       print *,'####### ABORTING ####### '
       stop
    else
       print *,'####### CONTINUING ####### '
       ErrorFlagged=.false.
    endif
  end subroutine ProcessException

  logical function WasThereError() result(answer)
    if(VerboseLevel) print *,'Error flagged = ',Errorflagged
    answer=ErrorFlagged
    return
  end function WasThereError

end module ErrorHandling
