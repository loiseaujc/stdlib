! Demonstrate usage of `getfile`
program example_getfile
  use stdlib_io, only: getfile
  use stdlib_string_type, only: string_type
  use stdlib_error, only: state_type
  implicit none

  character(*), parameter :: filename = "example.txt"
  type(string_type) :: filecontent
  type(state_type) :: err

  ! Read a file into a string
  call getfile(fileName, fileContent, err=err)

  if (err%error()) then
    print *, err%print()
  else
    print *, "Success! File "//fileName//" imported."
  end if
end program example_getfile
