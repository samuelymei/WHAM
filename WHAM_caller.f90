program WHAM_caller
  use WHAM
  implicit none
  integer(kind=4) :: fid
  character(len=60) :: datafile
  fid = 10
  open(10, file = datafile)
  call startWHAM(fid)
  close(10)
end program WHAM_caller
