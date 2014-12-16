program WHAM_caller
  use precision_m
  use WHAM
  implicit none
  integer(kind=4) :: narg, iarg
  character(len=20) :: arg
  integer(kind=4) :: fid
  character(len=60) :: metafile
  real(kind=fp_kind) :: tolerance = 1.D-4
  real(kind=fp_kind) :: temperature = 300.0
  character(len=200) :: command_help
  character(len=200):: thisProgram
  character(len=256) :: cmd
  logical :: enoughArg
  call get_command(cmd)
  call getarg(0, thisProgram)
  command_help = 'Usage: '//trim(thisProgram)//" -T temperature -t tolerance -f metafile"
  narg = iargc()
  iarg = 1
  do while (iarg <= narg)
    call getarg(iarg, arg)
    if( arg == '-T' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call getarg(iarg, arg)
      read(arg,*) temperature
    else if( arg == '-t' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call getarg(iarg, arg)
      read(arg,*) tolerance
    else if( arg == '-f' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call getarg(iarg, arg)
      metafile = arg
    else
      write(*,'(A)') command_help
      stop
    end if
    iarg = iarg + 1
  end do
  write(*,'(A,F10.3)')'Target temperature:', temperature
  write(*,'(A,E10.3)')'Tolerance:', tolerance
  write(*,'(A,A)') 'Metafile: ', metafile
  fid = 10
  open(fid, file = metafile, status = 'OLD')
  call startWHAM(fid, temperature, tolerance)
  call finalizeWHAM
  close(fid)
end program WHAM_caller

function enoughArg( narg, nargInNeed )
  implicit none
  logical :: enoughArg
  integer(kind=4), intent(in) :: narg, nargInNeed
  enoughArg = .true.
  if( narg < nargInNeed ) enoughArg = .false.
  return
end function enoughArg
