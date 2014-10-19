module react_coord_bin
  use precision_m
  implicit none
  private
  type :: reactCoordBin_info
    real(kind=fp_kind) :: energyBiasing
  end type reactCoordBin_info
  type ( reactCoordBin_info ), allocatable, public :: reactCoordBin(:)
  public :: readReactCoordBinInfo, deleteReactCoordBinInfo
contains
  subroutine readReactCoordBinInfo( fid, NumJ )
    implicit none
    integer(kind=4), intent(in) :: fid, NumJ
    allocate(reactCoordBin(NumJ))
    read(fid,*)reactCoordBin(:)%energyBiasing
  end subroutine readReactCoordBinInfo

  subroutine deleteReactCoordBinInfo
    implicit none
    deallocate(reactCoordBin)
  end subroutine deleteReactCoordBinInfo
end module react_coord_bin
