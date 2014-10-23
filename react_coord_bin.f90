module react_coord_bin
  use precision_m
  implicit none
  private
  type :: reactCoordBin_info
    real(kind=fp_kind) :: energyBiasing
  end type reactCoordBin_info
  type ( reactCoordBin_info ), allocatable, public :: reactCoordBin(:, :)
  public :: readReactCoordBinInfo, deleteReactCoordBinInfo
contains
  subroutine readReactCoordBinInfo( fid, NumW, NumJ )
    implicit none
    integer(kind=4), intent(in) :: fid, NumW, NumJ
    integer(kind=4) :: indexW
    allocate(reactCoordBin(NumJ, NumW))
    read(fid,*)(reactCoordBin(:, indexW)%energyBiasing, indexW = 1, NumW)
  end subroutine readReactCoordBinInfo

  subroutine deleteReactCoordBinInfo
    implicit none
    deallocate(reactCoordBin)
  end subroutine deleteReactCoordBinInfo
end module react_coord_bin
