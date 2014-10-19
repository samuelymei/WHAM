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
    integer(kind=4) :: IndexW
    allocate(reactCoordBin(NumJ, NumW))
    do IndexW = 1, NumW
      read(fid,*)reactCoordBin(:, IndexW)%energyBiasing
    end do
  end subroutine readReactCoordBinInfo

  subroutine deleteReactCoordBinInfo
    implicit none
    deallocate(reactCoordBin)
  end subroutine deleteReactCoordBinInfo
end module react_coord_bin
