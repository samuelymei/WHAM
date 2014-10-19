module snapshot
  use precision_m
  implicit none
  public 
  integer(kind=4) :: nSnapshots
  type :: snapshot_info
    real(kind=fp_kind) :: energyUnbiased ! Energy
    integer(kind=4) :: jReactCoordBin    ! Index of RC bin
    integer(kind=4) :: kEnergyBin        ! Index of Energy bin
    integer(Kind=4) :: iUnifiedBin       ! Index of unified bin
  end type snapshot_info
end module snapshot


