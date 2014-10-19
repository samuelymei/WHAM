module simulation
  use precision_m
  use snapshot
  use react_coord_bin
  use bin
  private
  integer(kind=4), public :: nSimulation

  type :: simulation_info
    real(kind=fp_kind) :: beta  ! inverse temperature of this simulation
    integer(kind=4) :: nSnapshots ! Number of snapshots in this simulation
    integer(kind=4) :: nEffectivenSnapshots ! snapshots outside the RC bins are removed
    real(kind=fp_kind) :: freeEnergy
    type ( snapshot_info ), pointer :: snapshots(:) ! point to the snapshots
    type ( bin_info ), pointer :: bins(:)
  end type simulation_info
  type ( simulation_info ), allocatable, public :: simulations(:) ! to save the information of all the simulations
  real(kind=fp_kind) :: energyMin, energyMax
  public :: readSimulationInfo, deleteSimulationInfo, initBins

contains

  subroutine readSimulationInfo(fid)
    use constant
    implicit none
    integer(kind=4), intent(in) :: fid
    integer(kind=4) :: IndexW, IndexS
    read(fid, *) nSimulation
    allocate(simulations(nSimulation))
    do IndexW = 1, nSimulation
      read(fid, *) simulations(IndexW)%beta, simulations(IndexW)%nSnapshots
      allocate(simulations(IndexW)%snapshots(simulations(IndexW)%nSnapshots))
      do IndexS = 1, simulations(IndexW)%nSnapshots
        read(fid,*) simulations(IndexW)%snapshots(IndexS)%jReactCoordBin, & 
                  & simulations(IndexW)%snapshots(IndexS)%energyUnbiased
      end do
    end do
    call energyMinMax
  end subroutine readSimulationInfo

  subroutine deleteSimulationInfo
    implicit none
    integer(kind=4) :: IndexW
    do IndexW = 1, nSimulation
      deallocate(simulations(IndexW)%snapshots)
    end do
    if(allocated(simulations))deallocate(simulations)
  end subroutine deleteSimulationInfo

  subroutine energyMinMax
    implicit none
    integer(kind=4) :: IndexW, IndexS
    energyMin = simulations(1)%snapshots(1)%energyUnbiased
    energyMax = simulations(1)%snapshots(1)%energyUnbiased
    do IndexW = 1, nSimulation
      do IndexS = 1, simulations(IndexW)%nSnapshots
        if(simulations(IndexW)%snapshots(IndexS)%energyUnbiased < energyMin) & 
           & energyMin = simulations(IndexW)%snapshots(IndexS)%energyUnbiased
        if(simulations(IndexW)%snapshots(IndexS)%energyUnbiased > energyMax) & 
           & energyMax = simulations(IndexW)%snapshots(IndexS)%energyUnbiased
      end do
    end do
  end subroutine energyMinMax

  subroutine initBins(NumRCbin, NumEbin, NumBin, T_target)
    use constant
    implicit none
    integer(kind=4), intent(in) :: NumRCbin, NumEbin
    integer(kind=4), intent(out) :: NumBin
    real(kind=fp_kind) :: T_target

    real(kind=fp_kind) :: energy(NumEbin)
    real(kind=fp_kind) :: deltaRC, deltaE
    integer(kind=4) :: IndexRCbin, IndexEbin
    real(kind=fp_kind) :: DELTA 
    integer(kind=4) :: IndexW, IndexB
    real(kind=fp_kind) :: beta_target
    real(kind=fp_kind) :: beta
    real(kind=fp_kind) :: biasingPotential

    DELTA = (energyMax - energyMin)/(100*NumEbin)
    deltaE = (energyMax - energyMin + DELTA)/NumEbin
    do IndexEbin = 1, NumEbin
      energy(IndexEbin) = energyMin + deltaE * (IndexEbin-0.5)
    end do

    NumBin = NumRCbin * NumEbin
    do IndexW = 1, nSimulation
      allocate(simulations(IndexW)%bins(NumBin))
    end do

    beta_target = 1.d0 / (kB*T_target)
    IndexB = 0
    do IndexRCbin = 1, NumRCbin
      do IndexEbin = 1, NumEbin
        IndexB = IndexB + 1 
        do IndexW = 1, nSimulation
          simulations(IndexW)%bins(IndexB)%energyBiasing = reactCoordBin(IndexRCbin)%energyBiasing
          simulations(IndexW)%bins(IndexB)%energy = energy(IndexEbin)
          simulations(IndexW)%bins(IndexB)%energyWidth = deltaE
          simulations(IndexW)%bins(IndexB)%biasingFactor = &
            & exp(-(simulations(IndexW)%beta-beta_target)*simulations(IndexW)%bins(IndexB)%energy) * &
            & exp(-simulations(IndexW)%beta*simulations(IndexW)%bins(IndexB)%energyBiasing)
          simulations(IndexW)%bins(IndexB)%jReactCoord = IndexRCbin
          simulations(IndexW)%bins(IndexB)%kEnergy = IndexEbin
        end do
      end do
    end do
    call buildHistogram(NumRCbin, NumEbin)
  end subroutine initBins

  subroutine buildHistogram(NumRCbin, NumEbin)
    implicit none
    integer(kind=4), intent(in) :: NumRCbin, NumEbin
    integer(kind=4) :: IndexW, IndexB, IndexS
    integer(kind=4) :: IndexRCbin, IndexEbin
    do IndexW = 1, nSimulation
      simulations(IndexW)%bins(:)%histogram = 0
      do IndexS = 1, simulations(IndexW)%nSnapshots
        if(simulations(IndexW)%snapshots(IndexS)%jReactCoordBin < 0 .or. &
         & simulations(IndexW)%snapshots(IndexS)%jReactCoordBin > NumRCbin ) cycle
        IndexEbin = ( simulations(IndexW)%snapshots(IndexS)%energyUnbiased - &
                     & simulations(IndexW)%bins(IndexB)%energy + &
                     & simulations(IndexW)%bins(IndexB)%energyWidth / 2.d0 ) &
                   & /simulations(IndexW)%bins(IndexB)%energyWidth + 1

        IndexB = ( IndexRCbin - 1 ) * NumEbin + IndexEbin
        simulations(IndexW)%bins(IndexB)%histogram = &
                   & simulations(IndexW)%bins(IndexB)%histogram + 1
      end do
      simulations(IndexW)%nEffectivenSnapshots = sum(simulations(IndexW)%bins(:)%histogram)
    end do
  end subroutine buildHistogram

end module simulation

