!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the calculation of the free energy from biased simulation like Umbrella sampling !
! and Temperature replica exchange molecular dynamics simualtion                              !
! Written by                                                                                  !
!                                        Ye Mei                                               !
!                                      10/17/2014                                             !
!                           East China Normal University                                      !
! Reference:                                                                                  !
!      J. Phys. Chem. B, 109, 6722-6731 (2005)                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
module WHAM
  use precision_m
  use simulation
  use react_coord_bin
  implicit none
  
  private

  integer(kind=4), public :: NumJ = 50  ! number of reaction coordinate bins
  integer(kind=4), public :: NumK = 100 ! number of energy bins
  integer(kind=4), public :: NumB 
  real(kind=fp_kind), parameter :: TOLERANCE = 1.0e-5 ! convergence criterion for free energy 
  integer(kind=4), parameter :: MaxITS = 1000  ! max number of iterations
  real(kind=fp_kind), public :: T_target = 300 ! target temperature

  real(kind=fp_kind) :: beta ! inverse temperature

  integer(kind=4) :: IndexW ! index of windows/simulations
  integer(kind=4) :: IndexS ! index of snapshot
  integer(kind=4) :: IndexK ! index of energy bin
  integer(kind=4) :: IndexJ ! index of reaction coordinate bin

  integer(kind=4) :: IndexB ! index of combined Bin

  real(kind=fp_kind), allocatable :: unbiasedDensity(:)
  
  public :: startWHAM
contains
  subroutine iteration
    implicit none
    real(kind=fp_kind) :: numerator(NumB), denominator(NumB)
    real(kind=fp_kind) :: freeenergyMin
    real(kind=fp_kind) :: freeenergyRMSD
    real(kind=fp_kind) :: freeenergyOld(nSimulation)
    logical :: converged

    simulations(:)%freeenergy = 1.d0   ! assign an initial guess of the free energy
    converged = .false.

    do while (.not. converged)
      freeenergyOld = simulations(:)%freeenergy
      do IndexB = 1, NumB
        numerator(IndexB) = 0.d0
        do IndexW = 1, nSimulation
          numerator(IndexB) = numerator(IndexB) + & 
             & simulations(IndexW)%bins(IndexB)%histogram
        end do
        denominator(IndexB) = 0.d0
        do IndexW = 1, nSimulation
          denominator(IndexB) = denominator(IndexB) + &
             & simulations(IndexW)%nEffectivenSnapshots * simulations(IndexW)%freeenergy * &
             & simulations(IndexW)%bins(IndexB)%biasingFactor
        end do
      end do
      unbiasedDensity = numerator / denominator
 
      simulations(:)%freeenergy = 0.d0
      do IndexW = 1, nSimulation
        do IndexB = 1, NumB
          simulations(IndexW)%freeenergy = simulations(IndexW)%freeenergy + &
            & simulations(IndexW)%bins(IndexB)%biasingFactor * unbiasedDensity(IndexB)
        end do    
      end do
      freeenergyMin = minval(simulations(:)%freeenergy) - 1
      simulations(:)%freeenergy = simulations(:)%freeenergy - freeenergyMin
      simulations(:)%freeenergy = 1.d0/ simulations(:)%freeenergy

      freeenergyRMSD =0.d0
      do IndexW = 1, nSimulation
        freeenergyRMSD = freeenergyRMSD + & 
          & (simulations(IndexW)%freeenergy - freeenergyOld(IndexW))**2
      end do
      freeenergyRMSD = sqrt(freeenergyRMSD/IndexW)
      if( freeenergyRMSD < TOLERANCE ) converged = .true.
    end do
  end subroutine iteration

  subroutine startWHAM(fid)
    implicit none
    integer(kind=4), intent(in) :: fid
    call readReactCoordBinInfo(fid, NumJ)
    call readSimulationInfo(fid)
    call initBins(NumJ, NumK, NumB, T_target)
    call iteration
    call deleteReactCoordBinInfo
    call deleteSimulationInfo
  end subroutine startWHAM
end module WHAM
