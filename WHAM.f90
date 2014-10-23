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
  integer(kind=4), public :: NumW  ! number of simulations
  integer(kind=4), public :: NumJ  ! number of reaction coordinate bins
  integer(kind=4), public :: NumK  ! number of energy bins
  integer(kind=4), public :: NumB 
  real(kind=fp_kind), parameter :: TOLERANCE = 1.0e-6 ! convergence criterion for free energy 
  integer(kind=4), parameter :: MaxITS = 1000  ! max number of iterations
  real(kind=fp_kind), public :: T_target = 300 ! target temperature

  real(kind=fp_kind) :: beta ! inverse temperature

  integer(kind=4) :: indexW ! index of windows/simulations
  integer(kind=4) :: indexS ! index of snapshot
  integer(kind=4) :: indexK ! index of energy bin
  integer(kind=4) :: indexJ ! index of reaction coordinate bin

  integer(kind=4) :: indexB ! index of combined Bin

  real(kind=fp_kind), allocatable :: unbiasedDensity(:)
  
  public :: startWHAM, finalizeWHAM
contains
  subroutine iteration
    implicit none
    real(kind=fp_kind) :: numerator(NumB), denominator(NumB)
    real(kind=fp_kind) :: freeenergyMin
    real(kind=fp_kind) :: freeenergyRMSD
    real(kind=fp_kind) :: freeenergyOld(nSimulation)
    integer(kind=4) :: iIteration
    logical :: converged

    allocate(unbiasedDensity(NumB))
    simulations(:)%freeenergy = 1.d0   ! assign an initial guess of the free energy
    converged = .false.
    do iIteration = 1, MAXITS
      freeenergyOld = simulations(:)%freeenergy
      do indexB = 1, NumB
        numerator(indexB) = 0.d0
        do indexW = 1, nSimulation
          numerator(indexB) = numerator(indexB) + & 
             & simulations(indexW)%bins(indexB)%histogram
        end do
        denominator(indexB) = 0.d0
        do indexW = 1, nSimulation
          denominator(indexB) = denominator(indexB) + &
             & simulations(indexW)%nEffectivenSnapshots * simulations(indexW)%freeenergy * &
             & simulations(indexW)%bins(indexB)%biasingFactor
        end do
      end do
      unbiasedDensity = numerator / denominator
 
      simulations(:)%freeenergy = 0.d0
      do indexW = 1, nSimulation
        do indexB = 1, NumB
          simulations(indexW)%freeenergy = simulations(indexW)%freeenergy + &
            & simulations(indexW)%bins(indexB)%biasingFactor * unbiasedDensity(indexB)
        end do    
      end do
      freeenergyMin = minval(simulations(:)%freeenergy) - 1
      simulations(:)%freeenergy = simulations(:)%freeenergy - freeenergyMin
      simulations(:)%freeenergy = 1.d0/ simulations(:)%freeenergy

      freeenergyRMSD =0.d0
      do indexW = 1, nSimulation
        freeenergyRMSD = freeenergyRMSD + & 
          & (simulations(indexW)%freeenergy - freeenergyOld(indexW))**2
      end do
      freeenergyRMSD = sqrt(freeenergyRMSD/indexW)
      write(6,'(A,1X,I4,A,E12.4)')'Iteration ', iIteration, ': RMSD of dimensionless free energy:', freeenergyRMSD
      if( freeenergyRMSD < TOLERANCE ) converged = .true.
      if( converged ) then
        exit
      end if
    end do
  end subroutine iteration

  subroutine startWHAM(fid)
    implicit none
    integer(kind=4), intent(in) :: fid
    read(fid,*)NumW, NumJ, NumK, T_target
    write(6,*)'Number of simulations:', NumW
    write(6,*)'Number of reaction coordinate bins:', NumJ
    write(6,*)'Number of Energy bins:', NumK
    write(6,*)'Target temperature:', T_target
    print*,'call readReactCoordBinInfo'
    call readReactCoordBinInfo(fid, NumW, NumJ)
    nSimulation = NumW
    print*,'call readSimulationInfo'
    call readSimulationInfo(fid)
    print*,'call initBins'
    call initBins(NumJ, NumK, NumB, T_target)
    print*,'iteration'
    call iteration
  end subroutine startWHAM

  subroutine finalizeWHAM
    implicit none
    deallocate(unbiasedDensity)
    call deleteReactCoordBinInfo
    call deleteSimulationInfo
  end subroutine finalizeWHAM
end module WHAM
