!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Container module for the Slater-Koster data
!>
!> This module contains the Slater-Koster tables. It decides, which one to call for which
!> species. It can be easily extended to contain different Slater-Koster schemes for different
!> species. At the moment, it handles only Slater-Koster data tabulated on an equidistant grid.
module slakocont
  use assert
  use accuracy
  use slakoeqgrid
  implicit none
  private

  public :: OSlakoCont, init, destruct
  public :: addTable, addVecTable, getMIntegrals, getCutoff, getSKIntegrals
  public :: getMPOrigins


  !> A specific Slater-Koster table implementation.
  type PSlaKo_
    integer :: iType = 0
    type(OSlakoEqGrid), allocatable :: pSlakoEqGrid
  end type PSlaKo_


  !> Container for Slater-Koster integrals for all pair-interactions
  type OSlakoCont
    private
    type(PSlaKo_), allocatable :: slakos(:,:)
    type(PMPoleOrigin_), pointer :: dipoleOrigin(:,:)
    integer :: nSpecies
    integer :: mInt
    real(dp) :: cutoff
    logical :: tDataOK
    logical :: tInit = .false.
  end type OSlakoCont

  !! Pointer to the table containing the origins of a multipole vector
  type PMPoleOrigin_
    integer :: iType = 0
    type(OOrigins), pointer :: pOrigin => null()
  end type PMPoleOrigin_

  !! Pointer to the table containing on-site integrals
  ! type POnSites_
  !   integer :: iType = 0
  !   type(OOrigins), pointer :: pOnSites => null()
  ! end type POnSites_

  !> Initialises SlakoCont
  interface init
    module procedure SlakoCont_init
  end interface init

  !!* Destroys the components of SlakoCont
  interface destruct
    module procedure SlakoCont_destruct
  end interface

  !> Adds a Slater-Koster table for a given diatomic pair to the container.
  interface addTable
    module procedure SlakoCont_addTableEqGrid
  end interface addTable

  !!* Adds a origin vector table for a given diatomic pair to the container.
  interface addVecTable
    module procedure SlakoCont_addVecTable
  end interface addVecTable


  !> Returns the maximal number of integrals needed for the interactions.
  interface getMIntegrals
    module procedure SlakoCont_getMIntegrals
  end interface getMIntegrals


  !> Returns the cutoff for all interactions
  interface getCutoff
    module procedure SlakoCont_getCutoff
  end interface getCutoff


  !> Returns the Slater-Koster integrals for a given distance for a given species pair.
  interface getSKIntegrals
    module procedure SlakoCont_getSKIntegrals
  end interface getSKIntegrals

  !!* Returns the origin vector table of the multipole expansion for a given
  !!* species pair.
  interface getMPOrigins
    module procedure SlakoCont_getMPOrigins
  end interface getMPOrigins

contains


  !> Initialises SlakoCont
  subroutine SlakoCont_init(self, nSpecies)

    !> SlakoCont instance
    type(OSlakoCont), intent(out) :: self

    !> Nr. of species in the system.
    integer, intent(in) :: nSpecies

    @:ASSERT(.not. self%tInit)

    self%nSpecies = nSpecies
    allocate(self%slakos(nSpecies, nSpecies))
    allocate(self%dipoleOrigin(nSpecies, nSpecies))
    ! INITALLOCATE_PARR(self%onSites, (nSpecies, nSpecies))
    self%mInt = 0
    self%cutoff = 0.0_dp
    self%tDataOK = .false.
    self%tInit = .true.

  end subroutine SlakoCont_init

  
  !!* Destroys the components of SlakoCont
  !!* @param self SlakoCont instance
  subroutine SlakoCont_destruct(self)
    type(OSlakoCont), intent(inout) :: self

    deallocate(self%slakos)
    deallocate(self%dipoleOrigin)
    ! DEALLOCATE_PARR(self%onSites)
    self%tInit = .false.
    self%tDataOK = .false.

  end subroutine SlakoCont_destruct

  !> Adds a Slater-Koster table for a given diatomic pair to the container.
  subroutine SlakoCont_addTableEqGrid(self, pTable, iSp1, iSp2)

    !> SlakoCont instance
    type(OSlakoCont), intent(inout) :: self

    !> Slater-Koster table to be added
    type(OSlakoEqGrid), allocatable, intent(inout) :: pTable

    !> Index of the first interacting species
    integer, intent(in) :: iSp1

    !> Index of the second interacting species
    integer, intent(in) :: iSp2

    @:ASSERT(self%tInit)
    self%mInt = max(self%mInt, getNIntegrals(pTable))
    self%cutoff = max(self%cutoff, getCutoff(pTable))
    self%slakos(iSp2, iSp1)%iType = 1
    call move_alloc(pTable, self%slakos(iSp2, iSp1)%pSlakoEqGrid)
    self%tDataOK = all(self%slakos(:,:)%iType /= 0)

  end subroutine SlakoCont_addTableEqGrid



  !!* Adds a vector table for a given diatomic pair to the container.
  !!* @param self SlakoCont instance
  !!* @param pTable Pointer to the vector table to be added
  !!* @param iSp1 Index of the first interacting species
  !!* @param iSp2 Index of the second interacting species
  subroutine SlakoCont_addVecTable(self, pOrigins, iSp1, iSp2)
    type(OSlakoCont), intent(inout) :: self
    type(OOrigins),  allocatable :: pOrigins
    integer, intent(in) :: iSp1, iSp2

    @:ASSERT(self%tInit)
    self%dipoleOrigin(iSp2, iSp1)%iType = 1
    self%dipoleOrigin(iSp2, iSp1)%pOrigin = pOrigins
    self%tDataOK = all(self%dipoleOrigin(:,:)%iType /= 0)
    self%mInt = max(self%mInt, getNIntegrals(pOrigins))

  end subroutine SlakoCont_addVecTable


  !> Returns the maximal number of integrals needed for describing any of the interactions in the
  !> container
  !>
  !> This subroutine is "pure", so that it can be used to determine the size of static arrays.
  pure function SlakoCont_getMIntegrals(self) result(mInt)

    !> SlakoCont instance
    type(OSlakoCont), intent(in) :: self

    !> Max. number of integrals.
    integer :: mInt

    !! Pure procedures can not contain any I/O, therefore the following assertion is commented out
    !@:ASSERT(self%tInit)
    mInt = self%mInt

  end function SlakoCont_getMIntegrals


  !> Returns the cutoff for all interactions
  function SlakoCont_getCutoff(self) result(cutoff)

    !> SlakoCont instance
    type(OSlakoCont), intent(in) :: self

    !> Cutoff of interaction
    real(dp) :: cutoff

    @:ASSERT(self%tInit)
    cutoff = self%cutoff

  end function SlakoCont_getCutoff


  !> Returns the Slater-Koster integrals for a given distance for a given species pair.
  subroutine SlakoCont_getSKIntegrals(self, sk, dist, sp1, sp2)

    !> SlakoCont instance
    type(OSlakoCont), intent(in) :: self

    !> Contains the integrals on exit
    real(dp), intent(out) :: sk(:)

    !> Distance of the two atoms
    real(dp), intent(in) :: dist

    !> Index of the first interacting species.
    integer, intent(in) :: sp1

    !> Index of the second interacting species.
    integer, intent(in) :: sp2

    @:ASSERT(self%tInit .and. self%tDataOK)
    call getSKIntegrals(self%slakos(sp2, sp1)%pSlakoEqGrid, sk, dist)

  end subroutine SlakoCont_getSKIntegrals

  !!* Returns the on-site integrals for a given species.
  !!* @param self SlakoCont instance
  !!* @param onSites Contains the integrals on exit
  !!* @param sp Index of the species.
  ! subroutine SlakoCont_getOnSites(self, onSites, sp)
  !   type(OSlakoCont), intent(in) :: self
  !   real(dp), intent(out) :: onSites(:)
  !   integer, intent(in) :: sp

  !   ASSERT(self%tInit .and. self%tDataOK)
  !   call getOnSites(self%slakos(sp, sp)%pOnSites, onSites)

  ! end subroutine SlakoCont_getOnSites



  !!* Returns the multipole expansion origin vectors (one for each interaction)
  !!* of a atom pair for a given distance.
  !!* @param self SlakoCont instance
  !!* @param origVecs Contains the integrals on exit
  !!* @param dist Distance of the two atoms
  !!* @param sp1 Index of the first interacting species.
  !!* @param sp2 Index of the second interacting species.
  subroutine SlakoCont_getMPOrigins(self, origVecs, centers, dist, sp1, sp2)
    type(OSlakoCont), intent(in) :: self
    real(dp), intent(out) :: origVecs(:,:)
    integer, intent(out) :: centers(:)
    real(dp), intent(in) :: dist
    integer, intent(in) :: sp1, sp2

    @:ASSERT(self%tInit .and. self%tDataOK)
    call getMPOrigins(self%dipoleOrigin(sp2, sp1)%pOrigin, origVecs, centers, &
        &dist)

  end subroutine SlakoCont_getMPOrigins


end module slakocont
