!!* Routines implementing the multipole expansion.
module Multipole
#:include 'common.fypp'
  use Accuracy
  use blasRoutines
  use CommonTypes
  use coulomb
  use Message
  use Periodic, only : TNeighbourList, getNrOfNeighbours
  use shortgamma
  implicit none
  private

  public :: TMPoleInit, OMultipole
  public :: init, getCutoff, updateCharges, updateCoords, updateDipoles
  public :: getEnergy, getShifts, getEnergyPerAtom, addGradientDC, addNumGradientDC

  !!* Input for the multipole module (information that comes from outside the
  !!* module to initialize it).
  type TMPoleInit
    type(TOrbitals) :: orb ! Information about orbitals
    real(dp), allocatable :: hubbU(:,:) ! hubb U for each atom
    !! The following variables are currently given only as a stopping mechanism:
    logical :: tThirdOrder, tGeoOpt, tOrbResolved, tPeriodic
  end type TMPoleInit



  !!* Internal variables of the multipole module.
  !!* @note Comments on notations:
  !!* 10, 11, 20: Notation convention for derivatives with respect to
  !!* coordinates of atom 1 or 2, see multipoles.pdf in documentation.
  type OMultipole
    integer :: nSpecies, nAtom, mShell
    real(dp) :: mCutoff ! Maximum cutoff value
    real(dp), allocatable :: UU(:) ! Hubbards and Derivatives
    real(dp), allocatable :: charges(:) ! Monopole charges
    real(dp), allocatable :: cutoffs(:,:) ! Gamma cutoffs
    integer, allocatable :: nNeigh(:,:), nNeighMax(:) ! Neighbor list
    integer, allocatable :: nNeighShort(:,:) ! Short range gamma neighbors
    real(dp), allocatable :: deltaDAtom(:,:) ! Negative dipole per atom
    !! Derivatives of the short-range part of gamma
    real(dp), allocatable :: shortGamma01(:,:,:), shortGamma10(:,:,:)
    real(dp), allocatable :: shortGamma20(:,:,:,:), shortGamma11(:,:,:,:)
    !! Coulomb potential derivative matrices.
    real(dp), allocatable :: invRMat01(:,:,:), invRMat10(:,:,:)
    real(dp), allocatable :: invRMat20(:,:,:,:), invRMat11(:,:,:,:)
    !! Complete gamma derivatives.
    real(dp), allocatable :: gamma01(:,:,:), gamma10(:,:,:)
    real(dp), allocatable :: gamma20(:,:,:,:), gamma11(:,:,:,:)
    !! Shift vectors:
    !!    o shift1: (gamma01_AC + gamma01_BC) . Delta d_C (scalar per atom)
    !!    o shift2: gamma10_AC . Delta q_C (3 dim per atom)
    real(dp), allocatable :: shift1(:), shift2(:,:) ! Gamma shifts
    logical :: tInitialised = .false. ! Initialization status of the module
    logical :: tCoordUp = .false. ! Update status of coords in geometry step
    logical :: tChargeUp = .false. ! Update status of charges in SCC step
    logical :: tDipoleUp = .false. ! Update status of dipoles in SCC step
  end type OMultipole



  !!* Initialization of the multipole module.
  interface init
    module procedure Multipole_init
  end interface init

  !!* Returns the maximum short range interaction cutoff.
  interface getCutoff
    module procedure Multipole_getCutoff
  end interface getCutoff

  !!* Updates the monopole charges for the module.
  interface updateCharges
    module procedure Multipole_updateCharges
  end interface updateCharges

  !!* Updates the coordinates of the multipole module.
  interface updateCoords
    module procedure Multipole_updateCoords
  end interface updateCoords

  !!* Updates the dipole polarizations for the multipole module.
  interface updateDipoles
    module procedure Multipole_updateDipoles
  end interface updateDipoles

  !!* Returns the multipole shift vectors per atom.
  interface getShifts
    module procedure Multipole_getShifts
  end interface getShifts

  !!* Returns the multipole energy per atom.
  interface getEnergyPerAtom
    module procedure Multipole_getEnergyPerAtom
  end interface getEnergyPerAtom

  !!* Returns the total multipole energy.
  interface getEnergy
    module procedure Multipole_getTotalEnergy
  end interface getEnergy

  !!* Returns the multipole gradient per atom.
  interface addGradientDC
    module procedure Multipole_addGradientDC
  end interface

  !!* Returns the multipole gradient per atom.
  interface addNumGradientDC
    module procedure Multipole_addNumGradientDC
  end interface


contains



  !!* Initializes instance.
  !!* @param inp Input file object
  subroutine Multipole_init(self, inp)
    type(OMultipole), intent(inout) :: self
    type(TMPoleInit), intent(in) :: inp

    integer :: iSp1, iSp2

    if (inp%tPeriodic) then
      !! TODO: Implement support for periodic systems, problem:
      !! Ewald summation of dipole moments//multipole-gamma contributions
      write (*,*) "Error: Multipoles in periodic systems is not implemented", &
          &" yet. Aborting..."
      stop
    elseif (maxval(inp%orb%angShell) > 1) then
      !! TODO: Extension to arbitrary number of orbitals
      write(*,*) "Error: Multipoles only possible up to p-orbitals (angular", &
          &" momentums up to 1). Aborting..."
      stop
    elseif (inp%tOrbResolved) then
      !! Orbital resolved DFTB: Not sure if we need to extend the formalism here!
      write (*,*) "Error: Orbital resolved calculations not possible", &
          &" with the multipole expansion. Aborting..."
      stop
    end if

    self%nAtom = size(inp%orb%nOrbAtom)
    self%nSpecies = size(inp%orb%nOrbSpecies)
    self%mShell = inp%orb%mShell

    @:ASSERT(.not. self%tInitialised)
    @:ASSERT(self%nSpecies == size(inp%hubbU, dim=2))
    @:ASSERT(self%mShell == size(inp%hubbU, dim=1))

    allocate(self%UU(self%nSpecies))
    self%UU(:) = inp%hubbU(1,:)

    !! Get cutoff
    !! @note: normal SCC cutoff for now, should be changed later if the
    !! multipole elements turn out to be longer reaching
    allocate(self%cutoffs(self%nSpecies, self%nSpecies))
    do iSp1 = 1, self%nSpecies
      do iSp2 = iSp1, self%nSpecies
        self%cutoffs(iSp2, iSp1) = &
            & expGammaCutoff(self%UU(iSp2), self%UU(iSp1))
        self%cutoffs(iSp1, iSp2) = self%cutoffs(iSp2, iSp1)
      end do
    end do
    self%mCutoff = maxval(self%cutoffs)

    allocate(self%nNeigh(self%nSpecies, self%nAtom))
    allocate(self%nNeighMax(self%nAtom))
    allocate(self%nNeighShort(self%nSpecies, self%nAtom))

    allocate(self%invRMat01(self%nAtom, self%nAtom, 3))
    allocate(self%invRMat10(self%nAtom, self%nAtom, 3))
    allocate(self%invRMat20(self%nAtom, self%nAtom, 3, 3))
    allocate(self%invRMat11(self%nAtom, self%nAtom, 3, 3))
    allocate(self%shortGamma01(self%nAtom, self%nAtom, 3))
    allocate(self%shortGamma10(self%nAtom, self%nAtom, 3))
    allocate(self%shortGamma20(self%nAtom, self%nAtom, 3, 3))
    allocate(self%shortGamma11(self%nAtom, self%nAtom, 3, 3))
    allocate(self%gamma01(self%nAtom, self%nAtom, 3))
    allocate(self%gamma10(self%nAtom, self%nAtom, 3))
    allocate(self%gamma20(self%nAtom, self%nAtom, 3, 3))
    allocate(self%gamma11(self%nAtom, self%nAtom, 3, 3))

    allocate(self%shift1(self%nAtom))
    allocate(self%shift2(self%nAtom, 3))
    allocate(self%charges(self%nAtom))

    allocate(self%deltaDAtom(self%nAtom, 3))

    self%tInitialised = .True.

  end subroutine Multipole_init



  !!* Returns real space cutoff.
  !!* @param self Instance.
  !!* @return Cutoff.
  function Multipole_getCutoff(self) result(cutoff)
    type(OMultipole), intent(inout) :: self
    real(dp) :: cutoff

    cutoff = self%mCutoff

  end function Multipole_getCutoff



  !!* Updates the number of neighbors for the multipole module (local).
  !!* @param self Instance.
  !!* @param specie Species for each atom
  !!* @param neighList Neighbor list for the atoms in the system.
  subroutine Multipole_updateNNeigh(self, species, neighList)
    type(OMultipole), intent(inout) :: self
    integer, intent(in) :: species(:)
    type(TNeighbourList), intent(in) :: neighList

    integer :: iAt1, iSp2

    self%nNeighShort(:,:) = 0 ! nSpecies, nAtom
    do iAt1 = 1, self%nAtom
      do iSp2 = 1, self%nSpecies
        self%nNeighShort(iSp2, iAt1) = &
            & getNrOfNeighbours(neighList, &
            & self%cutoffs(iSp2, species(iAt1)), iAt1)
      end do
    end do

  end subroutine Multipole_updateNNeigh



  !!* Updates data structures if there are changed coordinates for the instance.
  !!* @param self Instance.
  !!* @param neighList Neighbor list.
  !!* @param species Species for all atoms (nAllAtom).
  subroutine Multipole_updateCoords(self, coord, neighList, species)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(in) :: coord(:,:)
    type(TNeighbourList), intent(in) :: neighList
    integer, intent(in) :: species(:)

    @:ASSERT(self%tInitialised)

    self%invRMat01(:,:,:) = 0.0_dp
    self%invRMat10(:,:,:) = 0.0_dp
    self%invRMat20(:,:,:,:) = 0.0_dp
    self%invRMat11(:,:,:,:) = 0.0_dp

    call Multipole_updateNNeigh(self, species, neighList)
    self%nNeighMax(:) = maxval(self%nNeighShort, dim=1)

    !! First and second derivatives of the coulomb potential
    call invRPDipole(self%invRMat10, self%invRMat01, self%nAtom, coord)
    call invRPPDipole(self%invRMat20, self%invRMat11, self%nAtom, coord)

    !! Build gamma matrices
    call Multipole_initGammaP(self, coord, species, neighList%iNeighbour)
    call Multipole_initGammaPP(self, coord, species, neighList%iNeighbour)

    self%tCoordUp = .true.
    self%tChargeUp = .false.
    self%tDipoleUp = .false.

  end subroutine Multipole_updateCoords



  !!* Updates monopole charges of the multipole module
  !!* @param self Instance.
  !!* @param charges New charges
  subroutine Multipole_updateCharges(self, qq, q0)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(in) :: qq(:), q0(:)

    @:ASSERT(self%tInitialised)

    self%charges(:) = qq(:) - q0(:)

    self%tChargeUp = .true.

  end subroutine Multipole_updateCharges



  !!* Updates the polarizations of the dipoles.
  !!* @param self Multipole instance
  !!* @param iNeighbour Neighbor atom indices
  !!* @param dd Atomic dipole polarizations
  !!* @param d0 Reference dipole polarization (neutral atoms)
  subroutine Multipole_updateDipoles(self, iNeighbour, dd, d0)
    type(OMultipole), intent(inout) :: self
    integer, intent(in) :: iNeighbour(0:,:)
    real(dp), intent(in) :: dd(:,:)
    real(dp), intent(in), optional :: d0(:,:)

    @:ASSERT(self%tInitialised)
    @:ASSERT(size(self%deltaDAtom, dim=2) == 3)
    @:ASSERT(size(self%deltaDAtom, dim=1) == self%nAtom)
    @:ASSERT(size(dd, dim=1) == size(self%deltaDAtom, dim=1))
    @:ASSERT(size(dd, dim=2) == size(self%deltaDAtom, dim=2))
    if (present(d0)) then
      @:ASSERT(size(d0, dim=1) == size(self%deltaDAtom, dim=1))
      @:ASSERT(size(d0, dim=2) == size(self%deltaDAtom, dim=2))
    end if

    !! Build dipole fluctuations
    if (present(d0)) then
      self%deltaDAtom(:,:) = dd(:,:) - d0(:,:)
    else
      self%deltaDAtom(:,:) = dd(:,:)
    end if

    !! Update atomic dipole shift vectors
    call Multipole_buildShifts(self, iNeighbour)

    self%tDipoleUp = .true.

  end subroutine Multipole_updateDipoles



  !!* Set up the first derivative of gamma.
  !!* @param self Multipole instance
  !!* @param coord List of coordinates
  !!* @param species List of the species for each atom.
  !!* @param iNeighbour Index of neighboring atoms for each atom.
  subroutine Multipole_initGammaP(self, coord, species, iNeighbour)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(in) :: coord(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbour(0:,:)

    integer :: iAt1, iAt2, iNeigh, iSp1, iSp2
    real(dp) :: r1(3) ! Vector to atom 1
    real(dp) :: r2(3) ! Vector to atom 2
    real(dp) :: r12 ! Distance between atom 1 and 2
    real(dp) :: u1, u2 ! Hubbard parameters for atom 1 and 2
    real(dp) :: tmpExpGammaP

    !! Reallocate shortGammaP if it is unable to hold enough neighbors
    if (size(self%shortGamma01, dim=2) < maxval(self%nNeighShort) + 1 .or. &
        &size(self%shortGamma10, dim=2) < maxval(self%nNeighShort) + 1) then
      deallocate(self%shortGamma01)
      deallocate(self%shortGamma10)
      allocate(self%shortGamma01&
          &(0:maxval(self%nNeighShort), self%nAtom, 3))
      allocate(self%shortGamma10&
          &(0:maxval(self%nNeighShort), self%nAtom, 3))
    end if
    self%shortGamma01(:,:,:) = 0.0_dp
    self%shortGamma10(:,:,:) = 0.0_dp

    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      u1 = self%UU(iSp1)
      do iNeigh = 0, maxval(self%nNeighShort(:,iAt1))
        iAt2 = iNeighbour(iNeigh,iAt1)
        iSp2 = species(iAt2)
        u2 = self%UU(iSp2)
        !! Transform coordinates into (-1,0,1) -> y,z,x system
        r1(1) = coord(2,iAt1)
        r1(2) = coord(3,iAt1)
        r1(3) = coord(1,iAt1)
        r2(1) = coord(2,iAt2)
        r2(2) = coord(3,iAt2)
        r2(3) = coord(1,iAt2)
        r12 = sqrt(sum((r1(:)-r2(:))**2))
        if (r12 < tiny(0.0_dp)) then
          self%shortGamma10(iAt1,iAt2,:) = 0.0_dp
          self%shortGamma10(iAt2,iAt1,:) = 0.0_dp
          self%shortGamma01(iAt1,iAt2,:) = 0.0_dp
          self%shortGamma01(iAt2,iAt1,:) = 0.0_dp
        elseif (iNeigh <= self%nNeighShort(iSp2,iAt1)) then
          tmpExpGammaP = expGammaPrime(r12, u1, u2)
          self%shortGamma10(iAt1,iAt2,:) = &
              &tmpExpGammaP * (r1(:) - r2(:)) / r12
          self%shortGamma01(iAt1,iAt2,:) = &
              &tmpExpGammaP * (r2(:) - r1(:)) / r12
          !! Add (atom2, atom1) matrix elements
          if (iAt1 /= iAt2) then
            !! Sign change for tmpExpGammaP because the hubbard parameters
            !! are exchanged when atom indices are exchanged.
            self%shortGamma10(iAt2,iAt1,:) = &
                &-tmpExpGammaP * (r1(:) - r2(:)) / r12
            self%shortGamma01(iAt2,iAt1,:) = &
                &-tmpExpGammaP * (r2(:) - r1(:)) / r12
          end if
        end if
      end do
    end do

    !! Build gammaP from invRMatP and shortGammaP
    !! Minus because gammaP = invRMatP - expGammaP
    self%gamma10(:,:,:) = self%invRMat10(:,:,:) - self%shortGamma10(:,:,:)
    self%gamma01(:,:,:) = self%invRMat01(:,:,:) - self%shortGamma01(:,:,:)

  end subroutine Multipole_initGammaP



  !!* Set up the second derivative of gamma.
  !!* @param self Multipole instance
  !!* @param coord List of coordinates
  !!* @param species List of the species for each atom.
  !!* @param iNeighbour Index of neighboring atoms for each atom.
  !!* @note All limiting cases (r -> 0 and U1 -> U2) are covered in the
  !!* expGammaPrime2 routine in short_gamma.F90, all invRMat and extraTerm
  !!* values for these cases are set to 0.0_dp or produce errors at run-time.
  subroutine Multipole_initGammaPP(self, coord, species, iNeighbour)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(in), target :: coord(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbour(0:,:)

    integer :: iAt1, iAt2, iNeigh, iSp1, iSp2, iMM
    real(dp) :: r1(3) ! Vector to atom 1
    real(dp) :: r2(3) ! Vector to atom 2
    real(dp) :: r12 ! Distance between atom 1 and 2
    real(dp) :: u1, u2 ! Hubbard parameters for atom 1 and 2
    real(dp) :: tmpExpGammaPP

    !! Reallocate shortGammaPP, if it does not contain enough neighbors
    if (size(self%shortGamma20, dim=2) < self%nAtom .or. &
        &size(self%shortGamma11, dim=2) < self%nAtom) then
      deallocate(self%shortGamma20)
      deallocate(self%shortGamma11)
      allocate(self%shortGamma20 &
          &(self%nAtom, self%nAtom, 3, 3))
      allocate(self%shortGamma11 &
          &(self%nAtom, self%nAtom, 3, 3))
    end if
    self%shortGamma20(:,:,:,:) = 0.0_dp
    self%shortGamma11(:,:,:,:) = 0.0_dp

    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      u1 = self%UU(iSp1)
      do iNeigh = 0, maxval(self%nNeighShort(:,iAt1))
        iAt2 = iNeighbour(iNeigh,iAt1)
        iSp2 = species(iAt2)
        u2 = self%UU(iSp2)
        r1(1) = coord(2,iAt1)
        r1(2) = coord(3,iAt1)
        r1(3) = coord(1,iAt1)
        r2(1) = coord(2,iAt2)
        r2(2) = coord(3,iAt2)
        r2(3) = coord(1,iAt2)
        r12 = sqrt(sum((r1(:)-r2(:))**2))
        if (iNeigh <= self%nNeighShort(iSp2,iAt1)) then
          tmpExpGammaPP = expGammaPrime2(r12, u1, u2)
          do iMM = 1, 3
            self%shortGamma20(iAt2,iAt1,iMM,:) = &
                &tmpExpGammaPP * (r1(iMM) - r2(iMM)) * (r1(:) - r2(:))
            !! Add extra term for second derivative where both derivations are
            !! wrt the same coordinate (x, y, or z)
            self%shortGamma20(iAt2,iAt1,iMM,iMM) = &
                &self%shortGamma20(iAt2,iAt1,iMM,iMM) + &
                &extraTermGammaPrime2(r12, u1, u2)
            !! gamma11 = -gamma20
            self%shortGamma11(iAt2,iAt1,iMM,:) = &
                &-self%shortGamma20(iAt2,iAt1,iMM,:)
          end do
          if (iAt1 /= iAt2) then
            !! Complete the skew matrix (might be removed later)
            self%shortGamma11(iAt1,iAt2,:,:) = &
                &self%shortGamma11(iAt2,iAt1,:,:)
            self%shortGamma20(iAt1,iAt2,:,:) = &
                &self%shortGamma20(iAt2,iAt1,:,:)
          end if
        end if
      end do
    end do

    !! Build gammaP from invRMatP and shortGammaP
    !! Minus because gammaP = invRMatP - expGammaP
    self%gamma20(:,:,:,:) = self%invRMat20(:,:,:,:) - self%shortGamma20(:,:,:,:)
    self%gamma11(:,:,:,:) = self%invRMat11(:,:,:,:) - self%shortGamma11(:,:,:,:)

  end subroutine Multipole_initGammaPP



  !!* Add gradient component resulting from the derivative of the potential.
  !!* @param self Instance.
  !!* @param iNeighbour Neighbor list of atoms
  !!* @param species Specie for each atom.
  !!* @param derivs Gradient on exit.
  subroutine Multipole_addGradientDC(self, iNeighbour, species, derivs)
    type(OMultipole), intent(inout) :: self
    integer, intent(in) :: iNeighbour(0:,:)
    integer, intent(in) :: species(:)
    real(dp), intent(inout) :: derivs(:,:)

    integer :: iAt1, iAt2, iSp1, iSp2, iCC, ii!, iNeigh
    real(dp) :: tmpR1, tmpR2, tmpDerivs(3,self%nAtom)

    tmpDerivs(:,:) = 0.0_dp
    do iCC = 1, 3
      if (iCC == 1) then
        ii = 3
      else if (iCC == 2) then
        ii = 1
      else if (iCC == 3) then
        ii = 2
      end if
      do iAt1 = 1, self%nAtom
        iSp1 = species(iAt1)
        do iAt2 = 1, self%nAtom
          iSp2 = species(iAt2)
          !! Terms including the second derivatives of gamma.
          tmpR1 = sum(self%gamma20(iAt1,iAt2,:,ii) * &
              &(self%deltaDAtom(iAt1,:) * self%charges(iAt2) - &
              &self%deltaDAtom(iAt2,:) * self%charges(iAt1)))
          tmpDerivs(iCC,iAt2) = tmpDerivs(iCC,iAt2) + tmpR1
          ! tmpDerivs(iCC,iAt1) = tmpDerivs(iCC,iAt1) - tmpR1
        end do
      end do
    end do

    ! derivs(:,:) = derivs(:,:) + tmpDerivs(:,:)

    print *, "MD derivative gamma (analytic):"
    do iAt1 = 1, self%nAtom
      print *, tmpDerivs(:,iAt1)
    end do

  end subroutine Multipole_addGradientDC



  !! TEMP NUMERIC DC GRADIENT
  subroutine Multipole_addNumGradientDC(self, coord, species, neighList, derivs)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(in) :: coord(:,:)
    integer, intent(in) :: species(:)
    type(TNeighbourList), intent(in) :: neighList
    real(dp), intent(inout) :: derivs(:,:)

    integer :: iAt, iCC, iLR
    real(dp) :: geom(3,self%nAtom), energies(2), tmpDerivs(3,self%nAtom)
    real(dp), parameter :: deltaXDiff = epsilon(1.0_dp)**0.25_dp

    @:ASSERT(self%tInitialised)

    self%invRMat01(:,:,:) = 0.0_dp
    self%invRMat10(:,:,:) = 0.0_dp

    do iAt = 1, self%nAtom
      do iCC = 1, 3
        do iLR = 1, 2
          geom(:,:) = coord(:,:)
          geom(iCC,iAt) = geom(iCC,iAt) + real(3 * iLR - 3, dp) * deltaXDiff
          call invRPDipole(self%invRMat10, self%invRMat01, self%nAtom, geom)
          call Multipole_initGammaP(self, geom, species, neighList%iNeighbour)
          call Multipole_getTotalEnergy(self, energies(iLR))
        end do
        tmpDerivs(iCC,iAt) = (energies(2) - energies(1)) / (2.0_dp * deltaXDiff)
      end do
    end do

    derivs(:,:) = derivs(:,:) + tmpDerivs(:,:)

    print *, "MD derivative gamma (numeric):"
    do iAt = 1, self%nAtom
      print *, tmpDerivs(:,iAt)
    end do

  end subroutine Multipole_addNumGradientDC



  !!* Constructs dipole polarization shift for Multipole contributions.
  !!* @param self Multipole instance
  !!* @param iNeighbour Neighbor list of atoms
  subroutine Multipole_buildShifts(self, iNeighbour)
    type(OMultipole), target, intent(inout) :: self
    integer, intent(in) :: iNeighbour(0:,:)

    integer :: iAt1, iAt2, iNeigh

    @:ASSERT(self%tInitialised)
    @:ASSERT(self%tChargeUp)
    @:ASSERT(size(self%shift1) == size(self%shift2, dim=1))

    self%shift1(:) = 0.0_dp
    self%shift2(:,:) = 0.0_dp

    !! Adding the short-range part of the shift vectors.
    do iAt1 = 1, self%nAtom
      do iNeigh = 0, self%nNeighMax(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        !! Sign is taken care of in the initGammaP, so just summing up here.
        self%shift1(iAt2) = self%shift1(iAt2) + &
            &dot_product(self%gamma01(iAt2,iAt1,:), &
            &self%deltaDAtom(iAt1,:))
        self%shift2(iAt2,:) = self%shift2(iAt2,:) + self%charges(iAt1) * &
            &self%gamma10(iAt2,iAt1,:)
        if (iAt1 /= iAt2) then
          self%shift1(iAt1) = self%shift1(iAt1) + &
              &dot_product(self%gamma01(iAt1,iAt2,:), &
              &self%deltaDAtom(iAt2,:))
          self%shift2(iAt1,:) = self%shift2(iAt1,:) + self%charges(iAt2) * &
              &self%gamma10(iAt1,iAt2,:)
        end if
      end do
    end do

  end subroutine Multipole_buildShifts



  !!* Returns shifts per atom.
  !!* @param self Instance.
  !!* @param shift1 Shift vector (scalar per atom)
  !!* @param shift2 Spatial shift vector (3 dim per atom)
  subroutine Multipole_getShifts(self, shift1, shift2)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(out) :: shift1(:), shift2(:,:)

    @:ASSERT(size(shift1) == self%nAtom)
    @:ASSERT(all(shape(shift2) == (/self%nAtom, 3/)))

    shift1(:) = self%shift1(:)
    shift2(:,:) = self%shift2(:,:)

  end subroutine Multipole_getShifts



  !!* Returns energy per atom.
  !!* @param self Instance.
  !!* @param energyPerAtom Energy per atom.
  subroutine Multipole_getEnergyPerAtom(self, energyPerAtom)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(out) :: energyPerAtom(:)

    integer :: iAt

    @:ASSERT(size(energyPerAtom) == self%nAtom)

    do iAt = 1, self%nAtom
      energyPerAtom(iAt) = 0.5_dp * (self%shift1(iAt) * self%charges(iAt) + &
          &dot_product(self%shift2(iAt,:), self%deltaDAtom(iAt,:)))
    end do

  end subroutine Multipole_getEnergyPerAtom



  !!* Returns energy per atom.
  !!* @param self Instance.
  !!* @param energy Total multipole energy.
  subroutine Multipole_getTotalEnergy(self, energy)
    type(OMultipole), intent(inout) :: self
    real(dp), intent(out) :: energy

    integer :: iAt1, iAt2

    energy = 0.0_dp

    do iAt1 = 1, self%nAtom
      do iAt2 = 1, self%nAtom
        energy = energy + 0.5_dp * &
            &(dot_product(self%gamma10(iAt1,iAt2,:), &
            &self%deltaDAtom(iAt1,:)) * self%charges(iAt2) + &
            &dot_product(self%gamma01(iAt1,iAt2,:), &
            &self%deltaDAtom(iAt2,:)) * self%charges(iAt1))
      end do
    end do

  end subroutine Multipole_getTotalEnergy



end module Multipole
