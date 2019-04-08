!!* Routines that read first moment files and build the position matrix.
module positionmatrix
#:include 'common.fypp'
  use accuracy
  use commontypes
  use sknew
  use slakocont
  implicit none
  private

  public :: allocatePosMat, buildPosMat, buildPosMatPrime

  integer, parameter :: mAngRot_ = 1 ! Max. angular momentum for rotations
  integer, parameter :: mMoment = 1 ! Highest allowed moment of the overlap

contains

  !!* Allocates the position matrix.
  !!* @param posMtx  Position matrix (sparse array)
  !!* @param iNeighbor  Respective list of the neighbors
  !!* @param nNeighbor  Number of surrounding neighbors for each atom
  !!* @param species  Chemical species in the system
  !!* @param orb  Information about the orbitals for each species
  subroutine allocatePosMat(posMtx, iNeighbor, nNeighbor, species, orb)
    real(dp), pointer, intent(out) :: posMtx(:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb

    integer :: iAt1, iAt2, iNeigh1 ! Atom numbers
    integer :: nAtom ! Total number of atoms in the system
    integer :: nOrb1, nOrb2 ! Number of orbitals on atom 1 and 2
    integer :: nElem, nElemOld ! Number of matrix elements for the sparse posmat

    nAtom = size(iNeighbor, dim=2)

    if (associated(posMtx)) then
      nElemOld = size(posMtx, dim=2)
    else
      nElemOld = 0
    end if

    nElem = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbSpecies(species(iAt1))
      do iNeigh1 = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh1, iAt1)
        nOrb2 = orb%nOrbSpecies(species(iAt2))
        nElem = nElem + nOrb1 * nOrb2
      end do
    end do

    if (nElemOld /= nElem) then
      if (associated(posMtx)) then
        deallocate(posMtx)
      end if
      allocate(posMtx(nElem, 3))
    end if

  end subroutine allocatePosMat



  !!* Main routine that builds the sparse position matrix of the system.
  !!* @param posMtx1  Returns position matrix with dipole origin at 1st center.
  !!* @param posMtx2  Returns position matrix with dipole origin at 2nd center.
  !!* @param skOverCont  Container for SK overlap integrals.
  !!* @param fmOnSites  Table containing first moment on-site integrals
  !!* @param fmCont  Container for the first moment integrals.
  !!* @param coords  List of all atomic coordinates of the system
  !!* @param nNeighbor  Number of surrounding neighbors for each atom
  !!* @param iNeighbor  Respective list of the neighbors
  !!* @param iPair  Points to the starting index for the interactions between
  !!* two atoms in the sparse form.
  !!* @param species  Chemical species in the system
  !!* @param orb  Information about the orbitals for each species
  subroutine buildPosMat(posMtx1, posMtx2, skOverCont, fmOnSites, fmCont, coords, &
      &nNeighbor, iNeighbor, iPair, species, orb)
    real(dp), intent(out) :: posMtx1(:,:)
    real(dp), intent(out) :: posMtx2(:,:)
    type(OSlakoCont), intent(in) :: skOverCont
    real(dp), intent(in) :: fmOnSites(:,:)
    type(OSlakoCont), intent(in) :: fmCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb

    real(dp) :: dirCos(3) ! Directional cosines for one atom pair
    real(dp) :: dist ! Distance between two atoms
    real(dp) :: interSKOver(getMIntegrals(skOverCont)) ! Overlap interactions
    real(dp) :: interFirstMom(getMIntegrals(fmCont)) ! FM interactions
    real(dp) :: fmOrigins(3,getMIntegrals(fmCont)) ! FM Origin vectors
    integer :: fmCenter(getMIntegrals(fmCont)) ! FM dipole centers

    real(dp) :: XX1(orb%mOrb,orb%mOrb,3) ! Diatomic first moment block
    real(dp) :: XX2(orb%mOrb,orb%mOrb,3) ! Diatomic first moment block
    integer :: ind ! Position index in the sparse position matrix
    integer :: nAtom ! Total number of atoms
    integer :: nOrb1, nOrb2 ! Total number of orbitals on atom 1 and 2
    integer :: iAt1, iAt2, iSp1, iSp2, iNeigh1 ! Indexers
    integer :: iMM ! Indexer for the third dimension of the position matrix

    nAtom = size(nNeighbor)
    posMtx1(:,:) = 0.0_dp
    posMtx2(:,:) = 0.0_dp
    fmOrigins(:,:) = 0.0_dp
    fmCenter(:) = 0

    !! Loops to build on-site blocks.
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      ind = iPair(0,iAt1)
      dist = sqrt(sum(coords(:,iAt1)**2))
      call buildOnSiteBlock(XX1, fmOnSites, fmCont, iSp1, orb)
      do iMM = 1, 3
        posMtx1(ind+1:ind+nOrb1**2,iMM) = reshape(XX1(1:nOrb1,1:nOrb1,iMM), &
            &(/ nOrb1**2 /))
        posMtx2(ind+1:ind+nOrb1**2,iMM) = reshape(XX1(1:nOrb1,1:nOrb1,iMM), &
            &(/ nOrb1**2 /))
      end do
    end do

    !! Loops to build diatomic blocks.
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh1 = 1, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh1, iAt1)
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecies(iSp2)
        ind = iPair(iNeigh1,iAt1)
        dirCos(:) = coords(:,iAt2) - coords(:,iAt1)
        dist = sqrt(sum(dirCos(:)**2))
        dirCos(:) = dirCos(:) / dist
        call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
        call getSKIntegrals(fmCont, interFirstMom, dist, iSp1, iSp2)
        call getMPOrigins(fmCont, fmOrigins, fmCenter, dist, iSp1, iSp2)
        call buildOffSiteBlock(XX1, XX2, interSKover, interFirstMom, fmOrigins, &
            &fmCenter, iSp1, iSp2, dirCos, dist, orb)
        do iMM = 1, 3
          posMtx1(ind+1:ind+nOrb2*nOrb1,iMM) = reshape(XX1(1:nOrb2,1:nOrb1,iMM), &
              &(/ nOrb2 * nOrb1 /))
          posMtx2(ind+1:ind+nOrb2*nOrb1,iMM) = reshape(XX2(1:nOrb2,1:nOrb1,iMM), &
              &(/ nOrb2 * nOrb1 /))
        end do
      end do
    end do

    ! print *
    ! print *, "Position Matrix 1, y-Dipole"
    ! print *, posMtx1(:,1)
    ! print *
    ! print *, "Position Matrix 1, z-Dipole"
    ! print *, posMtx1(:,2)
    ! print *
    ! print *, "Position Matrix 1, x-Dipole"
    ! print *, posMtx1(:,3)

  end subroutine buildPosMat



  !!* Main routine that builds the diatomic block of the sparse first moment
  !!* tensor derivative.
  !!* @param iAt1 Atom 1 of the block
  !!* @param iAt2 Atom 2 of the block
  !!* @param pos1Prime  Returns pos. matrix with dipole origin at 1st center.
  !!* @param pos2Prime  Returns pos. matrix with dipole origin at 2nd center.
  !!* @param SS  Sparse overlap matrix
  !!* @param skOverCont  Container for SK overlap integrals.
  !!* @param fmCont  Container for the first moment integrals.
  !!* @param coords  List of all atomic coordinates of the system
  !!* @param nNeighbor  Number of surrounding neighbors for each atom
  !!* @param iNeighbor  Respective list of the neighbors
  !!* @param iPair  Points to the starting index for the interactions between
  !!* two atoms in the sparse form.
  !!* @param species  Chemical species in the system
  !!* @param orb  Information about the orbitals for each species
  subroutine buildPosMatPrime(iAt1, iAt2, pos1Prime, pos2Prime, skOverCont, &
      &fmCont, coords, nNeighbor, iNeighbor, iPair, species, orb)
    integer, intent(in) :: iAt1, iAt2
    real(dp), intent(out) :: pos1Prime(:,:,:)
    real(dp), intent(out) :: pos2Prime(:,:,:)
    type(OSlakoCont), intent(in) :: skOverCont
    type(OSlakoCont), intent(in) :: fmCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb

    real(dp) :: dirCos(3) ! Directional cosines for one atom pair
    real(dp) :: dist ! Distance between two atoms
    real(dp) :: interSKOver(getMIntegrals(skOverCont)) ! Overlap interactions
    real(dp) :: interFirstMom(getMIntegrals(fmCont)) ! FM interactions
    real(dp) :: fmOrigins(3,getMIntegrals(fmCont)) ! FM Origin vectors
    integer :: fmCenter(getMIntegrals(fmCont)) ! FM dipole centers

    real(dp) :: XX1(orb%mOrb,orb%mOrb,3) ! Diatomic first moment block
    real(dp) :: XX2(orb%mOrb,orb%mOrb,3) ! Diatomic first moment block
    real(dp) :: pos1Tmp(size(pos1Prime, dim=1),3,2,3)
    real(dp) :: pos2Tmp(size(pos2Prime, dim=1),3,2,3)
    integer :: ind ! Position index in the sparse position matrix
    integer :: nAtom ! Total number of atoms
    integer :: nOrb1, nOrb2 ! Total number of orbitals on atom 1 and 2
    integer :: iSp1, iSp2, iNeigh1 ! Indexers
    integer :: iCC ! Indexer for x,y,z-coordinates
    integer :: iMM ! Indexer for the third dimension of the position matrix
    integer :: iLR ! Loop for left and right derivatives

    real(dp), parameter :: deltaXDiff = epsilon(1.0_dp)**0.25_dp

    @:ASSERT(size(pos1Prime, dim=2) == 3)
    @:ASSERT(size(pos1Prime, dim=3) == 3)
    @:ASSERT(all(shape(pos1Prime) == shape(pos2Prime)))

    nAtom = size(nNeighbor)
    pos1Prime(:,:,:) = 0.0_dp
    pos2Prime(:,:,:) = 0.0_dp
    pos1Tmp(:,:,:,:) = 0.0_dp
    pos2Tmp(:,:,:,:) = 0.0_dp
    fmOrigins(:,:) = 0.0_dp
    fmCenter(:) = 0

    !! Build a diatomic block of the finite difference position matrices
    do iCC = 1, 3
      do iLR = 1, 2 ! 2 * iLR - 3 is either -1 or 1
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecies(iSp2)
        !! Calculate normalized directional cosine
        dirCos(:) = coords(:,iAt2) - coords(:,iAt1)
        dirCos(:) = shiftCoordinate(dirCos(:), iCC, iLR)
        dist = sqrt(sum(dirCos(:)**2))
        dirCos(:) = dirCos(:) / dist
        !! Gather necessary integral values
        call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
        call getSKIntegrals(fmCont, interFirstMom, dist, iSp1, iSp2)
        call getMPOrigins(fmCont, fmOrigins, fmCenter, dist, iSp1, iSp2)
        call buildOffSiteBlock(XX1, XX2, interSKover, interFirstMom, &
            &fmOrigins, fmCenter, iSp1, iSp2, dirCos, dist, orb)
        !! Calculate Position Matrices
        do iMM = 1, 3
          pos1Tmp(1:nOrb2*nOrb1,iMM,iLR,iCC) = reshape( &
              &XX1(1:nOrb2,1:nOrb1,iMM), (/ nOrb2 * nOrb1 /))
          pos2Tmp(1:nOrb2*nOrb1,iMM,iLR,iCC) = reshape( &
              &XX2(1:nOrb2,1:nOrb1,iMM), (/ nOrb2 * nOrb1 /))
        end do
      end do
    end do

    !! Calculate derivative by finite differences
    do iCC = 1, 3
      pos1Prime(:,:,iCC) = (pos1Tmp(:,:,2,iCC) - pos1Tmp(:,:,1,iCC)) / &
          &(2.0_dp * deltaXDiff)
      pos2Prime(:,:,iCC) = (pos2Tmp(:,:,2,iCC) - pos2Tmp(:,:,1,iCC)) / &
          &(2.0_dp * deltaXDiff)
    end do

  end subroutine buildPosMatPrime



  !!* Builds the on-site blocks of the position matrix.
  !!* @param XX  rectangular diatomic block matrix (dipole origin on center 1)
  !!* @param fmOnSites  Table containing first moment on-site integrals
  !!* @param fmCont  Container for the first moment integrals.
  !!* @param iSp  Chemical species integer for atom
  !!* @param orb  Information about the orbitals for each species
  subroutine buildOnSiteBlock(XX, fmOnSites, fmCont, iSp, orb)
    real(dp), intent(out) :: XX(:,:,:)
    real(dp), intent(in) :: fmOnSites(:,:)
    type(OSlakoCont), intent(in) :: fmCont
    integer, intent(in) :: iSp
    type(TOrbitals), intent(in) :: orb

    integer :: l1, l2, iSh1, iSh2, nOrb1, nOrb2
    integer :: iFM, iCol, iRow, iMM
    real(dp), target :: mat(mAngRot_*2+1,mAngRot_*2+1,-mMoment:mMoment)
    real(dp), pointer :: pMat(:,:,:) ! Pointer to select snippets from mat
    real(dp), pointer :: pOver(:) ! Pointer to needed overlap integrals
    real(dp), pointer :: pFirstMom(:) ! Pointer to needed first moments
    integer, pointer :: pCenter(:) ! Pointer to needed origin vectors
    !! The size of the over-array is more or less arbitrary, it just needs to be
    !! of size >= size(interSKOver) and passing skOverCont for that is overkill.
    real(dp), target :: over(getMIntegrals(fmCont))
    integer, target :: centers(getMIntegrals(fmCont))
    real(dp), target :: firstmom(getMIntegrals(fmCont))
    integer :: nIntegs(0:mAngRot_,0:mAngRot_) = reshape((/1, 2, 2, 4/), &
        &(/mAngRot_ + 1, mAngRot_ + 1/))

    @:ASSERT(maxval(orb%angShell) <=  mAngRot_)
    @:ASSERT(all(shape(XX) >= (/ orb%nOrbSpecies(iSp), orb%nOrbSpecies(iSp), 3 /)))

    XX(:,:,:) = 0.0_dp
    over(:) = 0.0_dp
    centers(:) = 0
    iFM = 1
    iCol = 1
    call prepareOnSites(firstmom, fmOnSites, iSp, orb, nIntegs)
    do iSh1 = 1, orb%nShell(iSp)
      l1 = orb%angShell(iSh1,iSp)
      nOrb1 = 2 * l1 + 1
      iRow = 1
      do iSh2 = 1, orb%nShell(iSp)
        l2 = orb%angShell(iSh2,iSp)
        nOrb2 = 2 * l2 + 1
        @:ASSERT(size(over) >= 1 + min(l1, l2))
        @:ASSERT(size(firstmom) >= iFM + nIntegs(l1,l2) - 1)
        pOver => over(1:1+min(l1,l2))
        pFirstMom => firstmom(iFM:iFM+nIntegs(l1,l2)-1)
        pCenter => centers(iFM:iFM+nIntegs(l1,l2)-1)
        if (l2 <= l1) then
          pMat => mat(1:2*l1+1,1:2*l2+1,:)
        else
          pMat => mat(1:2*l2+1,1:2*l1+1,:)
        end if
        call fillFMBlock(pMat, l1, l2, pFirstMom, pCenter, pOver, 0.0_dp, 1)
        do iMM = 1, 3
          if (l1 <= l2) then
            @:ASSERT(all(shape(XX(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM)) == &
                &shape(pMat(:,:,iMM))))
            XX(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM) = pMat(:,:,iMM)
          else
            @:ASSERT(all(shape(XX(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM)) == &
                &shape(transpose(pMat(:,:,iMM)))))
            XX(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM) = transpose(pMat(:,:,iMM))
          end if
        end do
        iFM = iFM + nIntegs(l1,l2)
        iRow = iRow + nOrb2
      end do
      iCol = iCol + nOrb1
    end do

  end subroutine buildOnSiteBlock



  !!* Builds a diatomic overlap block from reference integrals.n
  !!* @param XX1  rectangular diatomic block matrix (dipole origin on center 1)
  !!* @param XX2  rectangular diatomic block matrix (dipole origin on center 2)
  !!* @param interOver  Overlap table for the iSp1-iSp2 interaction
  !!* @param interFirstMom  First Moment table for the iSp1-iSp2 interaction
  !!* @param fmOrigins  Table containing the origin vectors for all interactions
  !!* @param fmCenter  Centeor of dipole origin
  !!* @param iSp1  Chemical species integer for atom 1
  !!* @param iSp2  Chemical species integer for atom 2
  !!* @param dirCos  Directional cosine vector (for transformation)
  !!* @param dist  Distance between the centers
  !!* @param orb  Orbital information of chemical species in the system
  !!* @caveat Only angular momenta up to "p" are allowed
  subroutine buildOffSiteBlock(XX1, XX2, interSKOver, interFirstMom, fmOrigins, &
      &fmCenter, iSp1, iSp2, dirCos, dist, orb)
    real(dp), intent(out) :: XX1(:,:,:)
    real(dp), intent(out) :: XX2(:,:,:)
    real(dp), intent(in), target :: interSKOver(:)
    real(dp), intent(in), target :: interFirstMom(:)
    real(dp), intent(in), target :: fmOrigins(:,:)
    integer, intent(in), target :: fmCenter(:)
    integer, intent(in) :: iSp1, iSp2
    real(dp), intent(in) :: dirCos(3)
    real(dp), intent(in) :: dist
    type(TOrbitals), intent(in) :: orb

    integer :: iCol, iRow ! Indexer for positions in the output overlap matrix
    integer :: iOver, iFM ! Position indices for integral lists
    integer :: iSh1, iSh2 ! Loop Indexer
    integer :: l1, l2, nOrb1, nOrb2
    integer :: iMM ! Indexer for the third dimension of the position matrix
    !! Wether a sign change is needed for exchange of indices in X2
    real(dp), target :: mat1(mAngRot_*2+1,mAngRot_*2+1,-mMoment:mMoment)
    real(dp), target :: mat2(mAngRot_*2+1,mAngRot_*2+1,-mMoment:mMoment)
    real(dp), pointer :: pMat1(:,:,:) ! Pointer to select snippets from mat
    real(dp), pointer :: pMat2(:,:,:) ! Pointer to select snippets from mat
    real(dp), pointer :: pOver(:) ! Pointer to needed overlap integrals
    real(dp), pointer :: pFirstMom(:) ! Pointer to needed first moments
    real(dp), pointer :: pOrigins(:,:) ! Pointer to needed origin vectors
    integer, pointer :: pCenter(:) ! Pointer to needed center of dipole origin
    integer :: nIntegs(0:mAngRot_,0:mAngRot_) = reshape((/1, 2, 2, 4/), &
        &(/mAngRot_ + 1, mAngRot_ + 1/))

    !! Safety stop while in development
    if (max(maxval(orb%angShell(:,iSp1)), maxval(orb%angShell(:,iSp2))) &
        > 1) then
      write (*,*) "No angular momenta larger than 1 (p-Orbitals) supported in &
          &multipole expansion."
      STOP
    end if

    @:ASSERT(maxval(orb%angShell) <=  mAngRot_)
    @:ASSERT(all(shape(XX1) >= (/ orb%nOrbSpecies(iSp1), &
        &orb%nOrbSpecies(iSp2), 3 /)))
    @:ASSERT(all(shape(XX2) >= (/ orb%nOrbSpecies(iSp1), &
        &orb%nOrbSpecies(iSp2), 3 /)))

    XX1(:,:,:) = 0.0_dp
    XX2(:,:,:) = 0.0_dp
    iOver = 1
    iFM = 1
    iCol = 1
    !! Build the first and second center diatomic block via shell-pair matrices.
    do iSh1 = 1, orb%nShell(iSp1)
      l1 = orb%angShell(iSh1,iSp1)
      nOrb1 = 2 * l1 + 1
      iRow = 1
      do iSh2 = 1, orb%nShell(iSp2)
        l2 = orb%angShell(iSh2,iSp2)
        nOrb2 = 2 * l2 + 1
        !! Check if the size of the first moment integrals array is large enough
        @:ASSERT(size(interSKOver) >= iOver + min(l1, l2))
        @:ASSERT(size(interFirstMom) >= iFM + nIntegs(l1,l2) - 1)
        pOver => interSKOver(iOver:iOver+min(l1,l2))
        pFirstMom => interFirstMom(iFM:iFM+nIntegs(l1,l2)-1)
        pOrigins => fmOrigins(:,iFM:iFM+nIntegs(l1,l2)-1)
        pCenter => fmCenter(iFM:iFM+nIntegs(l1,l2)-1)
        !! Selects orbital block size to fill
        if (l2 <= l1) then
          pMat1 => mat1(1:2*l1+1,1:2*l2+1,:)
          pMat2 => mat2(1:2*l1+1,1:2*l2+1,:)
        else
          pMat1 => mat1(1:2*l2+1,1:2*l1+1,:)
          pMat2 => mat2(1:2*l2+1,1:2*l1+1,:)
        end if
        call fillFMBlock(pMat1, l1, l2, pFirstMom, pCenter, pOver, dist, 1)
        call fillFMBlock(pMat2, l1, l2, pFirstMom, pCenter, pOver, dist, 2)
        call rotateFM(pMat1, dirCos(1), dirCos(2), dirCos(3))
        call rotateFM(pMat2, dirCos(1), dirCos(2), dirCos(3))
        do iMM = 1, 3
          !! For each dipole fill the diatomic block with integrals.
          if (l1 <= l2) then
            @:ASSERT(all(shape(XX1(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM)) == &
                &shape(pMat1(:,:,iMM))))
            @:ASSERT(all(shape(XX2(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM)) == &
                &shape(pMat2(:,:,iMM))))
            XX1(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM) = pMat1(:,:,iMM)
            XX2(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM) = pMat2(:,:,iMM)
          else
            @:ASSERT(all(shape(XX1(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM)) == &
                &shape(transpose(pMat1(:,:,iMM)))))
            @:ASSERT(all(shape(XX2(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM)) == &
                &shape(transpose(pMat2(:,:,iMM)))))
            XX1(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM) = &
                &transpose(pMat1(:,:,iMM))
            XX2(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1,iMM) = &
                &transpose(pMat2(:,:,iMM))
          end if
        end do
        iFM = iFM + nIntegs(l1,l2)
        iOver = iOver + min(l1, l2) + 1
        iRow = iRow + nOrb2
      end do
      iCol = iCol + nOrb1
    end do

    ! print *
    ! print *, "Position Matrix 1, y-Dipole"
    ! do iMM = 1, size(XX1, dim=2)
    !   print *, XX1(:,iMM,1)
    ! end do
    ! print *
    ! print *, "Position Matrix 1, z-Dipole"
    ! do iMM = 1, size(XX1, dim=2)
    !   print *, XX1(:,iMM,2)
    ! end do
    ! print *
    ! print *, "Position Matrix 1, x-Dipole"
    ! do iMM = 1, size(XX1, dim=2)
    !   print *, XX1(:,iMM,3)
    ! end do

  end subroutine buildOffSiteBlock



  !!* Prepare on-site integral array so it can be used by the block-building
  !!* routines.
  !!* @param firstmom  Contains the prepared array on exit
  !!* @param fmOnSites  Table containing first moment on-site integrals
  !!* @param iSp  Chemical species integer for atom
  !!* @param orb  Information about the orbitals for each species
  !!* @param nIntegs  Array indicating the amount of first moment integrals
  subroutine prepareOnSites(firstmom, fmOnSites, iSp, orb, nIntegs)
    real(dp), intent(out) :: firstmom(:)
    real(dp), intent(in) :: fmOnSites(:,:)
    integer, intent(in) :: iSp
    type(TOrbitals), intent(in) :: orb
    integer,  intent(in) :: nIntegs(0:,0:)

    integer :: iSh1, iSh2, iInteg, iOnSite
    integer :: l1, l2

    firstmom(:) = 0.0_dp
    iInteg = 1
    iOnSite = 0
    do iSh1 = 1, orb%nShell(iSp)
      l1 = orb%angShell(iSh1,iSp)
      do iSh2 = 1, orb%nShell(iSp)
        l2 = orb%angShell(iSh2,iSp)
        if (abs(l1 - l2) == 1) then ! <- Needs change with f-Orbitals!
          iOnSite = iOnSite + 1
          firstmom(iInteg:iInteg+nIntegs(l1,l2)-1) = fmOnSites(iOnSite,iSp)
        end if
          iInteg = iInteg + nIntegs(l1,l2)
      end do
    end do

  end subroutine prepareOnSites



  !!* Shifts the center of the dipole field for a given first moment element.
  !!* @param firstMom  First moment interaction integral in the reference system
  !!* @param shift  Factor for shifting the centers.
  !!* @param over  Overlap matrix element
  !!* @return Shifted center first moment interaction integral.
  function shiftCenter(firstMom, shift, over) result(shiftedFM)
    real(dp), intent(in) :: firstMom
    real(dp), intent(in) :: shift
    real(dp), intent(in) :: over
    real(dp) :: shiftedFM

    if (abs(shift) <= epsilon(0.0_dp) .or. &
        & abs(over) <= epsilon(0.0_dp)) then
      shiftedFM = firstMom
    else
      shiftedFM = firstMom - shift * over
    end if

  end function shiftCenter



  !!* Shifts the corresponding coordinate of the directional cosines with a
  !!* given coordinate quantum number.
  !!* @param dirCos Directional cosine vector to shift
  !!* @param iCC Quantum number of the coordinate to be shifted (1:3)
  !!* @param iLR Left or right derivative (1 -- left or 2 -- right)
  function shiftCoordinate(dirCos, iCC, iLR) result(shiftedDirCos)
    real(dp), intent(in) :: dirCos(3)
    integer, intent(in) :: iCC
    integer, intent(in) :: iLR
    real(dp) :: shiftedDirCos(3)
    real(dp), parameter :: deltaXDiff = epsilon(1.0_dp)**0.25_dp

    ! if (iCC == 1) then
    !   ii = 3
    ! else if (iCC == 2) then
    !   ii = 1
    ! else if (iCC == 3) then
    !   ii = 2
    ! end if

    shiftedDirCos(:) = dirCos(:)
    shiftedDirCos(iCC) = shiftedDirCos(iCC) + real(2 * iLR - 3, dp) * deltaXDiff

  end function shiftCoordinate



  !!* Selects the correct l-shell block and returns the
  !!* position matrix elements filled in it.
  !!* @param mat  The matrix to be filled with values
  !!* @param l1  Shell to select from the first atom
  !!* @param l2  Shell to select from the second atom
  !!* @param firstMom  First moment integrals taken from slko-table
  !!* @param fmCenter  Atom where the dipole center origin lies
  !!* @param over  Corresponding overlap integrals
  !!* @param dist  Internuclear distance
  !!* @param matCenter  Where the dipole center in the matrix block lies
  subroutine fillFMBlock(mat, l1, l2, firstMom, fmCenter, over, dist, &
      &matCenter)
    real(dp), intent(inout) :: mat(:,:,:)
    integer, intent(in) :: l1, l2
    real(dp), intent(in) :: firstMom(:)
    integer, intent(in) :: fmCenter(:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: dist
    integer, intent(in) :: matCenter

    real(dp) :: signFM, signS, shift(size(firstMom))
    integer :: iInt

    @:ASSERT(l1 <= mAngRot_)
    @:ASSERT(l2 <= mAngRot_)

    if (l1 > l2) then
      signS = (-1.0_dp)**(l1+l2)
      ! if (mod(l1, 2) /= 0 .and. mod(l2, 2) /= 0) then
      !   signFM = -1.0_dp
      ! else if (mod(l1, 2) == 0 .and. mod(l2, 2) /= 0) then
      !   signFM = 1.0_dp
      ! else if (mod(l1, 2) /= 0 .and. mod(l2, 2) == 0) then
      !   signFM = 1.0_dp
      ! else if (mod(l1, 2) == 0 .and. mod(l2, 2) == 0) then
      !   signFM = -1.0_dp
      ! else
      !   write (*,*) 'Error: Position matrix angular momenta incompatible!'
      ! end if
      signFM = (-1.0_dp)**(l1+l2+1) ! Not sure if that's correct for all cases
    else
      signS = 1.0_dp
      signFM = 1.0_dp
    end if

    do iInt = 1, size(firstMom)
      if (fmCenter(iInt) == matCenter) then
        shift(iInt) = 0.0_dp
      else if (fmCenter(iInt) == 1 .and. matCenter == 2) then
        shift(iInt) = dist
      else if (fmCenter(iInt) == 2 .and. matCenter == 1) then
        shift(iInt) = -dist
      end if
    end do

    select case (l1)
    case (0)
      select case (l2)
      case (0)
        call ss_FMBlock(mat, signFM * firstMom, shift, signS * over)
      case (1)
        call ps_FMBlock(mat, signFM * firstMom, shift, signS * over)
      ! case (2)
      !   call ds_FMBlock(mat, signFM * firstMom, shift, signS * over)
      end select
    case (1)
      select case (l2)
      case (0)
        call ps_FMBlock(mat, signFM * firstMom, shift, signS * over)
      case (1)
        call pp_FMBlock(mat, signFM * firstMom, shift, signS * over)
      ! case (2)
      !   call dp_FMBlock(mat, signFM * firstMom, shift, signS * over)
      end select
    ! case (2)
    !   select case (l2)
    !   case (0)
    !     call ds_FMBlock(mat, signFM * firstMom, shift, signS * over)
    !   case (1)
    !     call dp_FMBlock(mat, signFM * firstMom, shift, signS * over)
    !   case (2)
    !     call dd_FMBlock(mat, signFM * firstMom, shift, signS * over)
    !   end select
    end select

  end subroutine fillFMBlock



  !! SECTION: Diatomic first moment blocks
  !! Helper routines to set up diatomic block matrices with values from
  !! the integration frame.

  !!* Creates an diatomic s-s-Orbital first moment block.
  !!* @param mat  Rectangular block matrix describing the first moment
  !!* interaction between two orbitals
  !!* @param firstMom  First moment integrals for the orbital combination
  !!* @param shift  Factor for shifting the center of the dipole
  !!* @param over  Overlap integrals for the combination
  subroutine ss_FMBlock(mat, firstMom, shift, over)
    real(dp), intent(inout) :: mat(:,:,:)
    real(dp), intent(in) :: firstMom(:)
    real(dp), intent(in) :: shift(:)
    real(dp), intent(in) :: over(:)

    @:ASSERT(all(shape(mat) == (/ 1, 1, 3 /)))
    @:ASSERT(size(firstMom) == 1)
    @:ASSERT(size(shift) == 1)
    @:ASSERT(size(over) == 1)

    mat(:,:,:) = 0.0_dp

    !! third index 2 => z-Dipole (m =  0)
    mat(1,1,2) = shiftCenter(FirstMom(1), shift(1), over(1))

  end subroutine ss_FMBlock



  !!* Creates an diatomic s-p-Orbital first moment block.
  !!* @param mat  Rectangular block matrix describing the first moment
  !!* interaction between two orbitals
  !!* @param firstMom  First moment integrals for the orbital combination
  !!* @param shift  Factor for shifting the center of the dipole
  !!* @param over  Overlap integrals for the combination
  subroutine sp_FMBlock(mat, firstMom, shift, over)
    real(dp), intent(inout) :: mat(:,:,:)
    real(dp), intent(in) :: firstMom(:)
    real(dp), intent(in) :: shift(:)
    real(dp), intent(in) :: over(:)

    @:ASSERT(all(shape(mat) == (/ 1, 3, 3 /)))
    @:ASSERT(size(firstMom) == 2)
    @:ASSERT(size(shift) == 2)
    @:ASSERT(size(over) == 1)

    mat(:,:,:) = 0.0_dp

    !! third index 2 => y-Dipole (m = -1)
    mat(1,1,1) = FirstMom(2) ! <p1|1|s0>
    !! third index 2 => z-Dipole (m =  0)
    mat(1,2,2) = shiftCenter(FirstMom(1), shift(1), over(1)) ! <p0|0|s0>
    !! third index 3 => x-Dipole (m = +1)
    mat(1,3,3) = FirstMom(2) ! <p1|1|s0>

  end subroutine sp_FMBlock



  !!* Creates an diatomic p-s-Orbital first moment block.
  !!* @param mat  Rectangular block matrix describing the first moment
  !!* interaction between two orbitals
  !!* @param firstMom  First moment integrals for the orbital combination
  !!* @param shift  Factor for shifting the center of the dipole
  !!* @param over  Overlap integrals for the combination
  subroutine ps_FMBlock(mat, firstMom, shift, over)
    real(dp), intent(inout) :: mat(:,:,:)
    real(dp), intent(in) :: firstMom(:)
    real(dp), intent(in) :: shift(:)
    real(dp), intent(in) :: over(:)

    @:ASSERT(all(shape(mat) == (/ 3, 1, 3 /)))
    @:ASSERT(size(firstMom) == 2)
    @:ASSERT(size(shift) == 2)
    @:ASSERT(size(over) == 1)

    mat(:,:,:) = 0.0_dp

    !! third index 1 => y-Dipole (m = -1)
    mat(1,1,1) = FirstMom(2) ! <p1|1|s0>
    !! third index 2 => z-Dipole (m = 0)
    mat(2,1,2) = shiftCenter(FirstMom(1), shift(1), over(1)) ! <p0|0|s0>
    !! third index 3 => x-Dipole (m = 1)
    mat(3,1,3) = FirstMom(2) ! <p1|1|s0>

  end subroutine ps_FMBlock



  !!* Creates an diatomic p-p-Orbital first moment block.
  !!* @param mat  Rectangular block matrix describing the first moment
  !!* interaction between two orbitals
  !!* @param firstMom  First moment integrals for the orbital combination
  !!* @param shift  Factor for shifting the center of the dipole
  !!* @param over  Overlap integrals for the combination
  subroutine pp_FMBlock(mat, firstMom, shift, over)
    real(dp), intent(inout) :: mat(:,:,:)
    real(dp), intent(in) :: firstMom(:)
    real(dp), intent(in) :: shift(:)
    real(dp), intent(in) :: over(:)

    real(dp) :: rTmp

    @:ASSERT(all(shape(mat) == (/ 3, 3, 3 /)))
    @:ASSERT(size(firstMom) == 4)
    @:ASSERT(size(shift) == 4)
    @:ASSERT(size(over) == 2)

    mat(:,:,:) = 0.0_dp

    rTmp = shiftCenter(FirstMom(4), shift(4), over(2))
    !! third index 1 => y-Dipole (m = -1)
    mat(2,1,1) = FirstMom(2) ! <p1|1|p0>
    mat(1,2,1) = FirstMom(3) ! <p0|1|p1>
    !! third index 2 => z-Dipole (m = 0)
    mat(1,1,2) = rTmp ! <p1|0|p1>
    mat(2,2,2) = shiftCenter(FirstMom(1), shift(1), over(1)) ! <p0|0|p0>
    mat(3,3,2) = rTmp ! <p1|0|p1>
    !! third index 3 => x-Dipole (m = 1)
    mat(3,2,3) = firstMom(3) ! <p0|1|p1>
    mat(2,3,3) = firstMom(2) ! <p1|1|p0>

  end subroutine pp_FMBlock

  !! END SECTION: Diatomic first moment blocks

end module positionmatrix
