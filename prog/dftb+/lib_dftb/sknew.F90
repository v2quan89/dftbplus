!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains code to perform the sk rotations of matrix elements from the parameterization
!> orientation along <0,0,1> to the one needed in the actual calculation.
!>
!> To do: Transformations to give the derivatives with respect to ll, mm and nn.
!> Base on "Compact expression for the angular dependence of tight-binding hamiltonian matrix
!> elements", A. V. Podolskiy and P. Vogl, Phys. Rev.  B 69, 233101 (2004).
!>
!> Caveat: Only angular momenta up to f are currently allowed
module sknew
  use assert
  use accuracy
  use commontypes
  implicit none

  private

  public :: rotateH0, rotateFM

  !! individual rotation matrices only needed to be public, when testing.
  public :: rotYZ_p, rotYZ_d, rotYZ_f
  public :: rotXY_p, rotXY_d, rotXY_f
  ! public :: rotYZ_s, rotYZ_p, rotYZ_d, rotYZ_f
  ! public :: rotXY_s, rotXY_p, rotXY_d, rotXY_f

  ! Maximal angular momentum, for which rotations are present
  integer, parameter :: mAngRot_ = 3
  integer, parameter :: mOrbPerL_ = 2 * mAngRot_ + 1

  ! Angular momenta for rotation matrices.
  integer, parameter :: iRotFM_ = 1

  ! Treat sign convention as in the old rotation routines
  logical, parameter :: oldCompat_ = .true.

contains


  !> Driver for making the non-SCC hamiltonian or overlap matrices
  !! for a given diatomic block.
  !!
  !! \param hh  Rectangular matrix containing the resulting diatomic
  !!     matrix elements.
  !! \param skIntegs  Slater-Koster table for dimer of species i-j.
  !! \param ll  Directional cosine ll.
  !! \param mm  Directional cosine mm.
  !! \param nn  Directional cosine nn.
  !! \param iSp1  Chemical species of atom i.
  !! \param iSp2  Chemical species of atom j.
  !! \param orb  Information about the orbitals of chemical species in the
  !!     system.
  !!
  !> Driver for making the non-SCC hhamiltonian or overlap matrices for a given diatomic block
  !> Caveat: Only angular momenta up to f are currently allowed
  subroutine rotateH0(skIntegs, ll, mm, nn, iSp1, iSp2, orb, hh)
    real(dp), intent(in), target :: skIntegs(:)
    real(dp), intent(in) :: ll, mm, nn
    integer, intent(in) :: iSp1, iSp2
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(out) :: hh(:,:)

    real(dp) :: rotMtxPerL(mOrbPerL_, mOrbPerL_, 0:mAngRot_)
    integer :: maxL

    @:ASSERT(all(shape(hh) >= [ orb%nOrbSpecie(iSp2), orb%nOrbSpecie(iSp1) ]))

    maxL = max(maxval(orb%angShell(:,iSp1)), maxval(orb%angShell(:,iSp2)))
    call getRotationMatricesUpToL(ll, mm, nn, maxL, rotMtxPerL(:,:,0:maxL))
    call createDiatomicBlock(iSp1, iSp2, orb, skIntegs, rotMtxPerL, hh)

  end subroutine rotateH0


  !> Rotates a specific diatomic block of any position matrix into system
  !! coordinates.
  !! \param mat Position matrix to be rotated
  !! \param ll  Directional cosine ll.
  !! \param mm  Directional cosine mm.
  !! \param nn  Directional cosine nn.
  !! \caveat Only angular momenta up to $d$ are currently allewed.
  !! \note As opposed to rotateH0 this routine is strictly rotating, not
  !! building the whole matrix, see the positionmatrix module (might be
  !! changed in refacturing).
  subroutine rotateFM(mat, ll, mm, nn)
    real(dp), intent(inout) :: mat(:,:,:)
    real(dp), intent(in) :: ll, mm, nn

    real(dp) :: rotMtxM(2 * iRotFM_ + 1, 2 * iRotFM_ + 1)
    real(dp) :: refFrame(size(mat, dim=1), size(mat, dim=2), size(mat, dim=3))
    real(dp), target :: rotMtxPerL(mOrbPerL_, mOrbPerL_, 0:mAngRot_)
    real(dp), pointer :: pRotMtx1(:,:), pRotMtx2(:,:)
    integer :: l1, l2, maxL

    l1 = (size(mat, dim=1) - 1) / 2
    l2 = (size(mat, dim=2) - 1) / 2
    maxL = max(l1, l2)
    call getRotationMatricesUpToL(ll, mm, nn, maxL, rotMtxPerL(:,:,0:maxL))
    call getRotationMatrixForL(ll, mm, nn, iRotFM_, rotMtxM)
    refFrame(:,:,:) = mat(:,:,:)
    pRotMtx1 => rotMtxPerL(1:size(mat, dim=1),1:size(mat, dim=1), l1)
    pRotMtx2 => rotMtxPerL(1:size(mat, dim=2),1:size(mat, dim=2), l2)
    call rotate3(refFrame, pRotMtx1, pRotMtx2, rotMtxM, mat)

  end subroutine rotateFM



  !> Returns rotation matrices for a direction up to given angular momentum.
  !!
  !! \param uu  Direction cosine along x-axis.
  !! \param vv  Direction cosine along y-axis.
  !! \param ww  Direction cosine along z-axis.
  !! \param maxL  Maximal angular momentum to consider.
  !! \param rotMtxPerL  Matrix of the rotation matrices. Shape:
  !!     (2 * maxL + 1, 2 * maxL + 1, 0:maxL)
  !!
  subroutine getRotationMatricesUpToL(uu, vv, ww, maxL, rotMtxPerL)
    real(dp), intent(in) :: uu, vv, ww
    integer, intent(in) :: maxL
    real(dp), intent(out) :: rotMtxPerL(:,:,0:)

    integer :: iL, nOrb

    @:ASSERT(all(shape(rotMtxPerL) >= [ 2 * maxL + 1, 2 * maxL + 1, maxL + 1 ]))

    rotMtxPerL(:,:,:) = 0.0_dp
    do iL = 0, maxL
      nOrb = 2 * iL + 1
      call getRotationMatrixForL(uu, vv, ww, iL, rotMtxPerL(1:nOrb, 1:nOrb, iL))
    end do

  end subroutine getRotationMatricesUpToL


  !> Returns rotation matrices for a direction for a given angular momentum.
  !!
  !! \param uu  Direction cosine along x-axis.
  !! \param vv  Direction cosine along y-axis.
  !! \param ww  Direction cosine along z-axis.
  !! \param iL  Angular momentum to consider.
  !! \param rotMtxPerL  Rotatinon matrix. Shape: (2 * maxL + 1, 2 * maxL + 1)
  !!
  subroutine getRotationMatrixForL(uu, vv, ww, iL, rotMtx)
    real(dp), intent(in) :: uu, vv, ww
    integer, intent(in) :: iL
    real(dp), intent(out) :: rotMtx(:,:)

    @:ASSERT(all(shape(rotMtx) == [ 2 * iL + 1, 2 * iL + 1 ]))
    if (abs(vv) > abs(ww)) then
      select case (iL)
      case(0)
        ! call rotYZ_s(uu, vv, ww, rotMtx)
        rotMtx(1,1) = 1.0_dp
      case(1)
        call rotYZ_p(uu, vv, ww, rotMtx)
      case(2)
        call rotYZ_d(uu, vv, ww, rotMtx)
      case(3)
        call rotYZ_f(uu, vv, ww, rotMtx)
      end select
    else
      select case(iL)
      case(0)
        ! call rotXY_s(uu, vv, ww, rotMtx)
        rotMtx(1,1) = 1.0_dp
      case(1)
        call rotXY_p(uu, vv, ww, rotMtx)
      case(2)
        call rotXY_d(uu, vv, ww, rotMtx)
      case(3)
        call rotXY_f(uu, vv, ww, rotMtx)
      end select
    end if

  end subroutine getRotationMatrixForL


  !> Creates a diatomic block in the Hamiltonian/overlap.
  !!
  !! \param iSp1  Specie of 1st atom.
  !! \param iSp2  Specie of 2nd atom.
  !! \param orb  Orbital informations.
  !! \param skIntegs  SK-integrals in the reference frame between the atoms.
  !! \param rotMtxPerL  Rotation matrices up to the highest angular momentum
  !!     needed in the calculation.
  !! \param diatomicBlock  Diatomic block in the actual frame.
  !!
  subroutine createDiatomicBlock(iSp1, iSp2, orb, skIntegs, rotMtxPerL, &
      & diatomicBlock)
    integer, intent(in) :: iSp1, iSp2
    type(TOrbitals), intent(in) :: orb
    real(dp), target, intent(in) :: skIntegs(:), rotMtxPerL(:,:,0:)
    real(dp), intent(out) :: diatomicBlock(:,:)

    integer :: ind, iRow, iCol, iSh1, iSh2
    integer :: nOrb1, nOrb2, ang1, ang2
    real(dp), pointer :: pRotMtx1(:,:), pRotMtx2(:,:), pRotMtxM(:,:)
    real(dp), target :: diatomicBlockRefFrame(mOrbPerL_, mOrbPerL_, 1)
    real(dp), target :: diatomicBlockRotated(mOrbPerL_, mOrbPerL_, 1)
    real(dp), pointer :: pBlockRefFrame(:,:,:), pBlockRotated(:,:,:), pSK(:)

    diatomicBlock(:,:) = 0.0_dp
    ind = 1
    iCol = 1
    pRotMtxM => rotMtxPerL(:1,:1,0)
    do iSh1 = 1, orb%nShell(iSp1)
      ang1 = orb%angShell(iSh1, iSp1)
      nOrb1 = 2 * ang1 + 1
      iRow = 1
      pRotMtx1 => rotMtxPerL(1:nOrb1, 1:nOrb1, ang1)
      do iSh2 = 1, orb%nShell(iSp2)
        ang2 = orb%angShell(iSh2, iSp2)
        nOrb2 = 2 * ang2 + 1
        pRotMtx2 => rotMtxPerL(1:nOrb2, 1:nOrb2, ang2)
        pSK => skIntegs(ind : ind + min(ang1, ang2))
        pBlockRefFrame => diatomicBlockRefFrame(1:nOrb2, 1:nOrb1, :)
        pBlockRotated => diatomicBlockRotated(1:nOrb2, 1:nOrb1, :)
        call getDiatomicBlockInRefFrame(pSK, ang2, ang1, pBlockRefFrame(:,:,1))
        ! call rotate2(pBlockRefFrame, pRotMtx2, pRotMtx1, pBlockRotated)
        call rotate3(pBlockRefFrame, pRotMtx2, pRotMtx1, pRotMtxM, pBlockRotated)
        diatomicBlock(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = pBlockRotated(:,:,1)
        ind = ind + min(ang1, ang2) + 1
        iRow = iRow + nOrb2
      end do
      iCol = iCol + nOrb1
    end do

  end subroutine createDiatomicBlock


  !> Returns the diatomic block in the reference frame.
  !!
  !! \param integs  SK-integrals between orbitals with given angular momenta
  !!     (in the order sigma, pi, delta, ...)
  !! \param ang1  Angular momentum of first tesseral.
  !! \param ang2  Angular momentum of second tesseral.
  !! \param blockRefFrame  Diatomic block in the reference frame.
  !!
  subroutine getDiatomicBlockInRefFrame(integs, ang1, ang2, blockRefFrame)
    real(dp), intent(in) :: integs(:)
    integer, intent(in) :: ang1, ang2
    real(dp), intent(out) :: blockRefFrame(:,:)

    integer :: l1, l2, mm, row0, col0, parity1, parity2
    real(dp), target :: work(mOrbPerL_, mOrbPerL_)
    real(dp), pointer :: pBlock(:,:)
    logical :: transposed

    @:ASSERT(size(integs) == min(ang1, ang2) + 1)
    @:ASSERT(all(shape(blockRefFrame) == [ 2 * ang1 + 1, 2 * ang2 + 1 ]))

    l1 = min(ang1, ang2)
    l2 = max(ang1, ang2)
    pBlock => work(1:2*l1+1, 1:2*l2+1)
    pBlock(:,:) = 0.0_dp

    row0 = l1 + 1
    col0 = l2 + 1
    pBlock(row0, col0) = integs(1)
    do mm = 1, l1
      pBlock(row0 - mm, col0 - mm) = integs(mm + 1)
      pBlock(row0 + mm, col0 + mm) = integs(mm + 1)
    end do

    transposed = (ang1 > ang2)
    if (transposed) then
      blockRefFrame(:,:) = transpose(pBlock)
    else
      blockRefFrame(:,:) = pBlock
    end if

    ! Theoretically, the transposed block should be sign inverted if the orbital
    ! product has negative (odd) parity. However, the integrals in skIntegs have
    ! already a negative sign in those cases, so we invert the non-transposed
    ! matrix instead, to be sign compatible with the old routines.
    if (transposed .neqv. oldCompat_) then
      parity1 = 1 - 2 * mod(ang1, 2)
      parity2 = 1 - 2 * mod(ang2, 2)
      if (parity1 * parity2 < 0) then
        blockRefFrame(:,:) = -blockRefFrame
      end if
    end if

  end subroutine getDiatomicBlockInRefFrame


  !> Rotate a diatomic block with two indices (Hamiltonian/overlap)
  !!
  !! \param blockRefFrame  Diatomic block in the reference frame.
  !! \param rotMtx1  Rotation matrix for the tesserals of atom 1.
  !! \param rotMtx2  Rotation matrix for the tesserals of atom 2.
  !! \param blockRotated  Rotated block.
  !!
  subroutine rotate2(blockRefFrame, rotMtx1, rotMtx2, blockRotated)
    real(dp), intent(in) :: blockRefFrame(:,:), rotMtx1(:,:), rotMtx2(:,:)
    real(dp), intent(out) :: blockRotated(:,:)

    @:ASSERT(all(shape(blockRefFrame) == shape(blockRotated)))
    @:ASSERT(all(shape(rotMtx1) == [ size(blockRefFrame, dim=1), size(blockRefFrame, dim=1)]))
    @:ASSERT(all(shape(rotMtx2) == [ size(blockRefFrame, dim=2), size(blockRefFrame, dim=2)]))

    blockRotated(:,:) = matmul(matmul(transpose(rotMtx1), blockRefFrame), &
        & rotMtx2)

  end subroutine rotate2



  !> Rotate a diatomic block with three indices (first moment)
  !!
  !! \param blockRefFrame  Diatomic block in the reference frame.
  !! \param rotMtx1  Rotation matrix for the tesserals of atom 1.
  !! \param rotMtx2  Rotation matrix for the tesserals of atom 2.
  !! \param rotMtxM  Rotation matrix for the tesserals of order of the matrix.
  !! \param blockRotated  Rotated block.
  !!
  subroutine rotate3(blockRefFrame, rotMtx1, rotMtx2, rotMtxM, blockRotated)
    real(dp), intent(in) :: blockRefFrame(:,:,:)
    real(dp), intent(in) :: rotMtx1(:,:), rotMtx2(:,:), rotMtxM(:,:)
    real(dp), intent(out) :: blockRotated(:,:,:)

    real(dp) :: blockOrbRotated(size(blockRotated, dim=1), &
        &size(blockRotated, dim=2), size(blockRotated, dim=3))
    integer :: iMM, iMMP, mMM

    @:ASSERT(all(shape(blockRefFrame) == shape(blockRotated)))
    @:ASSERT(all(shape(rotMtx1) == [ size(blockRefFrame, dim=1), size(blockRefFrame, dim=1)]))
    @:ASSERT(all(shape(rotMtx2) == [ size(blockRefFrame, dim=2), size(blockRefFrame, dim=2)]))
    @:ASSERT(all(shape(rotMtxM) == [ size(blockRefFrame, dim=3), size(blockRefFrame, dim=3)]))

    mMM = size(blockRefFrame, dim=3)
    blockRotated(:,:,:) = 0.0_dp

    ! print *
    ! print *, "blockRefFrame"
    ! print *, blockRefFrame(:,:,:)

    !! Rotation of the orbitals portion of the matrix (first two indices)
    do iMM = 1, mMM
      call rotate2(blockRefFrame(:,:,iMM), rotMtx1, rotMtx2, &
          &blockOrbRotated(:,:,iMM))
    end do

    ! print *
    ! print *, "blockOrbRotated"
    ! print *, blockOrbRotated(:,:,:)

    !! Rotation of the third index of the matrix
    do iMM = 1, mMM
      do iMMP = 1, mMM
        blockRotated(:,:,iMM) = blockRotated(:,:,iMM) + rotMtxM(iMMP,iMM) * &
            blockOrbRotated(:,:,iMMP)
      end do
    end do

    ! print *
    ! print *, "BlockRotated"
    ! print *, blockRotated(:,:,:)

  end subroutine rotate3



  !> Rotation matrix from y,z rotation.
  !!
  ! subroutine rotYZ_s(uu, vv, ww, rot)
  !   real(dp), intent(in) :: uu, vv, ww
  !   real(dp), intent(out) :: rot(:,:)

  !   @:ASSERT(all(shape(rot) == [ 1, 1 ]))
  !   rot(1,1) = 1.0_dp

  ! end subroutine rotYZ_s


  !> Rotation matrix from derived from y,z rotation.
  !!
  !! \caveat It gets numerical unstable if ww is close to 1.0!
  !!
  subroutine rotYZ_p(uu, vv, ww, rot)
    real(dp), intent(in) :: uu, vv, ww
    real(dp), intent(out) :: rot(:,:)

    real(dp) :: w2c12

    @:ASSERT(all(shape(rot) == [ 3, 3 ]))

    w2c12 = sqrt(1.0_dp - ww * ww)

    rot(1,1) = uu / w2c12
    rot(2,1) = vv
    rot(3,1) = vv * ww / w2c12

    rot(1,2) = 0.0_dp
    rot(2,2) = ww
    rot(3,2) = -w2c12

    rot(1,3) = -vv / w2c12
    rot(2,3) = uu
    rot(3,3) = uu * ww / w2c12

  end subroutine rotYZ_p


  !> Rotation matrix from derived from y,z rotation.
  !!
  !! \caveat It gets numerical unstable if ww is close to 1.0!
  !!
  subroutine rotYZ_d(uu, vv, ww, rot)
    real(dp), intent(in) :: uu, vv, ww
    real(dp), intent(out) :: rot(:,:)

    real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
    real(dp) :: w2c, w2c12, u2

    @:ASSERT(all(shape(rot) == [ 5, 5 ]))

    w2c = 1.0_dp - ww * ww
    w2c12 = sqrt(w2c)
    u2 = uu * uu

    rot(1,1) = -ww + 2.0_dp * u2 * ww / w2c
    rot(2,1) = -w2c12 + 2.0_dp * u2 / w2c12
    rot(3,1) = uu * vv * sqrt3
    rot(4,1) = 2.0_dp * uu * vv * ww / w2c12
    rot(5,1) = -uu * vv + 2.0_dp * uu * vv / w2c

    rot(1,2) = -uu
    rot(2,2) = uu * ww / w2c12
    rot(3,2) = vv * ww * sqrt3
    rot(4,2) = -2.0_dp * vv * w2c12 + vv / w2c12
    rot(5,2) = -ww * vv

    rot(1,3) = 0.0_dp
    rot(2,3) = 0.0_dp
    rot(3,3) = 1.0_dp - 1.5_dp * w2c
    rot(4,3) = -sqrt3 * ww * w2c12
    rot(5,3) = 0.5_dp * sqrt3 * w2c

    rot(1,4) = vv
    rot(2,4) = -vv * ww / w2c12
    rot(3,4) = uu * ww * sqrt3
    rot(4,4) = -2.0_dp * uu * w2c12 + uu / w2c12
    rot(5,4) = -uu * ww

    rot(1,5) = -2.0_dp * uu * vv * ww / w2c
    rot(2,5) = -2.0_dp * uu * vv / w2c12
    rot(3,5) = sqrt3 * (-0.5_dp * w2c + u2)
    rot(4,5) = -ww * w2c12 + 2.0_dp * u2 * ww / w2c12
    rot(5,5) = 0.5_dp * w2c - 1.0_dp + u2 * (-1.0_dp + 2.0_dp / w2c)

  end subroutine rotYZ_d


  !> Rotation matrix from derived from y,z rotation.
  !!
  !! \caveat It gets numerical unstable if ww is close to 1.0!
  !!
  subroutine rotYZ_f(uu, vv, ww, rot)
    real(dp), intent(in) :: uu, vv, ww
    real(dp), intent(out) :: rot(:,:)

    real(dp), parameter :: sqrt6 = sqrt(6.0_dp)
    real(dp), parameter :: sqrt10 = sqrt(10.0_dp)
    real(dp), parameter :: sqrt15 = sqrt(15.0_dp)

    real(dp) :: w2c, w2c12, w2c32, u2

    @:ASSERT(all(shape(rot) == [ 7, 7 ]))

    w2c = 1.0_dp - ww * ww
    w2c12 = sqrt(w2c)
    w2c32 = w2c12**3
    u2 = uu * uu

    rot(1,1) = uu * (2.25_dp * w2c12 - 3.0_dp * (u2 + 1.0_dp) / w2c12 &
        & + 4.0_dp * u2 / w2c32)
    rot(2,1) = sqrt6 * uu * ww * (-1.5_dp + 2.0_dp * u2 / w2c)
    rot(3,1) = sqrt15 * uu * (-0.75_dp * w2c12 + u2 / w2c12)
    rot(4,1) = vv * sqrt10 * (-0.25_dp * w2c + u2)
    rot(5,1) = sqrt15 * vv * ww * (-0.25_dp * w2c12 + u2 / w2c12)
    rot(6,1) = sqrt6 * vv * (0.25_dp * w2c - 0.25_dp * (4.0_dp * u2 + 2.0_dp) &
        & + 2.0_dp * u2 / w2c)
    rot(7,1) = vv * ww * (0.25_dp * w2c12 - (u2 + 1.0_dp) / w2c12 &
        & + 4.0_dp * u2 / w2c32)

    rot(1,2) = sqrt6 * ww * (0.5_dp * w2c12 - u2 / w2c12)
    rot(2,2) = 2.0_dp * w2c - 1.0_dp + u2 * (-4.0_dp + 2.0_dp / w2c)
    rot(3,2) = sqrt10 * ww * (-0.5_dp * w2c12 + u2 / w2c12)
    rot(4,2) = uu * vv * ww * sqrt15
    rot(5,2) = sqrt10 * vv * uu * (-1.5_dp * w2c12 +  1.0_dp / w2c12)
    rot(6,2) = uu * vv * ww * (-3.0_dp + 2.0_dp / w2c)
    rot(7,2) = sqrt6 * vv * uu * (0.5_dp * w2c12 - 1.0_dp / w2c12)

    rot(1,3) = 0.25_dp * sqrt15 * uu * w2c12
    rot(2,3) = -0.5_dp * sqrt10 * uu * ww
    rot(3,3) = uu * (-1.25_dp * w2c12 + 1.0_dp / w2c12)
    rot(4,3) = sqrt6 * vv * (-1.25_dp * w2c + 1.0_dp)
    rot(5,3) = vv * ww * (-3.75_dp * w2c12 + 1.0_dp / w2c12)
    rot(6,3) = sqrt10 * vv * (0.75_dp * w2c - 0.5_dp)
    rot(7,3) = 0.25_dp * sqrt15 * vv * ww * w2c12

    rot(1,4) = 0.0_dp
    rot(2,4) = 0.0_dp
    rot(3,4) = 0.0_dp
    rot(4,4) = ww * (-2.5_dp * w2c + 1.0_dp)
    rot(5,4) = sqrt6 * (1.25_dp * w2c32 - w2c12)
    rot(6,4) = 0.5_dp *sqrt15 * w2c * ww
    rot(7,4) = -0.25_dp * sqrt10 * w2c32

    rot(1,5) = -0.25_dp * sqrt15 * vv *  w2c12
    rot(2,5) = 0.5_dp * sqrt10 * vv * ww
    rot(3,5) =  vv * (1.25_dp * w2c12 - 1.0_dp / w2c12)
    rot(4,5) = sqrt6 * uu * (-1.25_dp * w2c + 1.0_dp)
    rot(5,5) = uu * ww * (-3.75_dp * w2c12 + 1.0_dp / w2c12)
    rot(6,5) = sqrt10 * uu * (0.75_dp * w2c - 0.5_dp)
    rot(7,5) = 0.25_dp * sqrt15 * uu * ww * w2c12

    rot(1,6) = sqrt6 * uu * vv * ww / w2c12
    rot(2,6) = uu * vv * (4.0_dp - 2.0_dp / w2c)
    rot(3,6) = - sqrt10 * uu * vv * ww / w2c12
    rot(4,6) = sqrt15 * ww * (-0.5_dp * w2c + u2)
    rot(5,6) = sqrt10 * (0.75_dp * w2c32 &
        & - 0.25_dp * (6.0_dp * u2 + 2.0_dp) * w2c12 + u2 / w2c12)
    rot(6,6) = ww * (1.5_dp * w2c - 0.5_dp * (6.0_dp * u2 + 2.0_dp) &
        & + 2.0_dp * u2 / w2c)
    rot(7,6) = sqrt6 * (-0.25_dp * w2c32 &
        & + 0.25_dp * (2.0_dp * u2 + 2.0_dp) * w2c12 - u2 / w2c12)

    rot(1,7) = vv * (-0.75_dp * w2c12 &
        & + 0.25_dp * (12.0_dp * u2 + 4.0_dp) / w2c12 - 4.0_dp * u2 / w2c32)
    rot(2,7) = sqrt6 * vv * ww * (0.5_dp - 2.0_dp * u2 / w2c)
    rot(3,7) = sqrt15 * vv * (0.25_dp * w2c12 - u2 / w2c12)
    rot(4,7) = sqrt10 * uu * (-0.75_dp * w2c + u2)
    rot(5,7) = sqrt15 * uu * ww * (-0.75_dp * w2c12 + u2 / w2c12)
    rot(6,7) = sqrt6 * uu * (0.75_dp * w2c &
        & - 0.25_dp * (4.0_dp * u2 + 6.0_dp) + 2.0_dp * u2 / w2c)
    rot(7,7) = uu * ww * (0.75_dp * w2c12 &
        & - 0.25_dp * (4.0_dp * u2 + 12.0_dp) / w2c12 + 4.0_dp * u2 / w2c32)

  end subroutine rotYZ_f


  !> Rotation matrix from derived from x,y rotation.
  !!
  ! subroutine rotXY_s(uu, vv, ww, rot)
  !   real(dp), intent(in) :: uu, vv, ww
  !   real(dp), intent(out) :: rot(:,:)

  !   @:ASSERT(all(shape(rot) == [ 1, 1 ]))

  !   rot(1,1) = 1.0_dp

  ! end subroutine rotXY_s


  !> Rotation matrix from derived from x,y rotation.
  !!
  !! \caveat It gets numerical unstable if vv is close to 1.0!
  !!
  subroutine rotXY_p(uu, vv, ww, rot)
    real(dp), intent(in) :: uu, vv, ww
    real(dp), intent(out) :: rot(:,:)

    real(dp) :: v2c12

    @:ASSERT(all(shape(rot) == [ 3, 3 ]))

    v2c12 = sqrt(1.0_dp - vv * vv)
    rot(1,1) = v2c12
    rot(2,1) = vv

    rot(3,1) = 0.0_dp
    rot(1,2) = -vv * ww / v2c12
    rot(2,2) = ww

    rot(3,2) = -uu / v2c12
    rot(1,3) = -uu * vv / v2c12
    rot(2,3) = uu
    rot(3,3) = ww / v2c12

  end subroutine rotXY_p


  !> Rotation matrix from derived from x,y rotation.
  !!
  !! \caveat It gets numerical unstable if vv is close to 1.0!
  !!
  subroutine rotXY_d(uu, vv, ww, rot)
    real(dp), intent(in) :: uu, vv, ww
    real(dp), intent(out) :: rot(:,:)

    real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
    real(dp) :: v2c, v2c12

    @:ASSERT(all(shape(rot) == [ 5, 5 ]))

    v2c = 1.0_dp - vv * vv
    v2c12 = sqrt(v2c)

    rot(1,1) = ww
    rot(2,1) = uu * (2.0_dp * v2c12 - 1.0_dp / v2c12)
    rot(3,1) = uu * vv * sqrt3
    rot(4,1) = vv * ww / v2c12
    rot(5,1) = uu * vv

    rot(1,2) = -uu
    rot(2,2) = ww * (2.0_dp * v2c12 - 1.0_dp / v2c12)
    rot(3,2) = vv * ww * sqrt3
    rot(4,2) = -uu * vv / v2c12
    rot(5,2) = vv * ww

    rot(1,3) = uu * vv * ww * sqrt3 / v2c
    rot(2,3) = -vv * ww**2 * sqrt3 / v2c12
    rot(3,3) = 1.5_dp * ww**2 - 0.5_dp
    rot(4,3) = -uu * ww * sqrt3 / v2c12
    rot(5,3) = sqrt3 * (0.5_dp * (ww**2 + 1.0_dp) - ww**2 / v2c)

    rot(1,4) = vv * (1.0_dp - 2.0_dp * ww**2 / v2c)
    rot(2,4) = -2.0_dp * vv * uu * ww / v2c12
    rot(3,4) = uu * ww * sqrt3
    rot(4,4) = -v2c12 + 2.0_dp * ww**2 / v2c12
    rot(5,4) = uu * ww - 2.0_dp * uu * ww / v2c

    rot(1,5) = -uu * vv * ww / v2c
    rot(2,5) = vv * (-2.0_dp * v2c12 + ww**2 / v2c12)
    rot(3,5) = 0.5_dp * sqrt3 * (2.0_dp * v2c - ww**2 - 1.0_dp)
    rot(4,5) = uu * ww / v2c12
    rot(5,5) = v2c - 0.5_dp * ww**2 - 0.5_dp + ww**2 / v2c

  end subroutine rotXY_d


  !> Rotation matrix from derived from x,y rotation.
  !!
  !! \caveat It gets numerical unstable if vv is close to 1.0!
  !!
  subroutine rotXY_f(uu, vv, ww, rot)
    real(dp), intent(in) :: uu, vv, ww
    real(dp), intent(out) :: rot(:,:)

    real(dp), parameter :: sqrt6 = sqrt(6.0_dp)
    real(dp), parameter :: sqrt10 = sqrt(10.0_dp)
    real(dp), parameter :: sqrt15 = sqrt(15.0_dp)

    real(dp) :: v2c, v2c12, v2c32, w2

    @:ASSERT(all(shape(rot) == [ 7, 7 ]))

    v2c = 1.0_dp - vv * vv
    v2c12 = sqrt(v2c)
    v2c32 = v2c12**3
    w2 = ww * ww

    rot(1,1) = v2c32 + (-0.75_dp * w2 - 0.75_dp) * v2c12 + 1.5_dp * w2 / v2c12
    rot(2,1) = sqrt6 * uu * ww * (1.0_dp - 0.5_dp / v2c)
    rot(3,1) = sqrt15 * (v2c32 - 0.75_dp * (w2 + 1.0_dp) * v2c12 &
        & + 0.5_dp * w2 / v2c12)
    rot(4,1) = sqrt10 * vv * (v2c - 0.75_dp * w2 - 0.25_dp)
    rot(5,1) = 0.5_dp * sqrt15 * uu * vv * ww / v2c12
    rot(6,1) = sqrt6 * vv * (v2c - 0.25_dp * (3.0_dp * w2 + 1.0_dp) &
        & + 0.5_dp * w2 / v2c)
    rot(7,1) = 1.5_dp * uu * vv * ww / v2c12

    rot(1,2) = sqrt6 * uu * ww * (0.5_dp * v2c12 - 1.0_dp / v2c12)
    rot(2,2) = -2.0_dp * v2c + 4.0_dp * w2 + 1.0_dp - 2.0_dp * w2 / v2c
    rot(3,2) = sqrt10 * uu * ww * (1.5_dp * v2c12 - 1.0_dp / v2c12)
    rot(4,2) = sqrt15 * uu * vv * ww
    rot(5,2) = sqrt10 * vv * (-0.5_dp * v2c12 + w2 / v2c12)
    rot(6,2) = uu * vv * ww * (3.0_dp - 2.0_dp / v2c)
    rot(7,2) = sqrt6 * vv * (-0.5_dp * v2c12 + w2 / v2c12)

    rot(1,3) = sqrt15 * (0.25_dp * (w2 + 1.0_dp) * v2c12 &
        & - 0.5_dp * w2 / v2c12)
    rot(2,3) = sqrt10 * uu * ww * (-1.0_dp + 0.5_dp / v2c)
    rot(3,3) = (3.75_dp * w2 - 0.25_dp) * v2c12 - 2.5_dp * w2 / v2c12
    rot(4,3) = 0.25_dp * sqrt6 * (5.0_dp * w2 - 1.0_dp) * vv
    rot(5,3) = -2.5_dp * uu * vv * ww / v2c12
    rot(6,3) = sqrt10 * vv * (0.25_dp * (3.0_dp * w2 + 1.0_dp) &
        & -0.5_dp * w2 / v2c)
    rot(7,3) = -0.5_dp * sqrt15 * uu * vv * ww / v2c12

    rot(1,4) = sqrt10 * vv * ww * (-0.25_dp * (w2 + 3.0_dp) / v2c12 &
        & + w2 / v2c32)
    rot(2,4) = sqrt15 * uu * vv * w2 / v2c
    rot(3,4) = -0.25_dp * sqrt6 * vv * ww * (5.0_dp * w2 - 1.0_dp) / v2c12
    rot(4,4) = ww * (2.5_dp * w2 - 1.5_dp)
    rot(5,4) = -0.25_dp * sqrt6 * uu * (5.0_dp * w2 - 1.0_dp) / v2c12
    rot(6,4) = sqrt15 * ww * (0.5_dp * (w2 + 1.0_dp) - w2 / v2c)
    rot(7,4) = sqrt10 * uu * (-0.25_dp * (3.0_dp * w2 + 1.0_dp) / v2c12 &
        & + w2 / v2c32)

    rot(1,5) = sqrt15 * uu * vv * (-0.25_dp * (w2 + 1.0_dp) / v2c12 &
        & + w2 / v2c32)
    rot(2,5) = sqrt10 * vv * ww * (1.0_dp - 1.5_dp * w2 / v2c)
    rot(3,5) = -0.25_dp * uu * vv * (15.0_dp * w2 - 1.0_dp) / v2c12
    rot(4,5) = sqrt6 * uu * (1.25_dp * w2 - 0.25_dp)
    rot(5,5) = ww * (-2.5_dp * v2c12 &
        & - 0.25_dp * (-15.0_dp * w2 + 1.0_dp) / v2c12)
    rot(6,5) = sqrt10 * uu * (0.75_dp * w2 + 0.25_dp - 1.5_dp * w2 / v2c)
    rot(7,5) = sqrt15 * ww * (-0.5_dp * v2c12 &
        &+ (0.75_dp * w2 + 0.75_dp) / v2c12 - w2 / v2c32)

    rot(1,6) = sqrt6 * vv * ww * (-0.5_dp * v2c12 &
        & + (0.25_dp * w2 + 0.75_dp) / v2c12 - w2 / v2c32)
    rot(2,6) = vv * uu * (2.0_dp - 3.0_dp * w2 / v2c)
    rot(3,6) = sqrt10 * vv * ww * (-1.5_dp * v2c12 &
        & + (0.75_dp * w2 + 0.25_dp) / v2c12)
    rot(4,6) = sqrt15 * ww * (v2c - 0.5_dp * w2 - 0.5_dp)
    rot(5,6) = sqrt10 * uu * (-0.5_dp * v2c12 &
        & + (0.75_dp * w2 + 0.25_dp) / v2c12)
    rot(6,6) = ww * (3.0_dp * v2c - 1.5_dp * w2 - 3.5_dp + 3.0_dp * w2 / v2c)
    rot(7,6) = sqrt6 * uu * (-0.5_dp * v2c12 &
        & + (0.75_dp * w2 + 0.25_dp) / v2c12 - w2 / v2c32)

    rot(1,7) = uu * vv * (-v2c12 + (0.25_dp * w2 + 0.25_dp) / v2c12 &
        & - w2 / v2c32)
    rot(2,7) = sqrt6 * vv * ww * (-1.0_dp + 0.5_dp * w2 / v2c)
    rot(3,7) = sqrt15 * uu * vv * (-v2c12 + (0.25_dp * w2 + 0.25_dp) / v2c12)
    rot(4,7) = sqrt10 * uu * (v2c - 0.25_dp * w2 - 0.75_dp)
    rot(5,7) = sqrt15 * ww * (0.5_dp * v2c12 - (0.25_dp * w2 + 0.25_dp) / v2c12)
    rot(6,7) = sqrt6 * uu * (v2c - 0.25_dp * w2 - 0.75_dp + 0.5_dp * w2 / v2c)
    rot(7,7) = ww * (1.5_dp * v2c12 - (0.75_dp * w2 + 0.75_dp) / v2c12 &
        & + w2 / v2c32)

  end subroutine rotXY_f


end module sknew
