!!* Routines implementing transformation matrices for Hamiltonians as well as
!!* overlaps (and its moments) from the integration frame to the actual system.
module Transformation
#:include 'common.fypp'
  use accuracy
  use CommonTypes
  use Message
  implicit none
  private

  public :: rotate

  integer, parameter :: mAng = 2

  !!* Interface for rotation of matrices of different ranks
  interface rotate
    module procedure ij_rotate
    module procedure ijk_rotate
  end interface rotate

contains

  !!* Basic handler routine for 2-dimensional matrix rotations. Calls the
  !!* appropriate rotation prodecures depending on the angular momentum.
  !!* @param mat H or S matrix to fill
  !!* @param dirCos Directional cosine vector for rotation
  subroutine ij_rotate(mat, dirCos)
    real(dp), intent(out) :: mat(:,:)
    real(dp), intent(in) :: dirCos(3)

    real(dp), allocatable :: rot1(:,:), rot2(:,:)
    integer :: l1, l2

    @:ASSERT(size(dirCos) == 3)
    @:ASSERT(sqrt(sum(dirCos**2)) - 1.0_dp <= epsilon(0.0_dp))

    !! If the directional cosine vector is (0, 0, 1) we do not need to rotate.
    if (abs(dirCos(3)) - 1.0_dp <= -epsilon(0.0_dp)) then
      !! We derive the angular momenta from the size of the rotation matrices
      l1 = (size(mat, dim=1) - 1) / 2
      l2 = (size(mat, dim=2) - 1) / 2
      !! Size of the rotation matrices is the total number of possible magnetic
      !! quantum numbers for the given angular momentum, so 2*ll + 1.
      allocate(rot1(2*l1 + 1, 2*l1 + 1))
      allocate(rot2(2*l2 + 1, 2*l2 + 1))
      call getRotMat(rot1, dirCos, l1)
      call getRotMat(rot2, dirCos, l2)

      !! Sanity check for the matrices to be compatible.
      @:ASSERT(size(mat, dim=1) == size(rot1, dim=2))
      @:ASSERT(size(mat, dim=2) == size(rot2, dim=1))

      call rotateHS(mat, transpose(rot1), rot2) ! Execute rotation
    else if (abs(dirCos(3)) - 1.0_dp >= -epsilon(0.0_dp) .and. &
        &dirCos(3) < 0.0_dp) then
      mat(:,:) = -mat(:,:)
    end if
    !! No else condition: If we have a z unity vector (abs(dircos(3)) == 1),
    !! then we would need to rotate with the unity matrix (means: no rotation).

  end subroutine ij_rotate



  !!* Basic handler routine for 2-dimensional matrix rotations. Calls the
  !!* appropriate rotation prodecures depending on the angular momentum.
  !!* @param mat H or S matrix to fill
  !!* @param dirCos Directional cosine vector for rotation
  subroutine ijk_rotate(mat, dirCos)
    real(dp), intent(out) :: mat(:,:,-1:)
    real(dp), intent(in) :: dirCos(3)

    real(dp), allocatable :: rot1(:,:), rot2(:,:)
    integer :: l1, l2

    @:ASSERT(size(dirCos) == 3)
    @:ASSERT(sqrt(sum(dirCos**2)) - 1.0_dp <= epsilon(0.0_dp))
    @:ASSERT(mod(size(mat, dim=1) - 1, 2) == 0)
    @:ASSERT(mod(size(mat, dim=2) - 1, 2) == 0)

    !! If the directional cosine vector is (0, 0, 1) we do not need to rotate.
    if (abs(dirCos(3)) - 1.0_dp <= -epsilon(1.0_dp)) then
      !! We derive the angular momenta from the size of the rotation matrices
      l1 = (size(mat, dim=1) - 1) / 2
      l2 = (size(mat, dim=2) - 1) / 2
      !! Size of the rotation matrices is the total number of possible magnetic
      !! quantum numbers for the given angular momentum, so 2*ll + 1.
      allocate(rot1(2*l1 + 1, 2*l1 + 1))
      allocate(rot2(2*l2 + 1, 2*l2 + 1))
      call getRotMat(rot1, dirCos, l1)
      call getRotMat(rot2, dirCos, l2)

      !! Sanity check for the matrices to be compatible.
      @:ASSERT(size(mat, dim=1) == size(rot1, dim=2))
      @:ASSERT(size(mat, dim=2) == size(rot2, dim=1))

      call rotateFM(mat, transpose(rot1), rot2, dirCos) ! Execute rotation
    end if

  end subroutine ijk_rotate



  !!* Rotation matrix handler routine to deliver the right rotation matrices
  !!* for a given angular momentum.
  !!* @param rot The rotation matrix to be filled
  !!* @param dirCos Directional cosine vector for rotation
  !!* @param ll Angular momentum of the corresponding tesseral
  subroutine getRotMat(rot, dirCos, ll)
    real(dp), intent(out) :: rot(:,:)
    integer, intent(in) :: ll
    real(dp), intent(in) :: dirCos(3)

    character(lc) :: errorStr

    select case(ll)
    case(0)
      call getRotMatS(rot)
    case(1)
      call getRotMatP(rot, dirCos)
    case(2)
      call getRotMatD(rot, dirCos)
    case default
      write (errorStr, "(A,I2,A,I2)") "Incompatible angular momentum ", ll, &
          &"for matrix rotations. Maximum supported l = ", mAng
      call error(errorStr)
    end select

  end subroutine getRotMat



  !!* Rotates Hamiltonian or overlap matrices.
  !!* @param mat The matrix to be rotated
  !!* @param rot1 The first (transposed) rotation matrix
  !!* @param rot2 The second rotation matrix
  subroutine rotateHS(mat, rot1, rot2)
    real(dp), intent(inout) :: mat(:,:)
    real(dp), intent(in) :: rot1(:,:), rot2(:,:)

    mat = matmul(rot1, mat)
    mat = matmul(mat, rot2)

  end subroutine rotateHS



  !!* Rotates Hamiltonian or first moment matrices.
  !!* @param mat The matrix to be rotated
  !!* @param rot1 The first (transposed) rotation matrix
  !!* @param rot2 The second rotation matrix
  !!* @param dirCos Directional cosine vector
  subroutine rotateFM(mat, rot1, rot2, dirCos)
    real(dp), intent(inout) :: mat(:,:,-1:)
    real(dp), intent(in) :: rot1(:,:), rot2(:,:)
    real(dp), intent(in) :: dirCos(3)

    !! Tesseral decomposition rotation matrix for the first moment:
    real(dp), allocatable :: rotX1(:,:) ! Rotation matrix for m
    real(dp), allocatable :: tmpMat(:,:,:) ! Temporary matrix
    integer :: iMM ! m quantum number of the first moment matrix
    integer :: iMMP ! m prime quantum number (for summation)

    @:ASSERT(size(mat, dim=3) == 3)

    allocate(rotX1(-1:1, -1:1))
    allocate(tmpMat(size(mat, dim=1), size(mat, dim=2), -1:1))

    call getRotMat(rotX1, dirCos, 1)

    !! Rotation of the orbitals portion of the matrix
    do iMM = -1, 1
      mat(:,:,iMM) = matmul(rot1, mat(:,:,iMM))
      tmpMat(:,:,iMM) = matmul(mat(:,:,iMM), rot2)
    end do

    !! Rotation of the dipole vector portion of the matrix
    mat(:,:,:) = 0.0_dp ! Set zero so the output can be put in there
    do iMM = -1, 1
      do iMMP = -1, 1
        mat(:,:,iMM) = mat(:,:,iMM) + rotX1(iMMP,iMM) * &
            tmpMat(:,:,iMMP)
      end do
    end do

  end subroutine rotateFM



  !!* Rotates Hamiltonian or second moment matrices.
  !!* @param mat The matrix to be rotated
  !!* @param rot1 The first (transposed) rotation matrix
  !!* @param rot2 The second rotation matrix
  !!* @param dirCos Directional cosine vector
  !!* @note: THIS IS NOT THE REAL SUBROUTINE YET, JUST A DUMMY
  subroutine rotateSM(mat, rot1, rot2, dirCos)
    real(dp), intent(inout) :: mat(:,:)
    real(dp), intent(in) :: rot1(:,:), rot2(:,:)
    real(dp), intent(in) :: dirCos(3)

    !! Tesseral decomposition rotation matrices of the second moment:
    real(dp), allocatable :: rotQ0(:,:) ! For Q^0 elements
    real(dp), allocatable :: rotQ2(:,:) ! For Q^2 elements

    allocate(rotQ0(size(mat, dim=1), size(mat, dim=2)))
    allocate(rotQ2(size(mat, dim=1), size(mat, dim=2)))

    call getRotMat(rotQ0, dirCos, 0)
    call getRotMat(rotQ2, dirCos, 2)

    mat = matmul(rot1, mat)
    mat = matmul(mat, rot2)

  end subroutine rotateSM



  !!* Fills the given array with the rotation matrix for s-Orbital tesserals.
  !!* @param rot Rotation matrix
  !!* @param dirCos Directional cosine vector for rotation
  subroutine getRotMatS(rot)
    real(dp), intent(out) :: rot(:,:)

    @:ASSERT(all(shape(rot) == (/ 1, 1 /)))

    rot(1,1) = 1

  end subroutine getRotMatS



  !!* Fills the given array with the rotation matrix for p-Orbital tesserals.
  !!* @param rot Rotation matrix
  !!* @param dirCos Directional cosine vector for rotation
  subroutine getRotMatP(rot, dirCos)
    real(dp), intent(out) :: rot(:,:)
    real(dp), intent(in) :: dirCos(3)

    real(dp) :: tmpR

    @:ASSERT(all(shape(rot) == (/ 3, 3 /)))
    @:ASSERT(size(dirCos) == 3)

    !! Aliases for directional cosines w/o pointer indirection
    associate (uu => dirCos(1), vv => dirCos(2), ww => dirCos(3))

      tmpR = sqrt(1.0_dp - ww**2)

      rot(1,1) = uu / tmpR
      rot(2,1) = vv
      rot(3,1) = ww * vv / tmpR
      rot(1,2) = 0.0_dp
      rot(2,2) = ww
      rot(3,2) = - tmpR
      rot(1,3) = - vv  / tmpR
      rot(2,3) = uu
      rot(3,3) = ww * uu / tmpR

    end associate

  end subroutine getRotMatP



  !!* Fills the given array with the rotation matrix for d-Orbital tesserals.
  !!* @param rot Rotation matrix
  !!* @param dirCos Directional cosine vector for rotation
  !!* @note THIS ROUTINE HASN'T BEEN TESTED YET!
  subroutine getRotMatD(rot, dirCos)
    real(dp), intent(out) :: rot(:,:)
    real(dp), intent(in) :: dirCos(3)

    real(dp) :: tmpR1, tmpR2, tmpR3

    @:ASSERT(all(shape(rot) == (/ 5, 5 /)))
    @:ASSERT(size(dirCos) == 3)

    !! Aliases for directional cosines w/o pointer indirection
    associate (uu => dirCos(1), vv => dirCos(2), ww => dirCos(3))

      tmpR1 = sqrt(1.0_dp-ww**2)
      tmpR2 = ww**2 - 1
      tmpR3 = 2.0_dp * uu**2 + ww**2 - 1.0_dp

      !! 1.73205080756887729353_dp = sqrt(3)
      !! 0.86602540378443864676_dp = 0.5 * sqrt(3)
      rot(1,1) = -tmpR3 * ww / tmpR2
      rot(2,1) = tmpR3 / tmpR1
      rot(3,1) = vv * uu * 1.73205080756887729353_dp
      rot(4,1) = 2.0_dp * ww * vv * uu / tmpR1
      rot(5,1) = -(ww**2 + 1.0_dp) * vv * uu / tmpR2
      rot(1,2) = -uu
      rot(2,2) = ww * uu / tmpR1
      rot(3,2) = ww * vv * 1.73205080756887729353_dp
      rot(4,2) = (2.0_dp * ww**2 - 1) * vv / tmpR1
      rot(5,2) = -ww * vv
      rot(1,3) = 0.0_dp
      rot(2,3) = 0.0_dp
      rot(3,3) = 1.5_dp * ww**2 - 0.5_dp
      rot(4,3) = -tmpR1 * ww * 1.73205080756887729353_dp
      rot(5,3) = -tmpR2 * 0.86602540378443864676_dp
      rot(1,4) = vv
      rot(2,4) = -ww * vv / tmpR1
      rot(3,4) = ww * uu * 1.73205080756887729353_dp
      rot(4,4) = (2.0_dp * ww**2 - 1.0_dp) * uu / tmpR1
      rot(5,4) = -ww * uu
      rot(1,5) = 2.0_dp * ww * vv * uu / tmpR2
      rot(2,5) = -2.0_dp * vv * uu / tmpR1
      rot(3,5) = tmpR3 * 0.86602540378443864676_dp
      rot(4,5) = tmpR3 * ww / tmpR1
      rot(5,5) = -(2.0_dp * ww**2 * uu**2 + ww**4 + 2.0_dp * uu**2 - 1.0_dp) / &
          &(2.0_dp * tmpR1)

    end associate

  end subroutine getRotMatD

end module Transformation
