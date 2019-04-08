!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains Types and subroutine to build up and query a Slater-Koster table where the integrals are
!> specified on an equidistant grid.
module slakoeqgrid
  use assert
  use accuracy
  use interpolation
  use message
  implicit none
  private

  public :: OSlakoEqGrid, OOrigins, init
  public :: getSKIntegrals, getNIntegrals, getCutoff
  public :: skEqGridOld, skEqGridNew, getMPOrigins


  !> Represents an equally spaced Slater-Koster grid
  type OSlakoEqGrid
    private
    integer :: nGrid
    integer :: nInteg
    real(dp) :: dist
    real(dp), allocatable :: skTab(:,:)
    integer :: skIntMethod
    logical :: tInit = .false.
  end type OSlakoEqGrid

  !!* Origin vectors of SlakoEqGrid integrals for the multipole expansion
  !!* @note I don't really know why I put that here, I should move this some
  !!* other place at some point.
  type OOrigins
    private
    integer :: nInteg ! Number of integrals of the multipole table
    integer :: dims ! Number of dimensions of each origin vector (should be 3)
    real(dp), pointer :: origin(:,:) ! Table origin(dims, nInteg)
    integer, pointer :: center(:) ! Center table
    logical :: tInit = .false.
  end type OOrigins

  !!* On-site Integrals for the multipole expansion
  !!* @note I don't really know why I put that here, I should move this some
  !!* other place at some point.
  ! type OOnSites
  !   private
  !   integer :: nInteg ! Number of on-site integrals in the multipole table
  !   real(dp), pointer :: onSites(:) ! On-Site integrals
  !   logical :: tInit = .false.
  ! end type OOnSites


  !> Initialises SlakoEqGrid.
  interface init
    module procedure SlakoEqGrid_init
    module procedure Origins_init
    ! module procedure OnSites_init
  end interface init


  !> Returns the integrals for a given distance.
  interface getSKIntegrals
    module procedure SlakoEqGrid_getSKIntegrals
  end interface getSKIntegrals


  !> Returns the number of integrals the table contains
  interface getNIntegrals
    module procedure SlakoEqGrid_getNIntegrals
    module procedure Origins_getNIntegrals
  end interface getNIntegrals

  !!* Returns multipole origins for a given distance
  interface getMPOrigins
    module procedure Origins_getMPOrigins
  end interface getMPOrigins


  !> Returns the cutoff of the interaction.
  interface getCutoff
    module procedure SlakoEqGrid_getCutoff
  end interface getCutoff

  ! Interpolation methods

  !> Historical method
  integer, parameter :: skEqGridOld = 1

  !> Current method
  integer, parameter :: skEqGridNew = 2

  ! Nr. of grid points to use for the polynomial interpolation

  !> Historical choice
  integer, parameter :: nInterOld_ = 3

  !> Present choice
  integer, parameter :: nInterNew_ = 8

  ! Nr. of grid points on the right of the interpolated point.

  ! For an odd number of intervals, the number of right points should be bigger than the number of
  ! left points, to remain compatible with the old code.


  !> value nRightInterOld: floor(real(nInterOld_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterOld_ = 2

  !> value nRightInterNew: floor(real(nInterNew_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterNew_ = 4


  !> Displacement for deriving interpolated polynomials
  real(dp), parameter :: deltaR_ = 1e-5_dp

contains


  !> Initialises SlakoEqGrid.
  subroutine SlakoEqGrid_init(self, dist, table, skIntMethod)

    !> SlakoEqGrid instance.
    type(OSlakoEqGrid), intent(out) :: self

    !> Distance between the grid points.
    real(dp), intent(in) :: dist

    !> Slater-Koster table (first entry belongs to first grid point)
    real(dp), intent(in) :: table(:,:)

    !> Method for the interpolation between the entries.
    integer, intent(in) :: skintMethod

    @:ASSERT(.not. self%tInit)
    @:ASSERT(dist >= 0.0_dp)
    @:ASSERT(skIntMethod == skEqGridOld .or. skIntMethod == skEqGridNew)

    self%dist = dist
    self%nGrid = size(table, dim=1)
    self%nInteg = size(table, dim=2)
    allocate(self%skTab(self%nGrid, self%nInteg))
    self%skTab(:,:) = table(:,:)
    self%skIntMethod = skIntMethod
    self%tInit = .true.

  end subroutine SlakoEqGrid_init


  !!* Initialises Origin Vectors.
  !!* @param self Origins instance.
  !!* @param mpoleOrigins Origin vector table.
  !!* @param mpoleCenters Multipole centers table
  subroutine Origins_init(self, mpoleOrigins, mpoleCenters)
    type(OOrigins), intent(out) :: self
    real(dp), intent(in) :: mpoleOrigins(:,:)
    integer, intent(in) :: mpoleCenters(:)

    @:ASSERT(.not. self%tInit)
    @:ASSERT(size(mpoleOrigins, dim=2) == size(mpoleCenters))

    self%nInteg = size(mpoleCenters)
    self%dims = size(mpoleOrigins, dim=1)
    allocate(self%origin(self%dims, self%nInteg))
    allocate(self%center(self%nInteg))
    self%origin(:,:) = mpoleOrigins(:,:)
    self%center(:) = mpoleCenters(:)
    self%tInit = .true.

  end subroutine Origins_init



  !!* Initialises on-site integrals
  !!* @param self On-site integral instance.
  !!* @param table Table containing the integrals
  ! subroutine OnSites_init(self, table)
  !   type(OOnSites), intent(out) :: self
  !   real(dp), intent(in) :: table(:)

  !   @:ASSERT(.not. self%tInit)

  !   self%nInteg = size(table, dim=2)
  !   INITALLOCATE_PARR(self%onSites, (self%nInteg))
  !   self%onSites(:) = table(:)
  !   self%tInit = .true.

  ! end subroutine OnSites_init




  !> Returns the integrals for a given distance.
  subroutine SlakoEqGrid_getSKIntegrals(self, sk, dist)

    !> SlakoEqGrid instance.
    type(OSlakoEqGrid), intent(in) :: self

    !> Contains the interpolated integrals on exit
    real(dp), intent(out) :: sk(:)

    !> Distance for which the integrals should be interpolated.
    real(dp), intent(in) :: dist

    @:ASSERT(self%tInit)
    @:ASSERT(size(sk) >= self%nInteg)
    @:ASSERT(dist >= 0.0_dp)

    if (self%skIntMethod == skEqGridOld) then
      call SlakoEqGrid_interOld_(self, sk, dist)
    else
      call SlakoEqGrid_interNew_(self, sk, dist)
    end if

  end subroutine SlakoEqGrid_getSKIntegrals



  !!* Returns the multipole origin vectors for a given distance.
  !!* @param self SlakoEqGrid instance.
  !!* @param origVecs Contains the interpolated integrals on exit
  !!* @param dist Distance for which the integrals should be interpolated.
  subroutine Origins_getMPOrigins(self, origVecs, centers, dist)
    type(OOrigins), intent(in) :: self
    real(dp), intent(out) :: origVecs(:,:)
    integer, intent(out) :: centers(:)
    real(dp), intent(in) :: dist

    @:ASSERT(self%tInit)
    @:ASSERT(self%dims == 3)
    @:ASSERT(size(origVecs, dim=2) >= self%nInteg)
    @:ASSERT(size(centers) >= self%nInteg)
    @:ASSERT(dist >= 0.0_dp)

    origVecs(:,1:size(self%origin, dim=2)) = self%origin(:,:)
    centers(1:size(self%center)) = self%center(:)

  end subroutine Origins_getMPOrigins



  !!* Returns the on-site integrals.
  !!* @param self SlakoEqGrid instance.
  !!* @param onSites Contains the on-site integrals on exit
  ! subroutine OnSites_getIntegrals(self, onSites)
  !   type(OOrigins), intent(in) :: self
  !   real(dp), intent(out) :: onSites(:)

  !   @:ASSERT(self%tInit)
  !   @:ASSERT(size(origVecs) >= self%nInteg)

  !   onSites(:) = self%onSites(:)

  ! end subroutine OnSites_getIntegrals



  !> Returns the number of intgrals the table contains
  function SlakoEqGrid_getNIntegrals(self) result(nInt)

    !> SlakoEqGrid instance.
    type(OSlakoEqGrid), intent(in) :: self

    !> Number of integrals.
    integer :: nInt

    nInt = self%nInteg

  end function SlakoEqGrid_getNIntegrals

  !!* Returns the number of integrals the table contains
  !!* @param self Origins instance.
  !!* @return Number of integrals.
  function Origins_getNIntegrals(self) result(nInt)
    type(OOrigins), intent(in) :: self
    integer :: nInt

    nInt = self%nInteg

  end function Origins_getNIntegrals




  !> Returns the cutoff of the interaction.
  function SlakoEqGrid_getCutoff(self) result(cutoff)

    !>  SlakoEqGrid instance.
    type(OSlakoEqGrid), intent(in) :: self

    !> grid cutoff
    real(dp) :: cutoff

    cutoff = real(self%nGrid, dp) * self%dist
    if (self%skIntMethod == skEqGridOld) then
      cutoff = cutoff + distFudgeOld
    else
      cutoff = cutoff + distFudge
    end if

  end function SlakoEqGrid_getCutoff


  !> Inter- and extrapolation for SK-tables, new method.
  subroutine SlakoEqGrid_interNew_(self, dd, rr)

    !> SlakoEqGrid table on equiv. grid
    type(OSlakoEqGrid), intent(in) :: self

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> distance bewteen two atoms of interest
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterNew_), ya(nInterNew_), yb(self%nInteg,nInterNew_), y1, y1p, y1pp
    real(dp) :: incr, dr, rMax, y0(self%nInteg), y2(self%nInteg)
    integer :: leng, ind, iLast
    integer :: ii

    leng = self%nGrid
    incr = self%dist
    rMax = real(leng, dp) * incr + distFudge
    ind = floor(rr / incr)

    !! Sanity check, if SK-table contains enough entries
    if (leng < nInterNew_ + 1) then
      call error("SlakoEqGrid: Not enough points in the SK-table for &
          &interpolation!")
    end if

    dd(:) = 0.0_dp
    if (rr >= rMax) then
      !! Beyond last grid point + distFudge => no interaction
      dd(:) = 0.0_dp
    elseif (ind < leng) then
      !! Closer to origin than last grid point => polynomial fit
      iLast = min(leng, ind + nRightInterNew_)
      iLast = max(iLast, nInterNew_)
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(self%skTab(iLast-nInterNew_+1:iLast,:self%nInteg))
      dd(:self%nInteg) = polyInterUniform(xa, yb, rr)
    else
      !! Beyond the grid => extrapolation with polynomial of 5th order
      dr = rr - rMax
      iLast = leng
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(self%skTab(iLast-nInterNew_+1:iLast,:self%nInteg))
      y0 = polyInterUniform(xa, yb, xa(nInterNew_) - deltaR_)
      y2 = polyInterUniform(xa, yb, xa(nInterNew_) + deltaR_)
      do ii = 1, self%nInteg
        ya(:) = self%skTab(iLast-nInterNew_+1:iLast, ii)
        y1 = ya(nInterNew_)
        y1p = (y2(ii) - y0(ii)) / (2.0_dp * deltaR_)
        y1pp = (y2(ii) + y0(ii) - 2.0_dp * y1) / (deltaR_ * deltaR_)
        dd(ii) = poly5ToZero(y1, y1p, y1pp, dr, -1.0_dp * distFudge)
      end do
    end if

  end subroutine SlakoEqGrid_interNew_


  !> Inter- and extra-polation for SK-tables equivalent to the old DFTB code.
  subroutine SlakoEqGrid_interOld_(self, dd, rr)

    !> Data structure for SK interpolation
    type(OSlakoEqGrid), intent(in) :: self

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> distance bewteen two atoms of interest
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterOld_), yb(self%nInteg,nInterOld_),y0, y1, y2, y1p, y1pp
    real(dp) :: incr, dr
    integer :: leng, ind, mInd, iLast
    integer :: ii
    real(dp) :: r1, r2

    leng = self%nGrid
    incr = self%dist
    mInd = leng + floor(distFudgeOld/incr)
    ind = floor(rr / incr)

    !! Sanity check, if SK-table contains enough entries
    if (leng < nInterOld_ + 1) then
      call error("skspar: Not enough points in the SK-table for interpolation!")
    end if

    dd(:) = 0.0_dp
    if (ind < leng-1) then
      !! Distance closer than penultimate grid point => polynomial fit
      iLast = min(leng, ind + nRightInterOld_)
      iLast = max(iLast, nInterOld_)
      do ii = 1, nInterOld_
        xa(ii) = real(iLast - nInterOld_ + ii, dp) * incr
      end do
      yb = transpose(self%skTab(iLast-nInterOld_+1:iLast,:self%nInteg))
      dd(:self%nInteg) = polyInterUniform(xa, yb, rr)
    elseif (ind < leng) then
      !! Distance between penultimate and last grid point => free cubic spline
      dr = rr - real(leng - 1, dp) * incr
      do ii = 1, self%nInteg
        y0 = self%skTab(leng-2, ii)
        y1 = self%skTab(leng-1, ii)
        y2 = self%skTab(leng, ii)
        y1p = (y2 - y0) / (2.0_dp * incr)
        y1pp = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, y1p, y1pp, incr, y2, dr, dd(ii))
      end do
    elseif (ind < mInd - 1) then
      !! Extrapolation
      dr = rr - real(mInd - 1, dp) * incr
      do ii = 1, self%nInteg
        y0 = self%skTab(leng-2, ii)
        y1 = self%skTab(leng-1, ii)
        y2 = self%skTab(leng, ii)
        r1 = (y2 - y0) / (2.0_dp * incr)
        r2 = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, r1, r2, incr, y2, incr, yp=y1p, ypp=y1pp)
        dd(ii) = poly5ToZero(y2, y1p, y1pp, dr, &
            &-1.0_dp * real(mInd - leng -1, dp)*incr)
      end do
    else
      !! Dist. greater than tabulated sk range + distFudge => no interaction
      dd(:) = 0.0_dp
    end if

  end subroutine SlakoEqGrid_interOld_

end module slakoeqgrid
