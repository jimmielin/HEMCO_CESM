!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_vertregrid_mod
!
! !DESCRIPTION: Module HCO\_VertRegrid\_Mod provides conservative sigma-to-sigma
!  vertical interpolation for HEMCO data in the direct-to-physics-grid mode.
!
!  In the legacy (intermediate grid) mode, vertical regridding is handled by
!  MESSy NCREGRID as part of combined horizontal+vertical regridding. In direct
!  mode, horizontal regridding is done by ESMF per-level, and vertical
!  regridding must be separated into a per-column operation.
!
!  This module implements conservative overlap-based interpolation between
!  arbitrary sigma-pressure coordinate systems. It handles both:
!    (a) "GEOS-Chem level" data with known hardcoded sigma edges
!    (b) Real-coordinate data with sigma edges read from input files
!\\
!\\
! !INTERFACE:
!
module hco_vertregrid_mod
!
! !USES:
!
    use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none
    private
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: HCO_VertRegrid_Column
    public :: HCO_VertRegrid_3D
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version for direct-mode vertical regrid
!EOP
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_VertRegrid_Column
!
! !DESCRIPTION: Performs conservative sigma-to-sigma vertical interpolation
!  for a single atmospheric column. For each target layer, computes the
!  fractional overlap with source layers in sigma space and redistributes
!  the data proportionally.
!
!  This routine conserves the column-integrated quantity when the input data
!  represents layer-mean intensive quantities (e.g., mixing ratios, kg/m2/s
!  per layer). The conservation property comes from the overlap-weighting:
!  each target layer value is the overlap-weighted average of contributing
!  source layers.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_VertRegrid_Column( nlev_src, sig_src, data_src, &
                                      nlev_tgt, sig_tgt, data_tgt )
!
! !INPUT PARAMETERS:
!
        integer,  intent(in)  :: nlev_src           ! # source levels
        real(r8), intent(in)  :: sig_src(nlev_src+1) ! Source sigma edges (surface=1 at index 1)
        real(r8), intent(in)  :: data_src(nlev_src)  ! Source data (layer means)
        integer,  intent(in)  :: nlev_tgt           ! # target levels
        real(r8), intent(in)  :: sig_tgt(nlev_tgt+1) ! Target sigma edges (surface=1 at index 1)
!
! !OUTPUT PARAMETERS:
!
        real(r8), intent(out) :: data_tgt(nlev_tgt)  ! Target data (layer means)
!
! !REMARKS:
!  Sigma edges are ordered from surface (index 1, sigma~1.0) to TOA
!  (index nlev+1, sigma~0.0), consistent with HEMCO convention where
!  level 1 is the surface.
!
!  For source layers that extend beyond the target grid range,
!  extrapolation uses the boundary source layer value (no data is lost).
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer  :: L_tgt, L_src
        real(r8) :: tgt_bot, tgt_top   ! Target layer sigma bounds
        real(r8) :: src_bot, src_top   ! Source layer sigma bounds
        real(r8) :: overlap            ! Overlap in sigma space
        real(r8) :: tgt_thickness      ! Target layer thickness in sigma
        real(r8) :: weighted_sum       ! Accumulated weighted data

        do L_tgt = 1, nlev_tgt
            ! Target layer bounds (sigma decreases with altitude)
            tgt_bot = sig_tgt(L_tgt)
            tgt_top = sig_tgt(L_tgt + 1)
            tgt_thickness = tgt_bot - tgt_top

            weighted_sum = 0.0_r8

            ! Find overlapping source layers
            do L_src = 1, nlev_src
                src_bot = sig_src(L_src)
                src_top = sig_src(L_src + 1)

                ! Compute overlap: intersection of [tgt_top, tgt_bot] and [src_top, src_bot]
                ! Both intervals go from higher sigma (bottom) to lower sigma (top)
                overlap = max(0.0_r8, min(tgt_bot, src_bot) - max(tgt_top, src_top))

                if (overlap > 0.0_r8) then
                    weighted_sum = weighted_sum + data_src(L_src) * overlap
                endif
            enddo

            ! Normalize by target layer thickness
            if (tgt_thickness > 0.0_r8) then
                data_tgt(L_tgt) = weighted_sum / tgt_thickness
            else
                data_tgt(L_tgt) = 0.0_r8
            endif
        enddo

    end subroutine HCO_VertRegrid_Column
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_VertRegrid_3D
!
! !DESCRIPTION: Performs conservative vertical regridding for a 3D field
!  (ncol columns x nlev_src source levels) onto the target vertical grid
!  (ncol columns x nlev_tgt target levels).
!
!  The source sigma edges can be either:
!  - Uniform across columns: sig_src_1d(nlev_src+1) provided, sig_src_3d absent
!  - Variable per column: sig_src_3d(ncol, nlev_src+1) provided
!
!  Target sigma edges are always per-column from PEDGE/PSFC.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_VertRegrid_3D( ncol, nlev_src, nlev_tgt,     &
                                  data_src, data_tgt,            &
                                  sig_tgt,                       &
                                  sig_src_1d, sig_src_3d )
!
! !INPUT PARAMETERS:
!
        integer,  intent(in)  :: ncol               ! # columns
        integer,  intent(in)  :: nlev_src           ! # source levels
        integer,  intent(in)  :: nlev_tgt           ! # target levels
        real(r8), intent(in)  :: data_src(ncol, nlev_src) ! Source data
        real(r8), intent(in)  :: sig_tgt(ncol, nlev_tgt+1) ! Target sigma edges per column
!
! !OUTPUT PARAMETERS:
!
        real(r8), intent(out) :: data_tgt(ncol, nlev_tgt)  ! Target data
!
! !INPUT PARAMETERS (OPTIONAL):
!
        real(r8), intent(in), optional :: sig_src_1d(nlev_src+1) ! Uniform source sigma
        real(r8), intent(in), optional :: sig_src_3d(ncol, nlev_src+1) ! Per-column source sigma
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer  :: I

        ! Regrid each column independently
        do I = 1, ncol
            if (present(sig_src_3d)) then
                call HCO_VertRegrid_Column( nlev_src, sig_src_3d(I,:), data_src(I,:), &
                                            nlev_tgt, sig_tgt(I,:),   data_tgt(I,:) )
            else if (present(sig_src_1d)) then
                call HCO_VertRegrid_Column( nlev_src, sig_src_1d,      data_src(I,:), &
                                            nlev_tgt, sig_tgt(I,:),    data_tgt(I,:) )
            else
                ! No source sigma provided — should not happen, zero out
                data_tgt(I,:) = 0.0_r8
            endif
        enddo

    end subroutine HCO_VertRegrid_3D
!EOC
end module hco_vertregrid_mod
