#define VERIFY_(A) if(.not.HCO_ESMF_VRFY(A,subname,__LINE__)) stop -1
#define ASSERT_(A) if(.not.HCO_ESMF_ASRT(A,subname,__LINE__)) stop -1
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_esmf_regrid_cache
!
! !DESCRIPTION: Module HCO\_ESMF\_Regrid\_Cache provides ESMF-based regridding
!  infrastructure for the HEMCO direct-to-physics-grid mode. It manages a cache
!  of ESMF route handles for regridding from various input file grids
!  (rectilinear lat-lon) directly to the CAM physics mesh.
!
!  In the legacy (intermediate grid) mode, HEMCO reads data and regrids it
!  via MAP_A2A or MESSy NCREGRID onto a uniform rectilinear intermediate grid.
!  In direct mode, this module replaces that regridding step with ESMF
!  conservative regridding that goes directly from each input file's native
!  lat-lon grid to the CAM physics mesh (which may be unstructured).
!
!  Route handles are cached per unique input grid to avoid the cost of
!  repeated ESMF_FieldRegridStore calls. Typical HEMCO configurations use
!  5-10 unique input grids, so the cache is small.
!
!  This module also provides the HCO\_ESMF\_REGRID\_DIRECT subroutine, which
!  is called from hcoio\_read\_pio\_mod.F90 in place of MAP_A2A/MESSy when
!  direct mode is enabled.
!\\
!\\
! !INTERFACE:
!
module hco_esmf_regrid_cache
!
! !USES:
!
    use hco_esmf_wrappers
    use ESMF,            only: ESMF_Mesh, ESMF_Grid, ESMF_Field, ESMF_RouteHandle
    use ESMF,            only: ESMF_SUCCESS
    use shr_kind_mod,    only: r8 => shr_kind_r8
    use HCO_Types_Mod,   only: ListCont, hp, sp, dp
    use HCO_State_Mod,   only: HCO_State

    implicit none
    private
    save
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: HCO_RegridCache_Init
    public :: HCO_RegridCache_Cleanup
    public :: HCO_ESMF_REGRID_DIRECT
!
! !PUBLIC DATA:
!
    ! Flag indicating whether direct mode is active.
    ! Set by HEMCO_CESM during initialization, read by hcoio_read_pio_mod.
    logical, public     :: HcoDirectMode = .false.
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version for direct-mode regridding
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    integer, parameter :: MAX_CACHED_GRIDS = 20

    type :: RegridCacheEntry
        integer  :: nlon = 0            ! Input grid # longitudes
        integer  :: nlat = 0            ! Input grid # latitudes
        real(r8) :: lon0 = -999.0_r8    ! First longitude edge (for uniqueness)
        real(r8) :: lat0 = -999.0_r8    ! First latitude edge (for uniqueness)
        type(ESMF_Grid)        :: srcGrid
        type(ESMF_Field)       :: srcField2D
        type(ESMF_Field)       :: dstField2D
        type(ESMF_RouteHandle) :: rh2D
        logical  :: initialized = .false.
    end type RegridCacheEntry

    type(RegridCacheEntry)     :: cache(MAX_CACHED_GRIDS)
    integer                    :: nCached = 0

    ! Reference to the CAM physics mesh (set during init)
    type(ESMF_Mesh)     :: phys_mesh
    integer             :: phys_ncol = 0    ! Local # physics columns

contains
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_RegridCache_Init
!
! !DESCRIPTION: Initializes the regrid cache with a reference to the CAM
!  physics mesh. Must be called after HCO_Grid_ESMF_CreateCAM.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_RegridCache_Init( mesh, ncol, RC )
!
! !INPUT PARAMETERS:
!
        type(ESMF_Mesh), intent(in) :: mesh
        integer,         intent(in) :: ncol  ! Local # physics columns
!
! !OUTPUT PARAMETERS:
!
        integer, intent(out)        :: RC
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
        character(len=*), parameter :: subname = 'HCO_RegridCache_Init'

        phys_mesh = mesh
        phys_ncol = ncol
        nCached   = 0
        RC        = ESMF_SUCCESS

    end subroutine HCO_RegridCache_Init
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_RegridCache_GetRH
!
! !DESCRIPTION: Looks up or creates an ESMF route handle for regridding from
!  an input rectilinear lat-lon grid (described by its edge arrays) to the
!  CAM physics mesh. On cache miss, creates the ESMF grid, fields, and
!  route handle, and stores them in the cache.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_RegridCache_GetRH( nlon, nlat, LonEdge, LatEdge, &
                                      idx, RC )
!
! !USES:
!
        use cam_logfile,  only: iulog
        use spmd_utils,   only: masterproc

        use ESMF,         only: ESMF_GridCreate1PeriDim, ESMF_INDEX_GLOBAL
        use ESMF,         only: ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
        use ESMF,         only: ESMF_GridAddCoord, ESMF_GridGetCoord
        use ESMF,         only: ESMF_TYPEKIND_R8, ESMF_KIND_R8
        use ESMF,         only: ESMF_MESHLOC_ELEMENT
        use ESMF,         only: ESMF_ArraySpec, ESMF_ArraySpecSet
        use ESMF,         only: ESMF_FieldCreate, ESMF_FieldRegridStore
        use ESMF,         only: ESMF_REGRIDMETHOD_CONSERVE
        use ESMF,         only: ESMF_POLEMETHOD_NONE
        use ESMF,         only: ESMF_RouteHandleDestroy, ESMF_FieldDestroy, ESMF_GridDestroy
!
! !INPUT PARAMETERS:
!
        integer,  intent(in)  :: nlon, nlat       ! Input grid dimensions
        real(hp), intent(in)  :: LonEdge(nlon+1)  ! Longitude edges [deg]
        real(hp), intent(in)  :: LatEdge(nlat+1)  ! Latitude edges [deg]
!
! !OUTPUT PARAMETERS:
!
        integer,  intent(out) :: idx               ! Cache index for this grid
        integer,  intent(out) :: RC
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
        character(len=*), parameter :: subname = 'HCO_RegridCache_GetRH'
        integer  :: n, i, j
        integer  :: lbnd(2), ubnd(2)
        real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
        real(ESMF_KIND_R8), pointer :: coordX_E(:,:), coordY_E(:,:)
        type(ESMF_ArraySpec) :: arrayspec
        real(r8) :: lon0_in, lat0_in
        real(r8) :: dx, dy

        RC = ESMF_SUCCESS

        ! Composite key for cache lookup
        lon0_in = real(LonEdge(1), r8)
        lat0_in = real(LatEdge(1), r8)

        ! Check cache for existing entry
        do n = 1, nCached
            if (cache(n)%initialized .and. &
                cache(n)%nlon == nlon .and. cache(n)%nlat == nlat .and. &
                abs(cache(n)%lon0 - lon0_in) < 1.0e-6_r8 .and. &
                abs(cache(n)%lat0 - lat0_in) < 1.0e-6_r8) then
                ! Cache hit
                idx = n
                return
            endif
        enddo

        ! Cache miss — create new entry
        if (nCached >= MAX_CACHED_GRIDS) then
            if (masterproc) then
                write(iulog,*) "HEMCO RegridCache: WARNING - cache full (", MAX_CACHED_GRIDS, &
                               " entries). Reusing last slot."
            endif
            nCached = MAX_CACHED_GRIDS

            ! Destroy existing ESMF objects in the slot being overwritten
            if (cache(nCached)%initialized) then
                call ESMF_RouteHandleDestroy(cache(nCached)%rh2D, rc=RC)
                call ESMF_FieldDestroy(cache(nCached)%srcField2D, rc=RC)
                call ESMF_FieldDestroy(cache(nCached)%dstField2D, rc=RC)
                call ESMF_GridDestroy(cache(nCached)%srcGrid, rc=RC)
                cache(nCached)%initialized = .false.
            endif
        else
            nCached = nCached + 1
        endif
        idx = nCached

        cache(idx)%nlon = nlon
        cache(idx)%nlat = nlat
        cache(idx)%lon0 = lon0_in
        cache(idx)%lat0 = lat0_in

        if (masterproc) then
            write(iulog,*) "HEMCO RegridCache: Creating route handle for input grid ", &
                           nlon, "x", nlat, " (entry ", idx, ")"
        endif

        !-----------------------------------------------------------------------
        ! Create ESMF Grid for the input file (rectilinear, single-PE for now)
        ! Each PE creates the full input grid — ESMF handles the decomposition
        ! internally for FieldRegridStore.
        !-----------------------------------------------------------------------

        ! Compute grid center coordinates from edges
        dx = real(LonEdge(2) - LonEdge(1), r8)
        dy = real(LatEdge(2) - LatEdge(1), r8)

        ! Create 1-periodic-dim rectilinear grid on a single DE per PET
        ! (each PET holds the full grid — simplest decomposition for source data
        !  that is read in full by PIO on each PE)
        cache(idx)%srcGrid = ESMF_GridCreate1PeriDim(        &
            maxIndex=(/nlon, nlat/),                          &
            indexflag=ESMF_INDEX_GLOBAL,                      &
            rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Add center and corner coordinates
        call ESMF_GridAddCoord(cache(idx)%srcGrid, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridAddCoord(cache(idx)%srcGrid, &
                               staggerloc=ESMF_STAGGERLOC_CORNER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Fill center coordinates
        call ESMF_GridGetCoord(cache(idx)%srcGrid, coordDim=1, localDE=0, &
                               computationalLBound=lbnd, computationalUBound=ubnd, &
                               farrayPtr=coordX, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridGetCoord(cache(idx)%srcGrid, coordDim=2, localDE=0, &
                               farrayPtr=coordY, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Centers are midpoints of edges
        do j = lbnd(2), ubnd(2)
            do i = lbnd(1), ubnd(1)
                coordX(i, j) = real(LonEdge(i) + LonEdge(i+1), r8) * 0.5_r8
                coordY(i, j) = real(LatEdge(j) + LatEdge(j+1), r8) * 0.5_r8
            enddo
        enddo

        ! Fill corner coordinates
        call ESMF_GridGetCoord(cache(idx)%srcGrid, coordDim=1, localDE=0, &
                               computationalLBound=lbnd, computationalUBound=ubnd, &
                               farrayPtr=coordX_E, &
                               staggerloc=ESMF_STAGGERLOC_CORNER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridGetCoord(cache(idx)%srcGrid, coordDim=2, localDE=0, &
                               farrayPtr=coordY_E, &
                               staggerloc=ESMF_STAGGERLOC_CORNER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        do j = lbnd(2), ubnd(2)
            do i = lbnd(1), ubnd(1)
                coordX_E(i, j) = real(LonEdge(min(i, nlon+1)), r8)
                coordY_E(i, j) = real(LatEdge(min(j, nlat+1)), r8)
            enddo
        enddo

        !-----------------------------------------------------------------------
        ! Create source and destination fields
        !-----------------------------------------------------------------------

        ! Source field: 2D on input grid
        call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        cache(idx)%srcField2D = ESMF_FieldCreate(cache(idx)%srcGrid, arrayspec, &
            name='HCO_INPUT_SRC_2D', &
            staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Destination field: 1D on CAM physics mesh
        call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        cache(idx)%dstField2D = ESMF_FieldCreate(phys_mesh, arrayspec, &
            name='HCO_INPUT_DST_2D', &
            meshloc=ESMF_MESHLOC_ELEMENT, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        !-----------------------------------------------------------------------
        ! Create route handle (conservative regridding)
        !-----------------------------------------------------------------------
        call ESMF_FieldRegridStore(                                       &
            srcField=cache(idx)%srcField2D,                                &
            dstField=cache(idx)%dstField2D,                                &
            regridMethod=ESMF_REGRIDMETHOD_CONSERVE,                       &
            poleMethod=ESMF_POLEMETHOD_NONE,                               &
            routeHandle=cache(idx)%rh2D,                                   &
            srcTermProcessing=0,                                           &
            pipelineDepth=16, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        cache(idx)%initialized = .true.

        if (masterproc) then
            write(iulog,*) "HEMCO RegridCache: Route handle created for ", &
                           nlon, "x", nlat, " -> physics mesh"
        endif

    end subroutine HCO_RegridCache_GetRH
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ESMF_REGRID_DIRECT
!
! !DESCRIPTION: Main entry point for direct ESMF regridding. Called from
!  hcoio\_read\_pio\_mod.F90 in place of REGRID\_MAPA2A and HCO\_MESSY\_REGRID
!  when direct mode is enabled.
!
!  For 2D data: performs ESMF conservative horizontal regridding directly
!  from the input file grid to physics columns.
!
!  For 3D data: performs ESMF horizontal regridding per-level, then
!  sigma-to-sigma conservative vertical interpolation per-column.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_ESMF_REGRID_DIRECT( HcoState, NcArr, LonEdge, LatEdge, &
                                       SigEdge, Lct, IsModelLevel, RC )
!
! !USES:
!
        use cam_logfile,         only: iulog
        use spmd_utils,          only: masterproc

        use ESMF,                only: ESMF_FieldRegrid, ESMF_FieldGet
        use ESMF,                only: ESMF_TERMORDER_SRCSEQ
        use ESMF,                only: ESMF_KIND_R8

        use HCO_FileData_Mod,    only: FileData_ArrCheck
        use HCO_Error_Mod,       only: HCO_SUCCESS

        use hco_vertregrid_mod,  only: HCO_VertRegrid_3D
!
! !INPUT PARAMETERS:
!
        type(HCO_State), pointer       :: HcoState
        real(sp),        pointer       :: NcArr(:,:,:,:)   ! (nlon,nlat,nlev,ntime)
        real(hp),        pointer       :: LonEdge(:)       ! Input lon edges
        real(hp),        pointer       :: LatEdge(:)       ! Input lat edges
        real(hp),        pointer       :: SigEdge(:,:,:)   ! Input sigma edges (may be NULL)
        logical,         intent(in)    :: IsModelLevel      ! Data on GEOS-Chem levels?
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(ListCont),  pointer       :: Lct
        integer,         intent(inout) :: RC
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
        character(len=*), parameter :: subname = 'HCO_ESMF_REGRID_DIRECT'

        integer :: nlon, nlat, nlev, ntime
        integer :: NX, NZ
        integer :: cache_idx
        integer :: L, T, I, esmf_rc

        ! ESMF field data pointers
        real(ESMF_KIND_R8), pointer :: srcPtr(:,:)   ! Source field data
        real(ESMF_KIND_R8), pointer :: dstPtr(:)     ! Destination field data

        ! Intermediate arrays
        real(r8), allocatable :: hRegridded(:,:)     ! (ncol, nlev) after horiz regrid
        real(r8), allocatable :: data_tgt(:,:)       ! (ncol, NZ) after vert regrid
        real(r8), allocatable :: sig_tgt(:,:)        ! (ncol, NZ+1) target sigma edges
        real(r8), allocatable :: sig_src_1d(:)       ! Source sigma edges (1D, uniform)

        ! Hardcoded GEOS-Chem 72-level sigma edges (same as in hcoio_read_pio_mod.F90)
        real(hp) :: GC_72_EDGE_SIGMA(73) = (/ &
          1.000000E+00, 9.849998E-01, 9.699136E-01, 9.548285E-01, 9.397434E-01, 9.246593E-01, &
          9.095741E-01, 8.944900E-01, 8.794069E-01, 8.643237E-01, 8.492406E-01, 8.341584E-01, &
          8.190762E-01, 7.989697E-01, 7.738347E-01, 7.487007E-01, 7.235727E-01, 6.984446E-01, &
          6.733175E-01, 6.356319E-01, 5.979571E-01, 5.602823E-01, 5.226252E-01, 4.849751E-01, &
          4.473417E-01, 4.097261E-01, 3.721392E-01, 3.345719E-01, 2.851488E-01, 2.420390E-01, &
          2.055208E-01, 1.746163E-01, 1.484264E-01, 1.261653E-01, 1.072420E-01, 9.115815E-02, &
          7.748532E-02, 6.573205E-02, 5.565063E-02, 4.702097E-02, 3.964964E-02, 3.336788E-02, &
          2.799704E-02, 2.341969E-02, 1.953319E-02, 1.624180E-02, 1.346459E-02, 1.112953E-02, &
          9.171478E-03, 7.520355E-03, 6.135702E-03, 4.981002E-03, 4.023686E-03, 3.233161E-03, &
          2.585739E-03, 2.057735E-03, 1.629410E-03, 1.283987E-03, 1.005675E-03, 7.846040E-04, &
          6.089317E-04, 4.697755E-04, 3.602270E-04, 2.753516E-04, 2.082408E-04, 1.569208E-04, &
          1.184308E-04, 8.783617E-05, 6.513694E-05, 4.737232E-05, 3.256847E-05, 1.973847E-05, &
          9.869233E-06/)

        ! Get input array dimensions
        nlon  = size(NcArr, 1)
        nlat  = size(NcArr, 2)
        nlev  = size(NcArr, 3)
        ntime = size(NcArr, 4)
        NX    = HcoState%NX       ! = ncol (physics columns)
        NZ    = HcoState%NZ       ! = CAM vertical levels

        ! Get/create cached route handle for this input grid
        call HCO_RegridCache_GetRH(nlon, nlat, LonEdge, LatEdge, cache_idx, esmf_rc)
        ASSERT_(esmf_rc==ESMF_SUCCESS)

        ! Get pointers to the cached ESMF fields
        call ESMF_FieldGet(cache(cache_idx)%srcField2D, localDE=0, &
                           farrayPtr=srcPtr, rc=esmf_rc)
        ASSERT_(esmf_rc==ESMF_SUCCESS)

        call ESMF_FieldGet(cache(cache_idx)%dstField2D, localDE=0, &
                           farrayPtr=dstPtr, rc=esmf_rc)
        ASSERT_(esmf_rc==ESMF_SUCCESS)

        !-----------------------------------------------------------------------
        ! 2D data: horizontal regrid only
        !-----------------------------------------------------------------------
        if (Lct%Dct%Dta%SpaceDim == 2) then

            ! Ensure output array is allocated
            call FileData_ArrCheck(HcoState%Config, Lct%Dct%Dta, &
                                   NX, 1, ntime, RC)
            if (RC /= HCO_SUCCESS) return

            do T = 1, ntime
                ! Fill source field
                do I = 1, nlat
                    srcPtr(:, I) = real(NcArr(:, I, 1, T), r8)
                enddo

                ! Zero destination
                dstPtr(:) = 0.0_r8

                ! Regrid
                call ESMF_FieldRegrid(cache(cache_idx)%srcField2D,  &
                                      cache(cache_idx)%dstField2D,  &
                                      cache(cache_idx)%rh2D,        &
                                      termorderflag=ESMF_TERMORDER_SRCSEQ, &
                                      rc=esmf_rc)
                ASSERT_(esmf_rc==ESMF_SUCCESS)

                ! Store in HEMCO data container (ncol, 1) layout
                do I = 1, NX
                    Lct%Dct%Dta%V2(T)%Val(I, 1) = real(dstPtr(I), sp)
                enddo
            enddo

        !-----------------------------------------------------------------------
        ! 3D data: horizontal regrid per-level, then vertical regrid per-column
        !-----------------------------------------------------------------------
        else if (Lct%Dct%Dta%SpaceDim == 3) then

            ! Ensure output array is allocated (ncol, 1, NZ)
            call FileData_ArrCheck(HcoState%Config, Lct%Dct%Dta, &
                                   NX, 1, NZ, ntime, RC)
            if (RC /= HCO_SUCCESS) return

            ! Allocate intermediate arrays
            allocate(hRegridded(NX, nlev))

            ! For 3D with vertical regridding, allocate target arrays once
            if (nlev > 1) then
                allocate(data_tgt(NX, NZ))
                allocate(sig_tgt(NX, NZ+1))

                ! Compute target sigma edges from HcoState pressure edges
                ! sigma = PEDGE / PSFC where PSFC = PEDGE(:,:,1)
                do I = 1, NX
                    do L = 1, NZ+1
                        if (associated(HcoState%Grid%PEDGE%Val)) then
                            sig_tgt(I, L) = HcoState%Grid%PEDGE%Val(I, 1, L) &
                                          / HcoState%Grid%PEDGE%Val(I, 1, 1)
                        else
                            ! Fallback: uniform sigma spacing (should not happen)
                            sig_tgt(I, L) = 1.0_r8 - real(L-1, r8) / real(NZ, r8)
                        endif
                    enddo
                enddo
            endif

            do T = 1, ntime

                ! Step 1: Horizontal ESMF regrid for each input level
                do L = 1, nlev
                    ! Fill source field with this level's data
                    do I = 1, nlat
                        srcPtr(:, I) = real(NcArr(:, I, L, T), r8)
                    enddo

                    dstPtr(:) = 0.0_r8

                    call ESMF_FieldRegrid(cache(cache_idx)%srcField2D,  &
                                          cache(cache_idx)%dstField2D,  &
                                          cache(cache_idx)%rh2D,        &
                                          termorderflag=ESMF_TERMORDER_SRCSEQ, &
                                          rc=esmf_rc)
                    ASSERT_(esmf_rc==ESMF_SUCCESS)

                    hRegridded(:, L) = dstPtr(1:NX)
                enddo

                ! Step 2: Vertical regrid per-column from input levels to CAM levels
                if (nlev > 1) then

                    ! Determine source sigma edges and perform vertical regrid
                    if (IsModelLevel) then
                        ! GEOS-Chem level data: use hardcoded sigma edges
                        allocate(sig_src_1d(nlev + 1))
                        sig_src_1d(1:nlev+1) = real(GC_72_EDGE_SIGMA(1:nlev+1), r8)

                        call HCO_VertRegrid_3D(NX, nlev, NZ,       &
                                               hRegridded, data_tgt, &
                                               sig_tgt,              &
                                               sig_src_1d=sig_src_1d)
                        deallocate(sig_src_1d)

                    else if (associated(SigEdge)) then
                        ! Real-coordinate data: sigma from file
                        ! SigEdge is (nlon, nlat, nlev+1) on the input grid.
                        ! After horizontal regridding, use a representative profile.
                        ! Average the sigma edges across the input horizontal domain
                        ! (sigma is typically uniform across the domain for most datasets).
                        allocate(sig_src_1d(nlev + 1))
                        do L = 1, nlev + 1
                            sig_src_1d(L) = 0.0_r8
                            do I = 1, min(size(SigEdge, 1), nlon)
                                sig_src_1d(L) = sig_src_1d(L) + &
                                    real(SigEdge(I, 1, L), r8)
                            enddo
                            sig_src_1d(L) = sig_src_1d(L) / real(min(size(SigEdge, 1), nlon), r8)
                        enddo

                        call HCO_VertRegrid_3D(NX, nlev, NZ,       &
                                               hRegridded, data_tgt, &
                                               sig_tgt,              &
                                               sig_src_1d=sig_src_1d)
                        deallocate(sig_src_1d)
                    else
                        ! No sigma info — assume input levels map to model levels
                        ! (direct copy for as many levels as available)
                        data_tgt = 0.0_r8
                        do L = 1, min(nlev, NZ)
                            data_tgt(:, L) = hRegridded(:, L)
                        enddo
                    endif

                    ! Store in HEMCO data container (ncol, 1, NZ)
                    do L = 1, NZ
                        do I = 1, NX
                            Lct%Dct%Dta%V3(T)%Val(I, 1, L) = real(data_tgt(I, L), sp)
                        enddo
                    enddo

                else
                    ! Single-level 3D data: just store the horizontally-regridded data
                    do I = 1, NX
                        Lct%Dct%Dta%V3(T)%Val(I, 1, 1) = real(hRegridded(I, 1), sp)
                    enddo
                endif

            enddo ! T

            ! Cleanup
            if (nlev > 1) then
                deallocate(data_tgt)
                deallocate(sig_tgt)
            endif
            deallocate(hRegridded)

        endif ! SpaceDim

        RC = HCO_SUCCESS

    end subroutine HCO_ESMF_REGRID_DIRECT
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_RegridCache_Cleanup
!
! !DESCRIPTION: Destroys all cached ESMF objects to free memory.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_RegridCache_Cleanup( RC )
!
! !USES:
!
        use ESMF, only: ESMF_FieldDestroy, ESMF_GridDestroy, ESMF_RouteHandleDestroy
!
! !OUTPUT PARAMETERS:
!
        integer, intent(out) :: RC
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
        character(len=*), parameter :: subname = 'HCO_RegridCache_Cleanup'
        integer :: n, esmf_rc

        RC = ESMF_SUCCESS

        do n = 1, nCached
            if (cache(n)%initialized) then
                call ESMF_RouteHandleDestroy(cache(n)%rh2D, rc=esmf_rc)
                call ESMF_FieldDestroy(cache(n)%srcField2D, rc=esmf_rc)
                call ESMF_FieldDestroy(cache(n)%dstField2D, rc=esmf_rc)
                call ESMF_GridDestroy(cache(n)%srcGrid, rc=esmf_rc)
                cache(n)%initialized = .false.
            endif
        enddo
        nCached = 0

    end subroutine HCO_RegridCache_Cleanup
!EOC
end module hco_esmf_regrid_cache
