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
!  The hook contract is owned by HEMCO (hco\_directregrid\_mod.F90):
!  HCO\_RegridCache\_Init registers HCO\_ESMF\_REGRID\_DIRECT with
!  HCO\_DirectRegrid\_Register, and HEMCO's hcoio\_read\_pio\_mod.F90
!  dispatches to it through HCO\_DirectRegrid\_Run whenever HcoDirectMode
!  is enabled.
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
    use ESMF,            only: ESMF_SUCCESS, ESMF_FAILURE
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
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version for direct-mode regridding
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    integer, parameter :: MAX_CACHED_GRIDS = 20

    ! Tolerance [deg] for matching input grid edges against a cached entry.
    real(r8), parameter :: EDGE_TOL = 1.0e-6_r8

    type :: RegridCacheEntry
        integer  :: nlon = 0            ! Input grid # longitudes
        integer  :: nlat = 0            ! Input grid # latitudes
        ! Full edge arrays are the uniqueness key. Dims + first edge alone
        ! would collide for same-size, same-origin grids with different
        ! spacing (e.g. non-uniform latitude grids, half-polar variants).
        real(r8), allocatable :: lonEdges(:)   ! (nlon+1) [deg]
        real(r8), allocatable :: latEdges(:)   ! (nlat+1) [deg]
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

    ! Number of PETs in the HEMCO communicator. Cached at init time to
    ! drive the source-grid decomposition search in HCO_RegridCache_GetRH.
    integer             :: cached_nPET = 1

contains
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_RegridCache_Init
!
! !DESCRIPTION: Initializes the regrid cache with a reference to the CAM
!  physics mesh, and registers this module's regridder with HEMCO's
!  direct-regrid hook (which enables HcoDirectMode on the HEMCO side).
!  Must be called after HCO_Grid_ESMF_CreateCAM.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_RegridCache_Init( mesh, ncol, nPET, mpicom, RC )
!
! !USES:
!
        use HCO_DirectRegrid_Mod, only: HCO_DirectRegrid_Register
!
! !INPUT PARAMETERS:
!
        type(ESMF_Mesh), intent(in) :: mesh
        integer,         intent(in) :: ncol   ! Local # physics columns
        integer,         intent(in) :: nPET   ! # PETs in HEMCO communicator
        integer,         intent(in) :: mpicom ! HEMCO MPI communicator
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

        phys_mesh   = mesh
        phys_ncol   = ncol
        cached_nPET = max(nPET, 1)
        nCached     = 0
        RC          = ESMF_SUCCESS

        ! Enable HcoDirectMode and hand HEMCO the regridding entry point plus
        ! the communicator (used by HEMCO for collective point-source lookup).
        call HCO_DirectRegrid_Register(HCO_ESMF_REGRID_DIRECT, mpiComm=mpicom)

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
                                      idx, RC, msg_out )
!
! !USES:
!
        use cam_logfile,  only: iulog
        use spmd_utils,   only: masterproc

        use ESMF,         only: ESMF_GridCreate1PeriDim, ESMF_INDEX_GLOBAL
        use ESMF,         only: ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
        use ESMF,         only: ESMF_GridAddCoord, ESMF_GridGetCoord, ESMF_GridGet
        use ESMF,         only: ESMF_TYPEKIND_R8, ESMF_KIND_R8
        use ESMF,         only: ESMF_MESHLOC_ELEMENT
        use ESMF,         only: ESMF_ArraySpec, ESMF_ArraySpecSet
        use ESMF,         only: ESMF_FieldCreate, ESMF_FieldRegridStore
        use ESMF,         only: ESMF_REGRIDMETHOD_CONSERVE, ESMF_REGRIDMETHOD_BILINEAR
        use ESMF,         only: ESMF_POLEMETHOD_NONE, ESMF_POLEMETHOD_ALLAVG
        use ESMF,         only: ESMF_RouteHandleDestroy, ESMF_FieldDestroy, ESMF_GridDestroy
        use ESMF,         only: ESMF_RouteHandleIsCreated
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
        character(len=*), optional, intent(out) :: msg_out
!
! !REVISION HISTORY:
!  09 Apr 2026 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
        character(len=*), parameter :: subname = 'HCO_RegridCache_GetRH'
        integer  :: n, i, j
        integer  :: srcLocalDECount
        integer  :: lbnd(2), ubnd(2)
        integer  :: decompNx, decompNy
        real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
        real(ESMF_KIND_R8), pointer :: coordX_E(:,:), coordY_E(:,:)
        type(ESMF_ArraySpec) :: arrayspec
        ! ESMF_FieldRegridStore dummies for srcTermProcessing / pipelineDepth
        ! are intent(inout) - cannot pass literal constants.
        integer  :: srcTermProc_arg
        integer  :: pipelineDepth_arg

        RC = ESMF_SUCCESS

        ! Check cache for an existing entry: dims plus the full edge arrays
        ! must match.
        do n = 1, nCached
            if (.not. cache(n)%initialized) cycle
            if (cache(n)%nlon /= nlon .or. cache(n)%nlat /= nlat) cycle
            if (maxval(abs(cache(n)%lonEdges - real(LonEdge, r8))) > EDGE_TOL) cycle
            if (maxval(abs(cache(n)%latEdges - real(LatEdge, r8))) > EDGE_TOL) cycle
            ! Cache hit
            idx = n
            return
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
        if (allocated(cache(idx)%lonEdges)) deallocate(cache(idx)%lonEdges)
        if (allocated(cache(idx)%latEdges)) deallocate(cache(idx)%latEdges)
        allocate(cache(idx)%lonEdges(nlon+1), cache(idx)%latEdges(nlat+1))
        cache(idx)%lonEdges(:) = real(LonEdge, r8)
        cache(idx)%latEdges(:) = real(LatEdge, r8)

        if (masterproc) then
            write(iulog,*) "HEMCO RegridCache: Creating route handle for input grid ", &
                           nlon, "x", nlat, " (entry ", idx, ")"
        endif

        !-----------------------------------------------------------------------
        ! Create ESMF Grid for the input file (rectilinear). Each PE reads the
        ! full input grid via PIO; ESMF handles the regrid-weight computation
        ! across PETs internally via FieldRegridStore.
        !
        ! Choose an explicit source-grid decomposition: decompNx*decompNy DEs
        ! with every DE at least 2 cells wide in both dims (ESMF conservative
        ! regridding rejects width<2 DEs). This mirrors the intermediate-grid
        ! decomposition search in hco_esmf_grid (and the WACCM-X ionosphere
        ! interface it derives from), rather than relying on ESMF's default
        ! regDecomp of (petCount, 1), which produces zero-width DEs whenever
        ! nlon < petCount.
        !
        ! Large grids benefit dramatically from parallel weight computation -
        ! a single-DE decomp forces all O(M*N*P) work onto one PET, serializing
        ! ESMF_FieldRegridStore setup. PIO reads the source array in full on
        ! every PE, so no data redistribution is needed before the src fill -
        ! see the srcLocalDECount gating in HCO_ESMF_REGRID_DIRECT.
        !-----------------------------------------------------------------------
        decompNx = 0
        decompNy = 0
        do i = 2, min(nlon, cached_nPET)
            if (mod(cached_nPET, i) /= 0) cycle
            j = cached_nPET / i
            if (j > nlat) cycle
            if (nlon/i > 1 .and. nlat/j > 1) then
                decompNx = i
                decompNy = j
                exit
            endif
        enddo
        ! Fall back to a 1-D latitude decomposition, then to a single DE
        ! (small grids and prime PET counts land here; non-source-owning PETs
        ! are handled downstream via the srcLocalDECount guards).
        if (decompNx == 0) then
            if (cached_nPET <= nlat/2) then
                decompNx = 1
                decompNy = cached_nPET
            else
                decompNx = 1
                decompNy = 1
            endif
        endif

        cache(idx)%srcGrid = ESMF_GridCreate1PeriDim(        &
            maxIndex=(/nlon, nlat/),                          &
            regDecomp=(/decompNx, decompNy/),                 &
            indexflag=ESMF_INDEX_GLOBAL,                      &
            rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if (masterproc) then
            write(iulog,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
                "HEMCO RegridCache: ", nlon, "x", nlat, &
                " - decomposed into ", decompNx, "x", decompNy, &
                " DEs across ", cached_nPET, " PETs"
        endif

        ! Add center and corner coordinates (collective - all PETs call)
        call ESMF_GridAddCoord(cache(idx)%srcGrid, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridAddCoord(cache(idx)%srcGrid, &
                               staggerloc=ESMF_STAGGERLOC_CORNER, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! With a single-DE decomp, only one PET holds a local DE for the
        ! source grid. Non-owning PETs have localDECount == 0 and cannot
        ! call ESMF_GridGetCoord(localDE=0). Skip the coord fill on those
        ! PETs - they still participate collectively in
        ! ESMF_FieldCreate / ESMF_FieldRegridStore below.
        call ESMF_GridGet(cache(idx)%srcGrid, localDECount=srcLocalDECount, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if (srcLocalDECount > 0) then
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
        endif

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
        ! Create route handle with the CONSERVE method (correct for
        ! area-weighted flux regridding). CONSERVE computes cell areas from
        ! the corner coordinates and rejects degenerate configurations, e.g.
        ! very coarse index-like source grids (gc_layers.nc, 4x4) whose
        ! near-pole corner cells collapse. For those grids ONLY, fall back
        ! to BILINEAR: they carry index/scale-factor style data where the
        ! horizontal field is effectively constant, so non-conservative
        ! interpolation is acceptable. For any grid large enough to carry
        ! real flux data, a CONSERVE failure is a hard error - silently
        ! degrading base emissions to a non-conservative method would defeat
        ! the purpose of the direct regridding path.
        !
        ! srcTermProcessing / pipelineDepth are intent(inout) - pass via
        ! local scalars rather than literal constants.
        !-----------------------------------------------------------------------
        srcTermProc_arg   = 0
        pipelineDepth_arg = 16
        call ESMF_FieldRegridStore(                                       &
            srcField=cache(idx)%srcField2D,                                &
            dstField=cache(idx)%dstField2D,                                &
            regridMethod=ESMF_REGRIDMETHOD_CONSERVE,                       &
            poleMethod=ESMF_POLEMETHOD_NONE,                               &
            routeHandle=cache(idx)%rh2D,                                   &
            srcTermProcessing=srcTermProc_arg,                             &
            pipelineDepth=pipelineDepth_arg, rc=RC)

        if (RC /= ESMF_SUCCESS) then
            if (min(nlon, nlat) > 4) then
                if (masterproc) then
                    write(iulog,'(a,i0,a,i0,a)') &
                        "HEMCO RegridCache: ERROR - conservative regrid weight "// &
                        "generation failed for ", nlon, "x", nlat, " source grid."
                endif
                if (present(msg_out)) then
                    msg_out = subname//': ESMF_FieldRegridStore (CONSERVE) failed'// &
                              ' for a non-degenerate source grid; refusing to fall'// &
                              ' back to non-conservative regridding for flux data.'
                endif
                return
            endif

            if (masterproc) then
                write(iulog,'(a,i0,a,i0,a)') &
                    "HEMCO RegridCache: WARNING - CONSERVE regrid failed for ", &
                    nlon, "x", nlat, &
                    " source grid; falling back to BILINEAR (non-conservative)."
            endif

            ! Destroy any partial route handle from the failed attempt
            if (ESMF_RouteHandleIsCreated(cache(idx)%rh2D)) then
                call ESMF_RouteHandleDestroy(cache(idx)%rh2D, rc=RC)
            endif

            RC = ESMF_SUCCESS
            call ESMF_FieldRegridStore(                                   &
                srcField=cache(idx)%srcField2D,                            &
                dstField=cache(idx)%dstField2D,                            &
                regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                   &
                poleMethod=ESMF_POLEMETHOD_ALLAVG,                         &
                routeHandle=cache(idx)%rh2D, rc=RC)
            ASSERT_(RC==ESMF_SUCCESS)
        endif

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
! !DESCRIPTION: Main entry point for direct ESMF regridding. Registered with
!  HEMCO's hco\_directregrid\_mod hook at init time and dispatched from
!  hcoio\_read\_pio\_mod.F90 in place of REGRID\_MAPA2A and HCO\_MESSY\_REGRID.
!
!  For 2D data: performs ESMF conservative horizontal regridding directly
!  from the input file grid to physics columns.
!
!  For 3D data: performs ESMF horizontal regridding per-level, then
!  sigma-to-sigma conservative vertical interpolation per-column. Source
!  sigma edges (SigEdge, prepared by the HEMCO reader) are themselves
!  horizontally regridded onto the columns when they vary across the input
!  domain, so terrain-following inputs keep their column-local vertical
!  placement.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_ESMF_REGRID_DIRECT( HcoState, NcArr, LonEdge, LatEdge, &
                                       SigEdge, Lct, RC, msg_out )
!
! !USES:
!
        use cam_logfile,         only: iulog
        use spmd_utils,          only: masterproc

        use ESMF,                only: ESMF_FieldRegrid, ESMF_FieldGet
        use ESMF,                only: ESMF_TERMORDER_SRCSEQ
        use ESMF,                only: ESMF_KIND_R8
        use ESMF,                only: ESMF_GridGet

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
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(ListCont),  pointer       :: Lct
        integer,         intent(inout) :: RC
!
! !OUTPUT PARAMETERS:
!
        character(len=*), optional, intent(out) :: msg_out
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
        integer :: L, T, I, J, esmf_rc
        integer :: srcLocalDECount
        integer :: srcLo(2), srcHi(2)     ! srcPtr local computational bounds
        logical :: sig_uniform

        ! Threshold for NetCDF _FillValue detection at srcPtr fill time.
        ! Realistic emission fluxes peak O(1) kg/m^2/s and scale factors are
        ! O(1-10). _FillValue entries in source NetCDFs are typically
        ! +/-9.97e+36 (netCDF default) or +/-1e30. Treating |val| > 1e15 as
        ! "missing" is generous vs. any real physical value and avoids ESMF
        ! CONSERVE mixing huge fill values into neighboring cells. Primary
        ! masking happens at read time in HCOIO's CheckMissVal; this is
        ! defense-in-depth for NaN (used as _FillValue by some files) which
        ! slips past equality-based masking.
        real(r8), parameter :: FILL_THRESHOLD = 1.0e15_r8

        ! Tolerance for detecting horizontally-uniform source sigma edges
        ! (sigma is dimensionless, O(1e-5..1)).
        real(hp), parameter :: SIG_UNIFORM_TOL = 1.0e-10_hp

        ! ESMF field data pointers
        real(ESMF_KIND_R8), pointer :: srcPtr(:,:)   ! Source field data
        real(ESMF_KIND_R8), pointer :: dstPtr(:)     ! Destination field data

        ! Intermediate arrays
        real(r8), allocatable :: hRegridded(:,:)     ! (ncol, nlev) after horiz regrid
        real(r8), allocatable :: data_tgt(:,:)       ! (ncol, NZ) after vert regrid
        real(r8), allocatable :: sig_tgt(:,:)        ! (ncol, NZ+1) target sigma edges
        real(r8), allocatable :: sig_src_1d(:)       ! (nlev+1) uniform source sigma
        real(r8), allocatable :: sig_src_col(:,:)    ! (ncol, nlev+1) per-column source sigma

        ! Get input array dimensions
        nlon  = size(NcArr, 1)
        nlat  = size(NcArr, 2)
        nlev  = size(NcArr, 3)
        ntime = size(NcArr, 4)
        NX    = HcoState%NX       ! = ncol (physics columns)
        NZ    = HcoState%NZ       ! = CAM vertical levels

        ! Get/create cached route handle for this input grid. Propagate
        ! failures (e.g. the refused CONSERVE->BILINEAR fallback) via
        ! RC/msg_out so HEMCO's dispatch can report them.
        call HCO_RegridCache_GetRH(nlon, nlat, LonEdge, LatEdge, cache_idx, &
                                   esmf_rc, msg_out)
        if (esmf_rc /= ESMF_SUCCESS) then
            RC = ESMF_FAILURE
            return
        endif

        ! Only PETs owning a source-side DE may touch srcPtr. Guard src-side
        ! FieldGet so non-owning PETs don't trip "localDeCount <= 0" errors.
        ! ESMF_FieldRegrid is collective and handles src->dst PET
        ! communication internally, so all PETs must still call it.
        call ESMF_GridGet(cache(cache_idx)%srcGrid, localDECount=srcLocalDECount, &
                          rc=esmf_rc)
        ASSERT_(esmf_rc==ESMF_SUCCESS)

        nullify(srcPtr)
        if (srcLocalDECount > 0) then
            call ESMF_FieldGet(cache(cache_idx)%srcField2D, localDE=0, &
                               farrayPtr=srcPtr, rc=esmf_rc)
            ASSERT_(esmf_rc==ESMF_SUCCESS)
            ! Under ESMF_INDEX_GLOBAL the returned pointer carries global
            ! bounds, so srcPtr(J,I) with J,I in [srcLo..srcHi] addresses
            ! the correct slice of NcArr (every PE has NcArr in full from
            ! PIO). For single-DE decomp srcLo/srcHi span the whole grid;
            ! for a parallel decomp they span this PET's tile.
            srcLo = lbound(srcPtr)
            srcHi = ubound(srcPtr)
        endif

        ! Destination field is on the physics mesh - every PET with physics
        ! columns has a local DE and needs dstPtr.
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
                ! Fill this PET's source tile (srcLo..srcHi). NcArr is the
                ! global array on every PE, so indexing with global J,I
                ! picks the correct tile. Nested NaN + fill-threshold
                ! checks - Fortran does NOT guarantee short-circuit
                ! evaluation of .or., so the NaN test must gate the
                ! abs() call to avoid FPE under strict compiler flags.
                if (srcLocalDECount > 0) then
                    do I = srcLo(2), srcHi(2)
                        do J = srcLo(1), srcHi(1)
                            if (NcArr(J, I, 1, T) /= NcArr(J, I, 1, T)) then
                                srcPtr(J, I) = 0.0_r8
                            else if (abs(real(NcArr(J, I, 1, T), r8)) > FILL_THRESHOLD) then
                                srcPtr(J, I) = 0.0_r8
                            else
                                srcPtr(J, I) = real(NcArr(J, I, 1, T), r8)
                            endif
                        enddo
                    enddo
                endif

                ! Zero destination
                dstPtr(:) = 0.0_r8

                ! Regrid (collective - all PETs must call)
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

            ! For 3D data with a vertical dimension, source sigma edges must
            ! have been prepared by the HEMCO reader (model-level table or file
            ! coordinates); the target sigma comes from the HEMCO vertical grid.
            if (nlev > 1) then
                if (.not. associated(SigEdge)) then
                    RC = ESMF_FAILURE
                    if (present(msg_out)) then
                        msg_out = subname//': multi-level 3D input requires source'// &
                                  ' sigma edges (SigEdge) but none were provided by'// &
                                  ' the HEMCO reader.'
                    endif
                    deallocate(hRegridded)
                    return
                endif
                if (.not. associated(HcoState%Grid%PEDGE%Val)) then
                    RC = ESMF_FAILURE
                    if (present(msg_out)) then
                        msg_out = subname//': HcoState%Grid%PEDGE is not set - the'// &
                                  ' HEMCO vertical grid must be established (via'// &
                                  ' HCO_CalcVertGrid) before 3D data is read.'
                    endif
                    deallocate(hRegridded)
                    return
                endif

                allocate(data_tgt(NX, NZ))
                allocate(sig_tgt(NX, NZ+1))

                ! Compute target sigma edges from HcoState pressure edges:
                ! sigma = PEDGE / PSFC where PSFC = PEDGE(:,:,1)
                do I = 1, NX
                    do L = 1, NZ+1
                        sig_tgt(I, L) = HcoState%Grid%PEDGE%Val(I, 1, L) &
                                      / HcoState%Grid%PEDGE%Val(I, 1, 1)
                    enddo
                enddo

                ! Source sigma edges on the physics columns. Most inputs carry
                ! horizontally-uniform sigma (all GEOS-Chem model-level data,
                ! and most pressure-level files); detect that case and use a
                ! single shared profile. Otherwise (terrain-following
                ! coordinates), horizontally regrid each sigma edge level onto
                ! the columns so the vertical interpolation stays column-local.
                sig_uniform = .true.
                do L = 1, nlev + 1
                    if (maxval(SigEdge(:,:,L)) - minval(SigEdge(:,:,L)) > SIG_UNIFORM_TOL) then
                        sig_uniform = .false.
                        exit
                    endif
                enddo

                if (sig_uniform) then
                    allocate(sig_src_1d(nlev + 1))
                    do L = 1, nlev + 1
                        sig_src_1d(L) = real(SigEdge(1, 1, L), r8)
                    enddo
                else
                    allocate(sig_src_col(NX, nlev + 1))
                    do L = 1, nlev + 1
                        if (srcLocalDECount > 0) then
                            do I = srcLo(2), srcHi(2)
                                do J = srcLo(1), srcHi(1)
                                    srcPtr(J, I) = real(SigEdge(J, I, L), r8)
                                enddo
                            enddo
                        endif

                        dstPtr(:) = 0.0_r8

                        call ESMF_FieldRegrid(cache(cache_idx)%srcField2D,  &
                                              cache(cache_idx)%dstField2D,  &
                                              cache(cache_idx)%rh2D,        &
                                              termorderflag=ESMF_TERMORDER_SRCSEQ, &
                                              rc=esmf_rc)
                        ASSERT_(esmf_rc==ESMF_SUCCESS)

                        sig_src_col(:, L) = dstPtr(1:NX)
                    enddo
                endif
            endif

            do T = 1, ntime

                ! Step 1: Horizontal ESMF regrid for each input level
                do L = 1, nlev
                    ! Fill this PET's source tile (srcLo..srcHi). See the 2D
                    ! fill loop above for the NaN/fill-clamp rationale and
                    ! the global-vs-local bounds contract.
                    if (srcLocalDECount > 0) then
                        do I = srcLo(2), srcHi(2)
                            do J = srcLo(1), srcHi(1)
                                if (NcArr(J, I, L, T) /= NcArr(J, I, L, T)) then
                                    srcPtr(J, I) = 0.0_r8
                                else if (abs(real(NcArr(J, I, L, T), r8)) > FILL_THRESHOLD) then
                                    srcPtr(J, I) = 0.0_r8
                                else
                                    srcPtr(J, I) = real(NcArr(J, I, L, T), r8)
                                endif
                            enddo
                        enddo
                    endif

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

                    if (sig_uniform) then
                        call HCO_VertRegrid_3D(NX, nlev, NZ,       &
                                               hRegridded, data_tgt, &
                                               sig_tgt,              &
                                               sig_src_1d=sig_src_1d)
                    else
                        call HCO_VertRegrid_3D(NX, nlev, NZ,       &
                                               hRegridded, data_tgt, &
                                               sig_tgt,              &
                                               sig_src_3d=sig_src_col)
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
            if (allocated(sig_src_1d))  deallocate(sig_src_1d)
            if (allocated(sig_src_col)) deallocate(sig_src_col)
            if (allocated(data_tgt))    deallocate(data_tgt)
            if (allocated(sig_tgt))     deallocate(sig_tgt)
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
! !DESCRIPTION: Destroys all cached ESMF objects to free memory, and
!  deregisters the direct-regrid hook from HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_RegridCache_Cleanup( RC )
!
! !USES:
!
        use ESMF, only: ESMF_FieldDestroy, ESMF_GridDestroy, ESMF_RouteHandleDestroy
        use HCO_DirectRegrid_Mod, only: HCO_DirectRegrid_Reset
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
            if (allocated(cache(n)%lonEdges)) deallocate(cache(n)%lonEdges)
            if (allocated(cache(n)%latEdges)) deallocate(cache(n)%latEdges)
        enddo
        nCached = 0

        call HCO_DirectRegrid_Reset()

    end subroutine HCO_RegridCache_Cleanup
!EOC
end module hco_esmf_regrid_cache
