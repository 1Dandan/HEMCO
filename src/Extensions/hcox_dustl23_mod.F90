!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_dustl23_mod.F90
!
! !DESCRIPTION: Module hcox\_dustl23\_mod.F90 contains routines and
!  variables from Danny M. Leung's dust emission scheme.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_DustL23_mod
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCOX_TOOLS_MOD
  USE HCOX_State_MOD, ONLY : Ext_State
  USE HCO_State_MOD,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_DustL23_Run
  PUBLIC :: HCOX_DustL23_Init
  PUBLIC :: HCOX_DustL23_Final
!
! !REVISION HISTORY:
!  2 May 2024 - Dandan Zhang - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  ! MyInst is the extension-specific derived type. It should hold all module
  ! variables and arrays that are required to compute the emissions.
  ! For instance, if the extension relies on an input field read through the
  ! HEMCO configuration file (e.g. MY_INPUT_FIELD), the data array pointer
  ! to that field should be listed within the instance and NOT outside of it.
  ! This ensures that the same extension can be invoked in various instances,
  ! all of them potentially pointing to different data fields.
  TYPE :: MyInst
   ! Fields required by module
   INTEGER                         :: Instance
   INTEGER                         :: ExtNr            ! Extension num for DustL23
   INTEGER                         :: ExtNrAlk         ! Extension num for DustAlk
   INTEGER, ALLOCATABLE            :: HcoIDs(:)        ! tracer IDs for DustDead
   INTEGER, ALLOCATABLE            :: HcoIDsAlk(:)     ! tracer IDs for DustAlk
   INTEGER                         :: nSpc             ! # of species
   REAL(sp), ALLOCATABLE           :: SpcScl(:)        ! Species scale factors
   CHARACTER(LEN=31), ALLOCATABLE  :: SpcNames(:)
   INTEGER                         :: nSpcAlk          ! # of species
   CHARACTER(LEN=31), ALLOCATABLE  :: SpcNamesAlk(:)
   CHARACTER(LEN=61), ALLOCATABLE  :: SpcScalFldNme(:) ! Names of scale factor fields

   ! Other fields
   REAL(hp),  ALLOCATABLE            :: DMT_MIN(:)       ! Bin size min diameter [m]
   REAL(hp),  ALLOCATABLE            :: DMT_MAX(:)       ! Bin size max diameter [m]
   
   !---------------------------------------
   ! 2-D pointers pointing to netCDF arrays
   !---------------------------------------

   ! Land cover map
   REAL(hp), POINTER          :: A_bare    (:,:) => NULL() ! The fraction of barren and sparsely vegetated land cover [unitless]
   REAL(hp), POINTER          :: A_veg     (:,:) => NULL() ! The fraction of short vegetation land cover [unitless]
   REAL(hp), POINTER          :: A_water   (:,:) => NULL() ! The fraction of water bodies land cover [unitless]
   REAL(hp), POINTER          :: A_wetland (:,:) => NULL() ! The fraction of permenent wetland land cover [unitless]
   
   ! Soil textture map
   REAL(hp), POINTER          :: f_clay    (:,:) => NULL() ! The fraction of clay content in topmost soil [unitless]
   REAL(hp), POINTER          :: f_sand    (:,:) => NULL() ! The fraction of sand content in topmost soil [unitless]
   
   ! Surface roughness length due to rocks [m]
   REAL(hp), POINTER          :: roughness_r    (:,:) => NULL()
   
   ! Source function at 0.25 x 0.3125 degree
   REAL(hp), POINTER          :: SRCE_FUNC(:,:) => NULL()

   TYPE(MyInst), POINTER           :: NextInst => NULL()

  END TYPE MyInst

  ! Pointer to all instances
  TYPE(MyInst), POINTER            :: AllInst => NULL()


    !---------------------------------------
    ! MODULE PARAMETER
    !---------------------------------------
    INTEGER, PARAMETER   :: NBINS = 4        ! # of dust bins
    INTEGER, PARAMETER   :: TNSPEC = 5       ! # of dust bins + 1
    ! Fundamental physical constants
    REAL(hp),  PARAMETER   :: CST_VON_KRM      = 0.386_hp
    REAL(hp),  PARAMETER   :: SPC_HEAT_DRY_AIR = 1005.0_hp
    
    ! Other dust parameters
    REAL(hp),  PARAMETER   :: D_p        = 127.0e-6_hp ! Median diameter of soil particle [m]
    REAL(hp),  PARAMETER   :: rho_p      = 2650.0_hp ! Soil particle density [kg m-3]
    REAL(hp),  PARAMETER   :: rho_w      = 1000.0_hp ! Water density [kg m-3]
    REAL(hp),  PARAMETER   :: rho_a0     = 1.225_hp  ! Standard air density at sea level and 15 degrees C [kg m-3]
    REAL(hp),  PARAMETER   :: LAI_thr    = 0.3_hp    ! Threshold LAI [unitless]

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_dustl23_Run
!
! !DESCRIPTION: Subroutine HcoX\_dustl23\_Run is the driver routine
! for the HEMCO DustL23 extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_DustL23_Run( ExtState, HcoState, RC )
!
! !USES:
!   
    USE HCO_CALC_MOD,      ONLY : HCO_EvalFld, HCO_CalcEmis
    USE HCO_FLUXARR_MOD,   ONLY : HCO_EmisAdd
    USE HCO_CLOCK_MOD,     ONLY : HcoClock_Get
    USE HCO_CLOCK_MOD,     ONLY : HcoClock_First
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState    ! Module options
    TYPE(HCO_State), POINTER       :: HcoState    ! Hemco state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  2 May 2024 - Dandan Zhang - Revised from template for DustL23
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Local variables
    LOGICAL                :: ERR
    INTEGER                :: N
    ! Total dust emission flux [kg/m2/s]
    REAL(hp), TARGET       :: TFLUX(HcoState%NX, HcoState%NY)
    ! Flux array [kg/m2/s]
    REAL(hp), TARGET       :: FLUX(HcoState%NX, HcoState%NY, TNSPEC)
    ! Flux array for dust alkalinity [kg/m2/s]
    REAL(hp), TARGET       :: FLUX_ALK(HcoState%NX,HcoState%NY,TNSPEC)

    TYPE(MyInst), POINTER :: Inst => NULL()
    CHARACTER(LEN=255)    :: MSG, LOC

    !=================================================================
    ! HCOX_DustL23_RUN begins here!
    !=================================================================
    LOC = 'HCOX_DustL23_RUN (HCOX_DustL23_MOD.F90)'
    ! Return if extension disabled
    IF ( ExtState%DustL23 <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Get pointer to this instance. Varible Inst contains all module
    ! variables for the current instance. The instance number is
    ! ExtState%DustL23
    ! Get instance
    Inst => NULL()
    CALL InstGet ( ExtState%DustL23, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find DustL23 instance Nr. ', ExtState%DustL23
       CALL HCO_ERROR(MSG,RC)
       RETURN
    ENDIF

    !=================================================================
    ! Module code comes below
    !=================================================================
    CALL HCO_EvalFld( HcoState, 'DustL23_SRCE_FUNC', Inst%SRCE_FUNC, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_LandCover_bare', Inst%A_bare, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_LandCover_veg', Inst%A_veg, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_LandCover_water', Inst%A_water, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_LandCover_wetland', Inst%A_wetland, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_Soil_clay', Inst%f_clay, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_Soil_sand', Inst%f_sand, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'DustL23_roughness_r', Inst%roughness_r, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !=================================================================
    ! DustL23 Emission Scheme
    !=================================================================
    ! Error check
    ERR = .FALSE.

    CALL CAL_DUSTL23_EmisFlux(HcoState, ExtState, Inst, TFLUX, RC)
    ! Error check
    IF ( RC /= HCO_SUCCESS ) THEN
      ERR = .TRUE.
      RETURN
    ENDIF
    FLUX(:,:,1) = TFLUX
    FLUX(:,:,2) = TFLUX * 0.0766_hp
    FLUX(:,:,3) = TFLUX * 0.1924_hp
    FLUX(:,:,4) = TFLUX * 0.3491_hp
    FLUX(:,:,5) = TFLUX * 0.3819_hp

    ! Include DUST Alkalinity SOURCE, assuming an alkalinity
    ! of 4% by weight [kg].                  !tdf 05/10/08
    !tdf with 3% Ca, there's also 1% equ. Mg, makes 4%
    IF ( Inst%ExtNrAlk > 0 ) THEN
      FLUX_ALK = 0.04_hp * FLUX
    ENDIF

    ! Error check
    IF ( ERR ) THEN
      RC = HCO_FAIL
      RETURN
    ENDIF

   !=================================================================
   ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
   !=================================================================
    DO N = 1, TNSPEC
      IF ( Inst%HcoIDs(N) > 0 ) THEN
        ! Add to emissions array
        CALL HCO_EmisAdd( HcoState, FLUX(:,:,N), Inst%HcoIDs(N), RC, ExtNr=Inst%ExtNr )
        IF ( RC /= HCO_SUCCESS ) THEN
          WRITE(MSG,*) 'HCO_EmisAdd error: dust bin ', N
          CALL HCO_ERROR(MSG, RC )
          RETURN
        ENDIF
      ENDIF

      IF ( Inst%ExtNrAlk > 0 ) THEN
        IF ( Inst%HcoIDsAlk(N) > 0 ) THEN
          ! Add to dust alkalinity emissions array
          CALL HCO_EmisAdd( HcoState, FLUX_Alk(:,:,N), Inst%HcoIDsAlk(N), RC, ExtNr=Inst%ExtNrAlk )
           IF ( RC /= HCO_SUCCESS ) THEN
              WRITE(MSG,*) 'HCO_EmisAdd error: dust alk bin ', N
              CALL HCO_ERROR(MSG, RC )
              RETURN
           ENDIF

        ENDIF
     ENDIF
    ENDDO

    ! Cleanup
    Inst => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_DustL23_Run
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustL23_Init
!
! !DESCRIPTION: Subroutine HcoX\_DustL23\_Init initializes the HEMCO
! DUST\_L23 extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_DustL23_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_STATE_MOD,      ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod,    ONLY : GetExtSpcVal
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state
    INTEGER,          INTENT(INOUT) :: RC

! !REVISION HISTORY:
!  06 May 2024 - Dandan Zhang - Initial version for DustL23
!  25 Nov 2013 - C. Keller   - Now a HEMCO extension
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER          :: Inst => NULL()
    INTEGER                        :: ExtNr, N, AS
    CHARACTER(LEN=255)             :: MSG, LOC

    !=================================================================
    ! HCOX_DustL23_INIT begins here!
    !=================================================================
    LOC = 'HCOX_DustL23_INIT (HCOX_DUSTL23_MOD.F90)'

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Create instance for this simulation. Link instance number to the ExtState object
    ! for future reference to the instance. See InstCreate for more details.
    CALL InstCreate ( ExtNr, ExtState%DustL23, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( 'Cannot create DustL23 instance', RC )
       RETURN
    ENDIF

    ! Check for dust alkalinity option
    Inst%ExtNrAlk = GetExtNr( HcoState%Config%ExtList, 'DustAlk')

    ! Get species IDs.
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, Inst%HcoIDs, Inst%SpcNames, Inst%nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Get the dust alkalinity species defined for DustAlk option
    IF ( Inst%ExtNrAlk > 0 ) THEN
      CALL HCO_GetExtHcoID( HcoState, Inst%ExtNrAlk, Inst%HcoIDsAlk, Inst%SpcNamesAlk, Inst%nSpcAlk, RC)
      IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 3', RC, THISLOC=LOC )
        RETURN
      ENDIF
    ENDIF

    ! Sanity check
    IF ( Inst%nSpc /= TNSPEC ) THEN
      MSG = 'DustL23 model does not have four species!'
      CALL HCO_ERROR(MSG, RC )
      RETURN
    ENDIF

    ! There must be at least one species
    IF ( Inst%nSpc == 0 ) THEN
       CALL HCO_ERROR ( 'No DustL23 species specified', RC )
       RETURN
    ENDIF

    ! Determine scale factor to be applied to each species. This is 1.00
    ! by default, but can be set in the HEMCO configuration file via setting
    ! Scaling_<SpcName>.
    CALL GetExtSpcVal( HcoState%Config, Inst%ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'Scaling', 1.0_sp, Inst%SpcScl, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Get species mask fields
    CALL GetExtSpcVal( HcoState%Config, Inst%ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'ScaleField', HCOX_NOSCALE, Inst%SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 5', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN

       ! Write the name of the extension regardless of the verbose setting
       msg = 'Using HEMCO extension: DustL23 (dust emission scheme)'
       IF ( HCO_IsVerb( HcoState%Config%Err ) ) THEN
          CALL HCO_Msg( HcoState%Config%Err, sep1='-' ) ! with separator
       ELSE
          CALL HCO_Msg( msg, verb=.TRUE.              ) ! w/o separator
       ENDIF

       ! Write all other messages as debug printout only
       MSG = ' - use the following species (Name, HcoID, Scaling):'
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       DO N = 1, Inst%nSpc
          WRITE(MSG,*) TRIM(Inst%SpcNames(N)), ', ', Inst%HcoIDs(N), ', ', Inst%SpcScl(N)
          CALL HCO_MSG( HcoState%Config%Err, MSG)
          WRITE(MSG,*) 'Apply scale field: ', TRIM(Inst%SpcScalFldNme(N))
          CALL HCO_MSG( HcoState%Config%Err, MSG)
       ENDDO
    ENDIF

    !-----------------------------------------------------------------
    ! Init module arrays
    !-----------------------------------------------------------------
    ALLOCATE( Inst%A_bare( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%A_bare!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%A_bare = 0.0_hp

    ALLOCATE( Inst%A_veg( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%A_veg!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%A_veg = 0.0_hp

    ALLOCATE( Inst%A_water( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%A_water!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%A_water = 0.0_hp

    ALLOCATE( Inst%A_wetland( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%A_wetland!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%A_wetland = 0.0_hp

    ALLOCATE( Inst%SRCE_FUNC( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%SRCE_FUNC!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%SRCE_FUNC = 0.0_hp

    ALLOCATE( Inst%roughness_r( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%roughness_r!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%roughness_r = 0.0_hp
    
    ALLOCATE( Inst%f_clay( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%f_clay!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%f_clay = 0.0_hp

    ALLOCATE( Inst%f_sand( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS /= 0 ) THEN
        msg = 'Could not allocate Inst%f_sand!'
        CALL HCO_ERROR( msg, RC, thisLoc=loc )
        RETURN
    ENDIF
    Inst%f_sand = 0.0_hp

    ! Bin size min diameter [m]
    ALLOCATE( Inst%DMT_MIN( NBINS ), STAT=AS )
    IF ( AS /= 0 ) THEN
      CALL HCO_ERROR ( 'DMT_MIN', RC )
      RETURN
    ENDIF
    Inst%DMT_MIN(1) = 0.2e-6_hp
    Inst%DMT_MIN(2) = 2.0e-6_hp
    Inst%DMT_MIN(3) = 3.6e-6_hp
    Inst%DMT_MIN(4) = 6.0e-6_hp

    ! Bin size max diameter [m]
    ALLOCATE( Inst%DMT_MAX( NBINS ), STAT=AS )
    IF ( AS /= 0 ) THEN
      CALL HCO_ERROR ( 'DMT_MAX', RC )
      RETURN
    ENDIF
    Inst%DMT_MAX(1) = 2.0e-6_hp
    Inst%DMT_MAX(2) = 3.6e-6_hp
    Inst%DMT_MAX(3) = 6.0e-6_hp
    Inst%DMT_MAX(4) = 1.2e-5_hp

    ! Activate met fields used by this extension
    ExtState%T2M%DoUse      = .TRUE.
    ExtState%PS%DoUse       = .TRUE.
    ExtState%GWETTOP%DoUse  = .TRUE.
    ExtState%FRSNO%DoUse    = .TRUE.
    ExtState%USTAR%DoUse    = .TRUE.
    ExtState%PBLH%DoUse     = .TRUE.
    ExtState%HFLUX%DoUse    = .TRUE.
    ExtState%LAI%DoUse    = .TRUE.

    ! Cleanup
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_DustL23_Init
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustL23_Final
!
! !DESCRIPTION: Subroutine HcoX\_DustL23\_Final finalizes the HEMCO
! DustL23 extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_DustL23_Final ( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  06 May 2024 - Dandan Zhang - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! HCOX_DustL23_FINAL begins here!
    !=================================================================
    CALL InstRemove ( ExtState%DustL23 )

  END SUBROUTINE HCOX_DustL23_Final
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate adds a new instance to the list of
!  instances, assigns a unique instance number to this new instance, and
!  archives this instance number to output argument Instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst  => NULL()
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------
    Inst%SRCE_FUNC       => NULL()
    Inst%A_bare          => NULL()
    Inst%A_veg           => NULL()
    Inst%A_water         => NULL()
    Inst%A_wetland       => NULL()
    Inst%f_clay          => NULL()
    Inst%f_sand          => NULL()
    Inst%roughness_r     => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove removes an instance from the list of
! instances.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Initialize
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       !---------------------------------------------------------------------
       ! Deallocate fields of Inst before popping off from the list
       ! in order to avoid memory leaks (Bob Yantosca (17 Aug 2022)
       !---------------------------------------------------------------------
       IF ( ALLOCATED( Inst%HcoIDs ) ) THEN
          DEALLOCATE ( Inst%HcoIDs )
       ENDIF

       IF ( ALLOCATED( Inst%SpcScl ) ) THEN
          DEALLOCATE ( Inst%SpcScl )
       ENDIF

       IF ( ALLOCATED( Inst%SpcNames ) ) THEN
          DEALLOCATE ( Inst%SpcNames )
       ENDIF

       IF ( ALLOCATED( Inst%SpcScalFldNme ) ) THEN
          DEALLOCATE( Inst%SpcScalFldNme  )
       ENDIF
       
       ! ----------------------------------------------------------------
       ! Type specific initialization statements follow below
       ! ----------------------------------------------------------------
       IF ( ASSOCIATED( Inst%SRCE_FUNC ) ) THEN
          DEALLOCATE(Inst%SRCE_FUNC )
       ENDIF
       Inst%SRCE_FUNC => NULL()
       
       IF ( ASSOCIATED( Inst%A_bare ) ) THEN
        DEALLOCATE(Inst%A_bare )
       ENDIF
       Inst%A_bare => NULL()

       IF ( ASSOCIATED( Inst%A_veg ) ) THEN
        DEALLOCATE(Inst%A_veg )
       ENDIF
       Inst%A_veg => NULL()

       IF ( ASSOCIATED( Inst%A_water ) ) THEN
        DEALLOCATE(Inst%A_water )
       ENDIF
       Inst%A_water => NULL()

       IF ( ASSOCIATED( Inst%A_wetland ) ) THEN
        DEALLOCATE(Inst%A_wetland )
       ENDIF
       Inst%A_wetland => NULL()

       IF ( ASSOCIATED( Inst%f_clay ) ) THEN
        DEALLOCATE(Inst%f_clay )
       ENDIF
       Inst%f_clay => NULL()

       IF ( ASSOCIATED( Inst%f_sand ) ) THEN
        DEALLOCATE(Inst%f_sand )
       ENDIF
       Inst%f_sand => NULL()

       IF ( ASSOCIATED( Inst%roughness_r ) ) THEN
        DEALLOCATE(Inst%roughness_r )
       ENDIF
       Inst%roughness_r => NULL()

       ! ----------------------------------------------------------------
       ! Pop off instance from list
       ! ----------------------------------------------------------------
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
    ENDIF

    ! Free pointers before exiting
    PrevInst => NULL()
    Inst     => NULL()

   END SUBROUTINE InstRemove
  
  SUBROUTINE CAL_THR_FRIC_VEL(HcoState, rho_a, f_clay, f_sand, theta, &
                             u_star_ft0, u_star_ft, u_star_it, u_star_st, RC)
    ! Description: calculate threshold friction velocities

    !----------------
    ! Arguments
    !----------------
    TYPE(HCO_State), POINTER :: HcoState
    REAL(hp),  INTENT(IN)  :: rho_a(HcoState%NX, HcoState%NY) ! surface air density [kg m-3]
    REAL(hp),  INTENT(IN)  :: f_clay(HcoState%NX, HcoState%NY) ! Soil clay fraction [unitless]
    REAL(hp),  INTENT(IN)  :: f_sand(HcoState%NX, HcoState%NY) ! Soil sand fraction [unitless] -> used to calcualte soil porosity
    REAL(hp),  INTENT(IN)  :: theta(HcoState%NX, HcoState%NY) ! Volumetric soil water content [unitless]
    REAL(hp),  INTENT(OUT) :: u_star_ft0(HcoState%NX, HcoState%NY) ! Dry fluid thershold friction velocity [m s-1]
    REAL(hp),  INTENT(OUT) :: u_star_ft(HcoState%NX, HcoState%NY) ! Wet fluid thershold friction velocity [m s-1]
    REAL(hp),  INTENT(OUT) :: u_star_it(HcoState%NX, HcoState%NY) ! Dynamic fluid thershold friction velocity [m s-1]
    REAL(hp),  INTENT(OUT) :: u_star_st(HcoState%NX, HcoState%NY) ! Standardized wet fluid thershold friction velocity [m s-1]
    INTEGER, INTENT(INOUT) :: RC

    !-----------------
    ! Local variables
    !-----------------
    REAL(hp),  PARAMETER   :: A         = 0.0123_hp
    REAL(hp),  PARAMETER   :: gamma     = 1.65e-4_hp ! [kg s-2]
    REAL(hp),  PARAMETER   :: B_it      = 0.82_hp

    ! Gravimetric soil moisture [unitless]
    REAL(hp)               :: w(HcoState%NX, HcoState%NY)

    ! Threshols gravimetric soil moisture [unitlss]
    REAL(hp)               :: w_t(HcoState%NX, HcoState%NY)

    ! Factor by which threhold velocity increases due to soil wetness
    REAL(hp)               :: f_m(HcoState%NX, HcoState%NY)
    
    ! Soil porosity [unitless]
    REAL(hp)               :: phi(HcoState%NX, HcoState%NY)
    ! calculate soil porosity: phi = 0.489 - 0.126 * f_sand
    phi = 0.489_hp - 0.126_hp * f_sand

    !=================================================================
    ! CAL_THR_FRIC_VEL begins here!
    !=================================================================

    ! Dry fluid threshold velocity [m s-1]: 
    ! calculate u_star_ft0 = sqrt(A * (rho_p * g * D_p + gamma / D_p) / rho_a)
    u_star_ft0(:,:) =  SQRT(A * (rho_p * HcoState%Phys%g0 * D_p + gamma / D_p) / rho_a(:,:)) ! [m s-1]

    ! Factor by which soil wetness enhancing threhold friction velocity
    ! calculate f_m = sqrt (1 + 1.21 * ((100 * (w - w_t)) ** 0.68)) for w > w_t; and f_m = 1 for w <= w_t
    !! calculate w = rho_w / rho_b * theta = rho_w / (rho_p * ( 1 - phi)) * theta
    w = rho_w / (rho_p * (1.0_hp - phi)) * theta
    
    !! calculate w_t = 0.01 * a * (17 * f_clay + 14 * f_clay ** 2) where a is a tuning factor and was set to be 1 / f_clay
    !! and thus w_t = 0.01 * (17 + 14 * f_clay)
    w_t = 0.01_hp * (17.0_hp + 14.0_hp * f_clay)

    ! calculate f_m [unitless]
    WHERE ( w <= w_t )
      f_m = 1.0_hp
    ELSEWHERE
      f_m = SQRT(1.0_hp + 1.21_hp * ((100.0_hp * (w - w_t) ** 0.68_hp)))
    ENDWHERE

    ! Wet threshold friction velocity [m s-1]
    u_star_ft = u_star_ft0 * f_m

    ! Dynamic threshold friction velocity [m s-1]
    ! calculate u_star_it = B_it * u_star_ft0
    u_star_it = B_it * u_star_ft0

    ! Standardized wet fluid thershold friction velocity [m s-1]
    u_star_st = u_star_ft * SQRT(rho_a / rho_a0)
    
    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE CAL_THR_FRIC_VEL
  
  SUBROUTINE CAL_DRAG_PART(HcoState, z_0a, LAI, A_r, A_v, &
                          f_eff_r, f_eff_v, F_eff, RC)
    ! Description: calculate drag partioning effects due to rocks and vegetation

    !----------------
    ! Arguments
    !----------------
    TYPE(HCO_State), POINTER :: HcoState
    REAL(hp),  INTENT(IN)  :: z_0a(HcoState%NX, HcoState%NY) ! surface roughness length due to rocks [m]
    REAL(hp),  INTENT(IN)  :: LAI(HcoState%NX, HcoState%NY) ! Leaf area index [unitless]
    REAL(hp),  INTENT(IN)  :: A_r(HcoState%NX, HcoState%NY) ! The fraction of barren and sparsely vegetated land cover [unitless]
    REAL(hp),  INTENT(IN)  :: A_v(HcoState%NX, HcoState%NY) ! The fraction of short vegetation land cover [unitless]
    REAL(hp),  INTENT(OUT) :: f_eff_r(HcoState%NX, HcoState%NY) ! The drag partitioning effects due to rocks [unitless]
    REAL(hp),  INTENT(OUT) :: f_eff_v(HcoState%NX, HcoState%NY) ! The drag partitioning effects due to short vegetation [unitless]
    REAL(hp),  INTENT(OUT) :: F_eff(HcoState%NX, HcoState%NY) ! The total drag partitioning effects due to rocks and short vegetation [unitless]
    INTEGER, INTENT(INOUT) :: RC

    !-----------------
    ! Local variables
    !-----------------
    ! parameters
    REAL(hp),  PARAMETER   :: b1        = 0.7_hp
    REAL(hp),  PARAMETER   :: b2        = 0.8_hp
    REAL(hp),  PARAMETER   :: X         = 10.0_hp ! [m]
    REAL(hp),  PARAMETER   :: f0        = 0.32_hp
    REAL(hp),  PARAMETER   :: c         = 4.8_hp

    ! Derived quantities
    ! smooth roughness length which quantifies the roughness of a bed of fine soil particles in the absence of roughness elements [m]
    REAL(hp),  PARAMETER   :: z_0s      = 2.0_hp * D_p / 30.0_hp

    ! variables
    REAL(hp)               :: K(HcoState%NX, HcoState%NY)
    ! calculate K = 2 * (1 / f_v - 1) = 2 * (LAI_thr / LAI - 1)
    K = 2.0_hp * (LAI_thr / LAI - 1.0_hp)
    WHERE (K < 0.0_hp)
      K = 0.0_hp
    ENDWHERE

    ! Calculate drag partioning effects due to rocks:
    ! f_eff_r = 1 - ln(z_0a / z_0s) / ln(b1 * (X / z_0s) ** b2)
    f_eff_r = 1.0d0 - LOG(z_0a / z_0s) / LOG(b1 * (X / z_0s) ** b2)
    WHERE (f_eff_r < 0.0d0)
      f_eff_r = 0.0d0
    ELSEWHERE (f_eff_r > 1.0d0)
      f_eff_r = 1.0d0
    ENDWHERE
    
    ! calculate drag partioning effects due to vegetation:
    ! f_eff_v = (K + f0 * c) / (K + c)
    f_eff_v = (K + f0 * c) / (K + c)

    ! calculate the weighted-mean drag partioning effects due to rocks and vegetation
    F_eff = (A_r * (f_eff_r ** 3.0_hp) + A_v * (f_eff_v ** 3.0_hp)) ** (1.0_hp/3.0_hp)
    
    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE CAL_DRAG_PART

  SUBROUTINE CAL_INTERM_FACTOR(HcoState, u_star_ft, u_star_it, u_star_s, &
                               PBLH, rho_a, TS, u_star, HFLUX, eta, RC)
    ! Description: calculate intermittency factor due to turbulence

    !----------------
    ! Arguments
    !----------------
    TYPE(HCO_State), POINTER :: HcoState
    REAL(hp),  INTENT(IN) :: u_star_ft(HcoState%NX, HcoState%NY) ! Wet fluid thershold friction velocity [m s-1]
    REAL(hp),  INTENT(IN) :: u_star_it(HcoState%NX, HcoState%NY) ! Dynamic fluid thershold friction velocity [m s-1]
    REAL(hp),  INTENT(IN) :: u_star_s(HcoState%NX, HcoState%NY)  ! Soil friction velocity [m s-1]
    REAL(hp),  INTENT(IN) :: PBLH(HcoState%NX, HcoState%NY) ! Planetary boundary layer height [m]
    REAL(hp),  INTENT(IN) :: rho_a(HcoState%NX, HcoState%NY) ! Surface air density [kg m-3]
    REAL(hp),  INTENT(IN) :: TS(HcoState%NX, HcoState%NY) ! Temperature at 10 meter [K]
    REAL(hp),  INTENT(IN) :: u_star(HcoState%NX, HcoState%NY) ! Friction velocity [m s-1]
    REAL(hp),  INTENT(IN) :: HFLUX(HcoState%NX, HcoState%NY) ! Sensible heat flux [W m-2]
    REAL(hp),  INTENT(OUT) :: eta(HcoState%NX, HcoState%NY) ! Intermittency factor with values in [0,1] [unitless]
    INTEGER, INTENT(INOUT) :: RC

    !-----------------
    ! Local variables
    !-----------------
    ! parameters
    REAL(hp),  PARAMETER   :: z_sal        = 0.1_hp  ! Saltation height [m]
    REAL(hp),  PARAMETER   :: z_0a_c       = 1.0e-4_hp ! take as constant for simplicity [m]

    ! variables
    ! derived friction velocities at saltation height [m s-1]
    REAL(hp)               :: u_ft(HcoState%NX, HcoState%NY)
    REAL(hp)               :: u_it(HcoState%NX, HcoState%NY)
    REAL(hp)               :: u_s(HcoState%NX, HcoState%NY)

    REAL(hp)               :: L(HcoState%NX, HcoState%NY) ! Monin-bukhov length [m]
    REAL(hp)               :: sigma(HcoState%NX, HcoState%NY) ! standard deviation of u_s due to turbulence [m s-1]
    REAL(hp)               :: alpha(HcoState%NX, HcoState%NY)
    REAL(hp)               :: P_ft(HcoState%NX, HcoState%NY)
    REAL(hp)               :: P_it(HcoState%NX, HcoState%NY)
    
    u_ft = u_star_ft / CST_VON_KRM * LOG(z_sal / z_0a_c)
    u_it = u_star_it / CST_VON_KRM * LOG(z_sal / z_0a_c)
    u_s = u_star_s / CST_VON_KRM * LOG(z_sal / z_0a_c)
    
    ! calculate sigma for instantaneous soil friction velocity:
    ! sigma = u_star_s * (12 - 0.5 * PBLH / L) ** (1/3) for (12 - 0.5 * PBLH / L)>=0
    !! calculate the Monin-bukhov length L = - rho_a * cp * T * u_star ** 3 / (k * g * H)
    L = - rho_a * SPC_HEAT_DRY_AIR * TS * u_star ** 3.0_hp / (CST_VON_KRM * HcoState%Phys%g0 * HFLUX)
    sigma = u_star_s * ((12.0_hp - 0.5_hp * PBLH / L) ** (1.0_hp/3.0_hp))
    
    ! calculate the fluid thresholf crossing fraction: 
    ! alpha = (exp ((u_ft ** 2 - u_it ** 2 - 2 * u_s * (u_ft - u_it)) / (2 * sigma ** 2)) + 1) ** (-1)
    alpha = (EXP (((u_ft ** 2.0_hp) - (u_it ** 2.0_hp) - 2.0_hp * u_s * (u_ft - u_it)) / (2.0_hp * (sigma ** 2.0_hp))) + 1.0_hp) ** (-1.0_hp)

    ! calculate the cumulative probability that instananeous soil friction velocity does not exceed 
    ! the wet fluid threshold u_ft or the impact threshold u_it of P_ft or P_it
    ! P_ft = 0.5 * (1 + erf((u_ft - u_s) / (sqrt(2) * sigma))); P_it = 0.5 * (1 + erf((u_it - u_s) / (sqrt(2) * sigma)))
    P_ft = 0.5_hp * (1.0_hp + ERF((u_ft - u_s) / (SQRT(2.0_hp) * sigma)))
    P_it = 0.5_hp * (1.0_hp + ERF((u_it - u_s) / (SQRT(2.0_hp) * sigma)))

    ! calculate intermittency factor: eta = 1 - P_ft + alpha * (P_ft - P_it)
    eta = 1.0_hp - P_ft + alpha * (P_ft - P_it)
    ! if eta is out of range of [0,1], then skip eta multipling by making the value as 1
    WHERE ((eta < 0.0_hp) .or. (eta > 1.0_hp) .or. (sigma < 0))
      eta = 1.0_hp
    ENDWHERE

    ! Return w/ success
    RC = HCO_SUCCESS
  END SUBROUTINE CAL_INTERM_FACTOR

  SUBROUTINE CAL_DUSTL23_EmisFlux(HcoState, ExtState, Inst, DUST_EMIS_FLUX, RC)
    ! Description: calculate DustL23 total emission flux [kg m-2 s-1]
    !----------------
    ! Arguments
    !----------------
    TYPE(HCO_State), POINTER      :: HcoState    ! Hemco state
    TYPE(Ext_State), POINTER      :: ExtState    ! Module options
    TYPE(MyInst),    POINTER      :: Inst        ! Specific instances of DustL23
    REAL(hp),  INTENT(OUT)        :: DUST_EMIS_FLUX(HcoState%NX, HcoState%NY) ! Total dust emission flux [kg m-2 s-1]
    INTEGER, INTENT(INOUT)        :: RC

    ! Local variables
    REAL(hp)                :: u_star_ft0(HcoState%NX, HcoState%NY) ! Dry fluid thershold friction velocity [m s-1]
    REAL(hp)                :: u_star_ft(HcoState%NX, HcoState%NY) ! Wet fluid thershold friction velocity [m s-1]
    REAL(hp)                :: u_star_it(HcoState%NX, HcoState%NY) ! Dynamic fluid thershold friction velocity [m s-1]
    REAL(hp)                :: u_star_st(HcoState%NX, HcoState%NY) ! Standardized wet fluid thershold friction velocity [m s-1]

    REAL(hp)                :: f_eff_r(HcoState%NX, HcoState%NY) ! The drag partitioning effects due to rocks [unitless]
    REAL(hp)                :: f_eff_v(HcoState%NX, HcoState%NY) ! The drag partitioning effects due to short vegetation [unitless]
    REAL(hp)                :: F_eff(HcoState%NX, HcoState%NY) ! The total drag partitioning effects due to rocks and short vegetation [unitless]
    
    REAL(hp)                :: eta(HcoState%NX, HcoState%NY) ! Intermittency factor with values in [0,1] [unitless]
    REAL(hp)                :: DUST_EMIS_FLUX_Tmp(HcoState%NX, HcoState%NY) ! Total dust emission flux [kg m-2 s-1]

    CHARACTER(LEN=255)    :: SUBLOC
    ! empirical constants
    REAL(hp),  PARAMETER    :: C_tune         = 0.05_hp    ! [unitless]
    REAL(hp),  PARAMETER    :: C_d0           = 4.4e-5_hp  ! [unitless]
    REAL(hp),  PARAMETER    :: C_e            = 2.0_hp    ! [unitless]
    REAL(hp),  PARAMETER    :: C_kappa        = 2.7_hp   ! [unitless]

    ! other constants
    REAL(hp),  PARAMETER    :: u_star_st0     = 0.16_hp  ! [m s-1]
  
    REAL(hp)        :: rho_a(HcoState%NX, HcoState%NY) ! surface air density [kg m-3]
    REAL(hp)        :: TS(HcoState%NX, HcoState%NY)    ! Surface temperature [K]
    REAL(hp)        :: PS(HcoState%NX, HcoState%NY)    ! Surface pressure [Pa]
    REAL(hp)        :: LAI(HcoState%NX, HcoState%NY) ! Leaf area index [unitless]

    REAL(hp)        :: C_d(HcoState%NX, HcoState%NY)            ! Soil erodibility coefficient [unitless]
    REAL(hp)        :: f_bare(HcoState%NX, HcoState%NY)         ! [unitless]
    REAL(hp)        :: u_star_s(HcoState%NX, HcoState%NY)       ! Soil surface friction velocity [m s-1]
    REAL(hp)        :: f_clay_capped(HcoState%NX, HcoState%NY)  ! [unitless]
    REAL(hp)        :: kappa(HcoState%NX, HcoState%NY)          ! Fragmentaion exponent [unitless]
    REAL(hp)        :: u_star_t(HcoState%NX, HcoState%NY)       ! Thershold friction velocity used [m s-1]
    
    TS = ExtState%T2M%Arr%Val
    PS = ExtState%PS%Arr%Val * 100.0_hp ! convert hPa to Pa
    rho_a = PS * (HcoState%Phys%AIRMW * 1.0e-3_hp) / (HcoState%Phys%RSTARG * TS)

    
    ! LAI = SUM(ExtState%XLAI%Arr%Val, DIM=3)
    LAI = ExtState%LAI%Arr%Val

    SUBLOC = 'CAL_THR_FRIC_VEL'
    CALL CAL_THR_FRIC_VEL(HcoState, rho_a, Inst%f_clay, Inst%f_sand, ExtState%GWETTOP%Arr%Val, &
                          u_star_ft0, u_star_ft, u_star_it, u_star_st, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
      CALL HCO_ERROR( 'ERROR', RC, THISLOC=SUBLOC )
      RETURN
    ENDIF

    SUBLOC = 'CAL_DRAG_PART'
    CALL CAL_DRAG_PART(HcoState, Inst%roughness_r, LAI, Inst%A_bare, Inst%A_veg, &
                       f_eff_r, f_eff_v, F_eff, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
      CALL HCO_ERROR( 'ERROR', RC, THISLOC=SUBLOC )
      RETURN
    ENDIF
    
    u_star_s = ExtState%USTAR%Arr%Val * F_eff

    SUBLOC = 'CAL_INTERM_FACTOR'
    CALL CAL_INTERM_FACTOR(HcoState, u_star_ft, u_star_it, u_star_s, &
                           ExtState%PBLH%Arr%Val, rho_a, TS, ExtState%USTAR%Arr%Val, ExtState%HFLUX%Arr%Val, eta, RC)
    IF ( RC /= HCO_SUCCESS ) THEN
      CALL HCO_ERROR( 'ERROR', RC, THISLOC=SUBLOC )
      RETURN
    ENDIF

    ! calculate C_d = C_d0 * exp (- C_e * (u_star_st - u_star_st0) / u_star_st0)
    C_d = C_d0 * EXP (- C_e * (u_star_st - u_star_st0) / u_star_st0)

    ! calculate f_bare = 1 - LAI / LAI_thr for LAI <= LAI_thr, and f_bare = 0 for LAI > LAI_thr
    f_bare = (1.0_hp - Inst%A_water - Inst%A_wetland) * (1.0_hp - ExtState%FRSNO%Arr%Val) * (1.0_hp - LAI / LAI_thr)
    WHERE (LAI > LAI_thr)
      f_bare = 0.0_hp
    ENDWHERE

    f_clay_capped = Inst%f_clay
    WHERE (Inst%f_clay > 0.2_hp)
      f_clay_capped = 0.2_hp
    ENDWHERE

    kappa = C_kappa * (u_star_st - u_star_st0) / u_star_st0
    WHERE (kappa>3.0_hp)
      kappa = 3.0_hp
    ENDWHERE
    
    u_star_t = u_star_it  
    DUST_EMIS_FLUX_Tmp = eta * Inst%SRCE_FUNC * C_tune * C_d * f_bare * f_clay_capped * \
        rho_a * ((u_star_s ** 2.0_hp) - (u_star_t ** 2.0_hp)) / u_star_st * \
        ((u_star_s / u_star_t) ** kappa)
    WHERE ((DUST_EMIS_FLUX_Tmp < 0.0_hp) .or. (u_star_s .LE. u_star_t))
      DUST_EMIS_FLUX_Tmp = 0.0_hp
    ENDWHERE

    WHERE (DUST_EMIS_FLUX_Tmp > 0.0_hp)
      DUST_EMIS_FLUX = DUST_EMIS_FLUX_Tmp 
    ELSEWHERE
      DUST_EMIS_FLUX = 0.0_hp
    ENDWHERE

    ! Return w/ success
    RC = HCO_SUCCESS
  END SUBROUTINE CAL_DUSTL23_EmisFlux
!EOC
END MODULE HCOX_DustL23_Mod
