
        MODULE TIDEFACMODULE
            USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG
            IMPLICIT NONE

            TYPE TIDEFAC
                INTEGER(C_LONG),PRIVATE :: ptr
                CONTAINS
                    PROCEDURE            :: delete
                    PROCEDURE            :: generateLatitudeGrid
                    PROCEDURE            :: getInterpolationFactors
                    PROCEDURE            :: amplitudeSingle
                    PROCEDURE            :: frequencySingle
                    PROCEDURE            :: earthTideReductionFactorSingle
                    PROCEDURE            :: nodeFactorSingle
                    PROCEDURE            :: equilibriumArgumentSingle
                    PROCEDURE            :: nodefactorCorrectionSingle
                    PROCEDURE            :: astronomicArgumentSingle
                    PROCEDURE            :: amplitudeGrid
                    PROCEDURE            :: frequencyGrid
                    PROCEDURE            :: earthTideReductionFactorGrid
                    PROCEDURE            :: nodeFactorGrid
                    PROCEDURE            :: equilibriumArgumentGrid
                    PROCEDURE            :: nodefactorCorrectionGrid
                    PROCEDURE            :: astronomicArgumentGrid
                    PROCEDURE            :: addConstituent
                    PROCEDURE            :: referenceTime
                    PROCEDURE            :: setReferenceTime
                    PROCEDURE            :: calculateWithDt
                    PROCEDURE            :: calculateWithDate
                    PROCEDURE            :: calculateWithTwoDates
                    PROCEDURE            :: calculateGrid
                    GENERIC,PUBLIC       :: calculate => calculateWithDt,calculateWithDate,&
                                                         calculateWithTwoDates,calculateGrid
                    GENERIC,PUBLIC       :: amplitude => amplitudeSingle,amplitudeGrid
                    GENERIC,PUBLIC       :: frequency => frequencySingle,frequencyGrid
                    GENERIC,PUBLIC       :: earthTideReductionFactor => earthTideReductionFactorSingle,&
                                                earthTideReductionFactorGrid
                    GENERIC,PUBLIC       :: nodeFactor => nodeFactorSingle,nodeFactorGrid
                    GENERIC,PUBLIC       :: equilibriumArgument => equilibriumArgumentSingle,&
                                                equilibriumArgumentGrid
                    GENERIC,PUBLIC       :: nodeFactorCorrection => nodeFactorCorrectionSingle,&
                                                nodeFactorCorrectionGrid
                    GENERIC,PUBLIC       :: astronomicArgument => astronomicArgumentSingle,&
                                                astronomicArgumentGrid
            END TYPE TIDEFAC

            INTERFACE TIDEFAC
                PROCEDURE :: constructor
            END INTERFACE TIDEFAC

            INTERFACE
                INTEGER(C_LONG) FUNCTION c_createTidefac() BIND(C,NAME="createTidefac") RESULT(ptr)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG
                    IMPLICIT NONE
                END FUNCTION c_createTidefac

                SUBROUTINE c_deleteTidefac(ptr) BIND(C,NAME="deleteTidefac")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG
                    INTEGER(C_LONG),INTENT(IN) :: ptr
                END SUBROUTINE c_deleteTidefac

                SUBROUTINE c_purgeTideFac() BIND(C,NAME="purgeTidefac")
                    IMPLICIT NONE
                END SUBROUTINE c_purgeTideFac

                INTEGER(C_INT) FUNCTION c_generateLatitudeGrid(ptr,latmin,latmax,resolution) &
                        BIND(C,NAME="generateLatitudeGrid") RESULT(ierr)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT,C_DOUBLE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: latmin
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: latmax
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: resolution
                END FUNCTION c_generateLatitudeGrid

                INTEGER(C_INT) FUNCTION c_getInterpolationFactors(ptr,latitude,gridIndex,weight) &
                        BIND(C,NAME="getInterpolationFactors") RESULT(ierr)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_LONG,C_DOUBLE,C_INT
                    INTEGER(C_LONG),INTENT(IN),VALUE      :: ptr
                    REAL(C_DOUBLE),INTENT(IN),VALUE   :: latitude
                    INTEGER(C_INT),INTENT(OUT)        :: gridIndex
                    REAL(C_DOUBLE),INTENT(OUT)        :: weight
                END FUNCTION c_getInterpolationFactors
                    
                SUBROUTINE c_setReferenceTime(ptr,year,month,day,hour,minute,second) &
                        BIND(C,NAME="setReferenceTime")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: year,month,day,hour,minute,second
                END SUBROUTINE c_setReferenceTime

                SUBROUTINE c_referenceTime(ptr,year,month,day,hour,minute,second) &
                        BIND(C,NAME="referenceTime")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE :: ptr
                    INTEGER(C_INT),INTENT(OUT)   :: year,month,day,hour,minute,second
                END SUBROUTINE c_referenceTime

                INTEGER(C_INT) FUNCTION c_addConstituent(ptr,harmonicName) BIND(C,NAME="addConstituent") &
                        RESULT(ierr)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_CHAR,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE      :: ptr
                    CHARACTER(KIND=C_CHAR),INTENT(IN) :: harmonicName
                END FUNCTION c_addConstituent

                REAL(8) FUNCTION c_amplitudeSingle(ptr,idx) BIND(C,NAME="amplitude") RESULT(amp)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_amplitudeSingle
                
                REAL(8) FUNCTION c_amplitudeGrid(ptr,grid,idx) BIND(C,NAME="amplitudeGrid") RESULT(amp)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_amplitudeGrid
                
                REAL(8) FUNCTION c_frequencySingle(ptr,idx) BIND(C,NAME="frequency") RESULT(freq)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_frequencySingle
                
                REAL(8) FUNCTION c_frequencyGrid(ptr,grid,idx) BIND(C,NAME="frequencyGrid") RESULT(freq)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                END FUNCTION c_frequencyGrid
                
                REAL(8) FUNCTION c_earthTideReductionFactorSingle(ptr,idx) &
                        BIND(C,NAME="earthTideReductionFactor") RESULT(etrf)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_earthTideReductionFactorSingle
                
                REAL(8) FUNCTION c_earthTideReductionFactorGrid(ptr,grid,idx) &
                        BIND(C,NAME="earthTideReductionFactorGrid") RESULT(etrf)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_earthTideReductionFactorGrid
                
                REAL(8) FUNCTION c_nodeFactorSingle(ptr,idx) BIND(C,NAME="nodeFactor") RESULT(nf)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_nodeFactorSingle
                
                REAL(8) FUNCTION c_nodeFactorGrid(ptr,grid,idx) BIND(C,NAME="nodeFactorGrid") RESULT(nf)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_nodeFactorGrid
                
                REAL(8) FUNCTION c_equilibriumArgumentSingle(ptr,idx) BIND(C,NAME="equilibriumArgument") RESULT(ea)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_equilibriumArgumentSingle
                
                REAL(8) FUNCTION c_equilibriumArgumentGrid(ptr,grid,idx) BIND(C,NAME="equilibriumArgumentGrid") RESULT(ea)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_equilibriumArgumentGrid
                
                REAL(8) FUNCTION c_nodefactorCorrectionSingle(ptr,idx) BIND(C,NAME="nodefactorCorrection") RESULT(ea)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_nodefactorCorrectionSingle

                REAL(8) FUNCTION c_nodefactorCorrectionGrid(ptr,grid,idx) BIND(C,NAME="nodefactorCorrectionGrid") RESULT(ea)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_nodefactorCorrectionGrid
                
                REAL(8) FUNCTION c_astronomicArgumentSingle(ptr,idx) BIND(C,NAME="astronomicArgument") RESULT(aa)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_astronomicArgumentSingle
                
                REAL(8) FUNCTION c_astronomicArgumentGrid(ptr,grid,idx) BIND(C,NAME="astronomicArgumentGrid") RESULT(aa)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: grid
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_astronomicArgumentGrid

                SUBROUTINE c_calculateWithDt(ptr,dt,lat) BIND(C,NAME="calculateWithDt")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_DOUBLE
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE :: ptr
                    REAL(C_DOUBLE),VALUE         :: dt
                    REAL(C_DOUBLE),VALUE         :: lat
                END SUBROUTINE c_calculateWithDt

                SUBROUTINE c_calculateWithDate(ptr,year,month,day,hour,minute,second,lat) BIND(C,NAME="calculateWithDate")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT,C_DOUBLE
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: year,month,day,hour,minute,second
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: lat
                END SUBROUTINE c_calculateWithDate
                
                SUBROUTINE c_calculateWithTwoDates(ptr,year1,month1,day1,hour1,minute1,second1,&
                                                       year2,month2,day2,hour2,minute2,second2,&
                                                       lat) BIND(C,NAME="calculateWithTwoDates")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_INT,C_DOUBLE
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: year1,month1,day1,hour1,minute1,second1
                    INTEGER(C_INT),INTENT(IN),VALUE :: year2,month2,day2,hour2,minute2,second2
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: lat
                END SUBROUTINE c_calculateWithTwoDates

                SUBROUTINE c_calculateGrid(ptr,dt) BIND(C,NAME="calculateGrid")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG,C_DOUBLE
                    IMPLICIT NONE
                    INTEGER(C_LONG),INTENT(IN),VALUE    :: ptr
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: dt
                END SUBROUTINE c_calculateGrid

            END INTERFACE

            CONTAINS
            
                FUNCTION constructor() RESULT(this)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG
                    IMPLICIT NONE
                    TYPE(TIDEFAC) :: this
                    this%ptr = c_createTideFac()
                END FUNCTION constructor

                SUBROUTINE delete(this)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_LONG
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    call c_deleteTidefac(this%ptr)
                END SUBROUTINE delete
                    
                REAL(8) FUNCTION amplitudeSingle(this,idx) RESULT(amp)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    amp = c_amplitudeSingle(this%ptr,idx)
                END FUNCTION amplitudeSingle

                REAL(8) FUNCTION amplitudeGrid(this,grid,idx) RESULT(amp)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    amp = c_amplitudeGrid(this%ptr,grid,idx)
                END FUNCTION amplitudeGrid
                
                REAL(8) FUNCTION frequencySingle(this,idx) RESULT(freq)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    freq = c_frequencySingle(this%ptr,idx)
                END FUNCTION frequencySingle

                REAL(8) FUNCTION frequencyGrid(this,grid,idx) RESULT(freq)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    freq = c_frequencyGrid(this%ptr,grid,idx)
                END FUNCTION frequencyGrid

                REAL(8) FUNCTION earthTideReductionFactorSingle(this,idx) RESULT(etrf)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    etrf = c_earthTideReductionFactorSingle(this%ptr,idx)
                END FUNCTION earthTideReductionFactorSingle
                
                REAL(8) FUNCTION earthTideReductionFactorGrid(this,grid,idx) RESULT(etrf)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    etrf = c_earthTideReductionFactorGrid(this%ptr,grid,idx)
                END FUNCTION earthTideReductionFactorGrid

                REAL(8) FUNCTION nodeFactorSingle(this,idx) RESULT(nf)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    nf = c_nodeFactorSingle(this%ptr,idx)
                END FUNCTION nodeFactorSingle
                
                REAL(8) FUNCTION nodeFactorGrid(this,grid,idx) RESULT(nf)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    nf = c_nodeFactorGrid(this%ptr,grid,idx)
                END FUNCTION nodeFactorGrid
                
                REAL(8) FUNCTION equilibriumArgumentSingle(this,idx) RESULT(ea)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    ea = c_equilibriumArgumentSingle(this%ptr,idx)
                END FUNCTION equilibriumArgumentSingle
                
                REAL(8) FUNCTION equilibriumArgumentGrid(this,grid,idx) RESULT(ea)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    ea = c_equilibriumArgumentGrid(this%ptr,grid,idx)
                END FUNCTION equilibriumArgumentGrid
                
                REAL(8) FUNCTION nodefactorCorrectionSingle(this,idx) RESULT(aa)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    aa = c_nodefactorCorrectionSingle(this%ptr,idx)
                END FUNCTION nodefactorCorrectionSingle
                
                REAL(8) FUNCTION nodefactorCorrectionGrid(this,grid,idx) RESULT(aa)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    aa = c_nodefactorCorrectionGrid(this%ptr,grid,idx)
                END FUNCTION nodefactorCorrectionGrid
                
                REAL(8) FUNCTION astronomicArgumentSingle(this,idx) RESULT(aa)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    aa = c_astronomicArgumentSingle(this%ptr,idx)
                END FUNCTION astronomicArgumentSingle
                
                REAL(8) FUNCTION astronomicArgumentGrid(this,grid,idx) RESULT(aa)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: grid
                    INTEGER,INTENT(IN)           :: idx
                    aa = c_astronomicArgumentGrid(this%ptr,grid,idx)
                END FUNCTION astronomicArgumentGrid

                INTEGER FUNCTION addConstituent(this,harmonicName) RESULT(ierr)
                    USE,INTRINSIC                :: ISO_C_BINDING,ONLY:C_LONG,C_NULL_CHAR
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    CHARACTER(*)                 :: harmonicName
                    ierr = c_addConstituent(this%ptr,harmonicName//C_NULL_CHAR)
                END FUNCTION addConstituent

                INTEGER FUNCTION referenceTime(this,year,month,day,hour,minute,second) &
                        RESULT(IERR)
                    USE,INTRINSIC                :: ISO_C_BINDING,ONLY:C_LONG,C_INT
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER(C_INT),INTENT(OUT)   :: year,month,day,hour,minute,second
                    CALL c_referenceTime(this%ptr,year,month,day,hour,minute,second)
                    ierr = 0
                END FUNCTION referenceTime

                INTEGER FUNCTION setReferenceTime(this,year,month,day,hour,minute,second) &
                        RESULT(IERR)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)    :: year,month,day,hour,minute,second
                    CALL c_setReferenceTime(this%ptr,year,month,day,hour,minute,second)
                    ierr = 0
                END FUNCTION setReferenceTime

                INTEGER FUNCTION calculateWithDt(this,dt,latitude) RESULT(ierr)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    REAL(8),INTENT(IN) :: dt
                    REAL(8),INTENT(IN) :: latitude
                    CALL c_calculateWithDt(this%ptr,dt,latitude)
                END FUNCTION calculateWithDt
                
                INTEGER FUNCTION calculateWithDate(this,year,month,day,hour,minute,second,latitude) RESULT(ierr)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN) :: year,month,day,hour,minute,second
                    REAL(8),INTENT(IN) :: latitude
                    CALL c_calculateWithDate(this%ptr,year,month,day,hour,minute,second,latitude)
                END FUNCTION calculateWithDate
                
                INTEGER FUNCTION calculateWithTwoDates(this,year1,month1,day1,hour1,minute1,second1,&
                                                            year2,month2,day2,hour2,minute2,second2,&
                                                            latitude) RESULT(ierr)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN) :: year1,month1,day1,hour1,minute1,second1
                    INTEGER,INTENT(IN) :: year2,month2,day2,hour2,minute2,second2
                    REAL(8),INTENT(IN) :: latitude
                    CALL c_calculateWithTwoDates(this%ptr,year1,month1,day1,hour1,minute1,second1,&
                                                          year2,month2,day2,hour2,minute2,second2,latitude)
                END FUNCTION calculateWithTwoDates

                INTEGER FUNCTION calculateGrid(this,dt) RESULT(ierr)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    REAL(8),INTENT(IN)           :: dt
                    CALL c_calculateGrid(this%ptr,dt)
                END FUNCTION calculateGrid

                INTEGER FUNCTION generateLatitudeGrid(this,latmin,latmax,resolution) RESULT(ierr)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    REAL(8),INTENT(IN)           :: latmin,latmax,resolution
                    ierr = c_generateLatitudeGrid(this%ptr,latmin,latmax,resolution)
                END FUNCTION generateLatitudeGrid

                INTEGER FUNCTION getInterpolationFactors(this,latitude,gridIndex,weight) RESULT(ierr)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    REAL(8),INTENT(IN)           :: latitude
                    INTEGER,INTENT(OUT)          :: gridIndex
                    REAL(8),INTENT(OUT)          :: weight
                    ierr = c_getInterpolationFactors(this%ptr,latitude,gridIndex,weight)
                END FUNCTION getInterpolationFactors

                SUBROUTINE purgeTidefac() 
                    IMPLICIT NONE
                    CALL c_purgeTidefac()
                END SUBROUTINE purgeTidefac

        END MODULE TIDEFACMODULE
