
        MODULE TIDEFACMODULE
            USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR
            IMPLICIT NONE

            TYPE TIDEFAC
                TYPE(C_PTR),PRIVATE :: ptr
                CONTAINS
                    PROCEDURE            :: amplitude
                    PROCEDURE            :: frequency
                    PROCEDURE            :: earthTideReductionFactor
                    PROCEDURE            :: nodeFactor
                    PROCEDURE            :: equilibriumArgument
                    PROCEDURE            :: nodefactorCorrection
                    PROCEDURE            :: astronomicArgument
                    PROCEDURE            :: addConstituent
                    PROCEDURE            :: referenceTime
                    PROCEDURE            :: setReferenceTime
                    PROCEDURE            :: calculateWithDt
                    PROCEDURE            :: calculateWithDate
                    PROCEDURE            :: calculateWithTwoDates
                    GENERIC,PUBLIC       :: calculate => calculateWithDt,calculateWithDate,&
                                                         calculateWithTwoDates
            END TYPE TIDEFAC

            INTERFACE TIDEFAC
                PROCEDURE :: constructor
            END INTERFACE TIDEFAC

            INTERFACE
                TYPE(C_PTR) FUNCTION c_createTidefac() BIND(C,NAME="createTidefac") RESULT(ptr)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR
                    IMPLICIT NONE
                END FUNCTION c_createTidefac

                SUBROUTINE c_purgeTideFac() BIND(C,NAME="purgeTidefac")
                    IMPLICIT NONE
                END SUBROUTINE c_purgeTideFac

                SUBROUTINE c_setReferenceTime(ptr,year,month,day,hour,minute,second) &
                        BIND(C,NAME="setReferenceTime")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: year,month,day,hour,minute,second
                END SUBROUTINE c_setReferenceTime

                SUBROUTINE c_referenceTime(ptr,year,month,day,hour,minute,second) &
                        BIND(C,NAME="referenceTime")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE :: ptr
                    INTEGER(C_INT),INTENT(OUT)   :: year,month,day,hour,minute,second
                END SUBROUTINE c_referenceTime

                INTEGER(C_INT) FUNCTION c_addConstituent(ptr,harmonicName) BIND(C,NAME="addConstituent") &
                        RESULT(ierr)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_CHAR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE      :: ptr
                    CHARACTER(KIND=C_CHAR),INTENT(IN) :: harmonicName
                END FUNCTION c_addConstituent

                REAL(8) FUNCTION c_amplitude(ptr,idx) BIND(C,NAME="amplitude") RESULT(amp)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_amplitude
                
                REAL(8) FUNCTION c_frequency(ptr,idx) BIND(C,NAME="frequency") RESULT(freq)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_frequency
                
                REAL(8) FUNCTION c_earthTideReductionFactor(ptr,idx) &
                        BIND(C,NAME="earthTideReductionFactor") RESULT(etrf)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_earthTideReductionFactor
                
                REAL(8) FUNCTION c_nodeFactor(ptr,idx) BIND(C,NAME="nodeFactor") RESULT(nf)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_nodeFactor
                
                REAL(8) FUNCTION c_equilibriumArgument(ptr,idx) BIND(C,NAME="equilibriumArgument") RESULT(ea)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_equilibriumArgument
                
                REAL(8) FUNCTION c_nodefactorCorrection(ptr,idx) BIND(C,NAME="nodefactorCorrection") RESULT(ea)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_nodefactorCorrection
                
                REAL(8) FUNCTION c_astronomicArgument(ptr,idx) BIND(C,NAME="astronomicArgument") RESULT(aa)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: idx
                END FUNCTION c_astronomicArgument

                SUBROUTINE c_calculateWithDt(ptr,dt,lat) BIND(C,NAME="calculateWithDt")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_DOUBLE
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE :: ptr
                    REAL(C_DOUBLE),VALUE         :: dt
                    REAL(C_DOUBLE),VALUE         :: lat
                END SUBROUTINE c_calculateWithDt

                SUBROUTINE c_calculateWithDate(ptr,year,month,day,hour,minute,second,lat) BIND(C,NAME="calculateWithDate")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT,C_DOUBLE
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: year,month,day,hour,minute,second
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: lat
                END SUBROUTINE c_calculateWithDate
                
                SUBROUTINE c_calculateWithTwoDates(ptr,year1,month1,day1,hour1,minute1,second1,&
                                                       year2,month2,day2,hour2,minute2,second2,&
                                                       lat) BIND(C,NAME="calculateWithTwoDates")
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR,C_INT,C_DOUBLE
                    IMPLICIT NONE
                    TYPE(C_PTR),INTENT(IN),VALUE    :: ptr
                    INTEGER(C_INT),INTENT(IN),VALUE :: year1,month1,day1,hour1,minute1,second1
                    INTEGER(C_INT),INTENT(IN),VALUE :: year2,month2,day2,hour2,minute2,second2
                    REAL(C_DOUBLE),INTENT(IN),VALUE :: lat
                END SUBROUTINE c_calculateWithTwoDates

            END INTERFACE

            CONTAINS
            
                FUNCTION constructor() RESULT(this)
                    USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_PTR
                    IMPLICIT NONE
                    TYPE(TIDEFAC) :: this
                    this%ptr = c_createTideFac()
                END FUNCTION constructor
                    
                REAL(8) FUNCTION amplitude(this,idx) RESULT(amp)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    amp = c_amplitude(this%ptr,idx)
                END FUNCTION amplitude
                
                REAL(8) FUNCTION frequency(this,idx) RESULT(freq)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    freq = c_frequency(this%ptr,idx)
                END FUNCTION frequency
                
                REAL(8) FUNCTION earthTideReductionFactor(this,idx) RESULT(etrf)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    etrf = c_earthTideReductionFactor(this%ptr,idx)
                END FUNCTION earthTideReductionFactor

                REAL(8) FUNCTION nodeFactor(this,idx) RESULT(nf)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    nf = c_nodeFactor(this%ptr,idx)
                END FUNCTION nodeFactor
                
                REAL(8) FUNCTION equilibriumArgument(this,idx) RESULT(ea)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    ea = c_equilibriumArgument(this%ptr,idx)
                END FUNCTION equilibriumArgument
                
                REAL(8) FUNCTION nodefactorCorrection(this,idx) RESULT(aa)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    aa = c_nodefactorCorrection(this%ptr,idx)
                END FUNCTION nodefactorCorrection
                
                REAL(8) FUNCTION astronomicArgument(this,idx) RESULT(aa)
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    INTEGER,INTENT(IN)           :: idx
                    aa = c_astronomicArgument(this%ptr,idx)
                END FUNCTION astronomicArgument

                INTEGER FUNCTION addConstituent(this,harmonicName) RESULT(ierr)
                    USE,INTRINSIC                :: ISO_C_BINDING,ONLY:C_PTR,C_NULL_CHAR
                    IMPLICIT NONE
                    CLASS(TIDEFAC),INTENT(INOUT) :: this
                    CHARACTER(*)                 :: harmonicName
                    ierr = c_addConstituent(this%ptr,harmonicName//C_NULL_CHAR)
                END FUNCTION addConstituent

                INTEGER FUNCTION referenceTime(this,year,month,day,hour,minute,second) &
                        RESULT(IERR)
                    USE,INTRINSIC                :: ISO_C_BINDING,ONLY:C_PTR,C_INT
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

                SUBROUTINE purgeTidefac() 
                    IMPLICIT NONE
                    CALL c_purgeTidefac()
                END SUBROUTINE purgeTidefac

        END MODULE TIDEFACMODULE
