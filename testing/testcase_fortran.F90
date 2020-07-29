

        PROGRAM TIDEFAC_FORTRAN
            USE TIDEFACMODULE
            IMPLICIT NONE
            TYPE(TIDEFAC) :: tide
            REAL(8) :: correctAmp(8) = (/ 0.24289000000000002D0, &
                                       0.11342D0, 0.046544999999999996D0, &
                                       0.0307625D0, 0.14204250000000002D0, &
                                       0.1008475D0, 0.0193135D0, 0.04708D0 /)
            REAL(8) :: correctEq(8) = (/ 241.82802987977877D0, 359.88858705269575D0, &
                                         166.2219210330479D0, 97.65123191473414D0, &
                                         318.8245122078793D0, 278.40076786446616D0, &
                                         202.89697390456098D0, 50.63607437499934D0 /)
            REAL(8) :: correctEtrf(8) = (/ 0.693D0,0.693D0,0.693D0,0.693D0,0.736D0,&
                                           0.695D0,0.695D0,0.706D0 /)
            REAL(8) :: correctFreq(8) = (/ 0.00014049900478554353D0,&
                                           0.00014538592669112763D0,&
                                           0.000137881010907552030D0,&
                                           0.000145909525466725930D0,&
                                           7.295476273336297D0*10e-05,6.754424205218055D0*10e-05,&
                                           6.492624817418905D0*10e-05,7.26056968829641D0*10e-05 /)
            REAL(8) :: correctNodeFac(8) = (/ 1.0088910309243924D0, 0.9994881350641318D0, &
                                              1.010033667119077D0, 0.9493529816818805D0, &
                                              0.9865471373144313D0, 0.9761350226537572D0, &
                                              0.9748772611466947D0, 1.0016946305592558D0 /)

            REAL(8) :: tol = 100D0*EPSILON(1D0)
            INTEGER :: I
            INTEGER :: ierr
            INTEGER :: year,month,day,hour,minute,second
            INTEGER :: yearq,monthq,dayq,hourq,minuteq,secondq

            year = 2011
            month = 11
            day = 1
            hour = 0
            minute = 0
            second = 0

            tide = TideFac()
            ierr = tide%setReferenceTime(2011,11,1,0,0,0)

            ierr = tide%referenceTime(yearq,monthq,dayq,hourq,minuteq,secondq)
            IF( year.ne.yearq.or.month.ne.monthq.or.day.ne.dayq.or.hour.ne.hourq.or. &
                minute.ne.minuteq.or.second.ne.secondq)THEN
                WRITE(*,'(A)') "ERROR: Date was not set/read correctly"
                CALL EXIT(1)
            ENDIF

            ierr = tide%addConstituent("M2")
            ierr = tide%addConstituent("S2")
            ierr = tide%addConstituent("N2")
            ierr = tide%addConstituent("K2")
            ierr = tide%addConstituent("K1")
            ierr = tide%addConstituent("O1")
            ierr = tide%addConstituent("Q1")
            ierr = tide%addConstituent("P1")
            ierr = tide%calculate(0D0,29.7046D0)

            DO I = 1,8
                IF(tide%amplitude(i)-correctAmp(i).GT.tol)THEN
                    WRITE(*,'(A)') "ERROR: Amplitude incorrect"
                    WRITE(*,*) I,tide%amplitude(i),correctAmp(i)
                    CALL EXIT(1)
                ENDIF
                IF(tide%frequency(i)-correctFreq(i).GT.tol)THEN
                    WRITE(*,'(A)') "ERROR: Frequency incorrect"
                    WRITE(*,*) I,tide%amplitude(i),correctAmp(i)
                    CALL EXIT(1)
                ENDIF
                IF(tide%earthTideReductionFactor(i)-correctEtrf(i).GT.tol)THEN
                    WRITE(*,'(A)') "ERROR: ETRF incorrect"
                    WRITE(*,*) I,tide%earthTideReductionFactor(i),correctEtrf(i)
                    CALL EXIT(1)
                ENDIF
                IF(tide%equilibriumArgument(i)-correctEq(i).GT.tol)THEN
                    WRITE(*,'(A)') "ERROR: Equilibrium Argument incorrect"
                    WRITE(*,*) I,tide%equilibriumArgument(i),correctEq(i)
                    CALL EXIT(1)
                ENDIF
                IF(tide%nodeFactor(i)-correctNodeFac(i).GT.tol)THEN
                    WRITE(*,'(A)') "ERROR: Node factor incorrect"
                    WRITE(*,*) I,tide%nodeFactor(i),correctNodeFac(i)
                    CALL EXIT(1)
                ENDIF
            ENDDO

            CALL purgeTidefac()

        END PROGRAM TIDEFAC_FORTRAN


