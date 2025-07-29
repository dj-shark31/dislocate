MODULE FitModule

  USE disloc_elasticity_ani, ONLY : max_nd
  USE dDipoleModule, ONLY : max_nDDipole
  SAVE

  ! Variable=.true. if the corresponding parameters cDislo(:,:)
  !   or bDislo(:,:) have to be fitted
  LOGICAL, dimension(1:3,1:max_nd), private :: fit_cDislo
  LOGICAL, dimension(1:3,1:max_nd), private :: fit_bDislo
  LOGICAL, dimension(1:3,1:max_nDDipole), private :: fit_c1DDipole
  LOGICAL, dimension(1:3,1:max_nDDipole), private :: fit_c2DDipole
  LOGICAL, dimension(1:3,1:max_nDDipole), private :: fit_bDDipole

  ! Fractional tolerance in the function value such that failure to 
  ! decrease by more than this amount on one iteration signals doneness.
  REAL(kind(0.d0)) :: ftol

  ! Threshold for data points that are fitted
  REAL(kind(0.d0)) :: threshold

CONTAINS

  SUBROUTINE InitFit()

    IMPLICIT NONE

    fit_cDislo(:,:) = .FALSE.
    fit_bDislo(:,:) = .FALSE.
    fit_c1DDipole(:,:) = .FALSE.
    fit_c2DDipole(:,:) = .FALSE.
    fit_bDDipole(:,:) = .FALSE.
    fTol = 1.d-2
    threshold = 0.d0

  END SUBROUTINE InitFit

  SUBROUTINE ReadFit(inp)

    IMPLICIT NONE
    INTEGER, intent(in) :: inp  ! Input unit

    INTEGER, dimension(1:3,1:max_nd) :: cDislo, bDislo
    INTEGER, dimension(1:3,1:max_nDDipole) :: c1DDipole, c2DDipole, bDDipole

    NAMELIST /fit/ cDislo, bDislo, &
        c1DDipole, c2DDipole, bDDipole, &
        fTol, threshold

    ! Initialization
    cDislo(:,:)=0
    bDislo(:,:)=0
    c1DDipole(:,:)=0
    c2DDipole(:,:)=0
    bDDipole(:,:)=0

    ! Read input data
    READ(inp,nml=fit)

    ! cDislo(:,:)=0 => Variable not fitted
    fit_cDislo(:,:) = cDislo(:,:).NE.0
    fit_bDislo(:,:) = bDislo(:,:).NE.0
    fit_c1DDipole(:,:) = c1DDipole(:,:).NE.0
    fit_c2DDipole(:,:) = c2DDipole(:,:).NE.0
    fit_bDDipole(:,:) = bDDipole(:,:).NE.0

  END SUBROUTINE ReadFit

  SUBROUTINE PrintFit(out)

    USE disloc_elasticity_ani, ONLY : nd
    USE dDipoleModule, ONLY : nDDipole
    IMPLICIT NONE

    INTEGER, intent(in) :: out

    CHARACTER(len=15) :: phrase
    INTEGER :: nDim, n

    ! Number of free variables
    nDim = nDimFit()

    WRITE(out,*)
    WRITE(out,'(a,i0)') 'Fit: number of free variables: ', nDim
    WRITE(out,'(a)') '  Dislocation definition'
    IF (nd.GT.0) THEN
            DO n=1, nd
               IF ( ( Count(fit_cDislo(:,n)) + Count(fit_bDislo(:,n)) ).EQ.0) THEN
                       WRITE(out,'(a,i0,a)') '    dislo ', n, ': no free variable'
               ELSE
                       phrase=''
                       IF (fit_cDislo(1,n)) WRITE(phrase,'(2a)') Trim(phrase), ' x'
                       IF (fit_cDislo(2,n)) WRITE(phrase,'(2a)') Trim(phrase), ' y'
                       IF (fit_cDislo(3,n)) WRITE(phrase,'(2a)') Trim(phrase), ' z'
                       IF (fit_bDislo(1,n)) WRITE(phrase,'(2a)') Trim(phrase), ' bx'
                       IF (fit_bDislo(2,n)) WRITE(phrase,'(2a)') Trim(phrase), ' by'
                       IF (fit_bDislo(3,n)) WRITE(phrase,'(2a)') Trim(phrase), ' bz'
                       WRITE(out,'(a,i0,3a)') '    dislo ', n, ':', phrase, ' free'
               END IF
            END DO
    END IF
    IF (nDDipole.GT.0) THEN
            DO n=1, nDDipole
               IF ( ( Count(fit_c1DDipole(:,n)) + Count(fit_c2DDipole(:,n)) &
                        + Count(fit_bDDipole(:,n)) ).EQ.0) THEN
                       WRITE(out,'(a,i0,a)') '    dipole ', n, ': no free variable'
               ELSE
                       phrase=''
                       IF (fit_c1DDipole(1,n)) WRITE(phrase,'(2a)') Trim(phrase), ' x1'
                       IF (fit_c1DDipole(2,n)) WRITE(phrase,'(2a)') Trim(phrase), ' y1'
                       IF (fit_c1DDipole(3,n)) WRITE(phrase,'(2a)') Trim(phrase), ' z1'
                       IF (fit_c2DDipole(1,n)) WRITE(phrase,'(2a)') Trim(phrase), ' x2'
                       IF (fit_c2DDipole(2,n)) WRITE(phrase,'(2a)') Trim(phrase), ' y2'
                       IF (fit_c2DDipole(3,n)) WRITE(phrase,'(2a)') Trim(phrase), ' z2'
                       IF (fit_bDDipole(1,n))  WRITE(phrase,'(2a)') Trim(phrase), ' bx'
                       IF (fit_bDDipole(2,n))  WRITE(phrase,'(2a)') Trim(phrase), ' by'
                       IF (fit_bDDipole(3,n))  WRITE(phrase,'(2a)') Trim(phrase), ' bz'
                       WRITE(out,'(a,i0,3a)') '    dipole ', n, ':', phrase, ' free'
               END IF
            END DO
    END IF
    WRITE(out,'(a,g20.12)') '  Fractional tolerance in the function value: ftol = ', ftol
    WRITE(out,'(a,g20.12)') '  Threshold for values of the fitted function: threshold = ', threshold
    WRITE(out,*)

  END SUBROUTINE PrintFit

  FUNCTION nDimFit() RESULT(nDim)
    ! Count number of free variables in the fit

    USE disloc_elasticity_ani, ONLY : nd
    USE dDipoleModule, ONLY : nDDipole
    IMPLICIT NONE
    INTEGER :: nDim

    nDim = Count(fit_cDislo(:,1:nd)) + Count(fit_bDislo(:,1:nd)) &
        + Count(fit_c1DDipole(:,1:nDDipole)) + Count(fit_c2DDipole(:,1:nDDipole)) &
        + Count(fit_bDDipole(:,1:nDDipole))

  END FUNCTION nDimFit

  SUBROUTINE Fit_PrintX(x, nDim, out)

    IMPLICIT NONE
    INTEGER, intent(in) :: nDim
    REAL(kind(0.d0)), dimension(1:nDim), intent(in) :: x
    INTEGER, intent(in) :: out

    INTEGER :: n

    WRITE(out,'(a,i0)') '  nDim = ', nDim
    DO n=1, nDim
       WRITE(out,'(a,i0,a,g20.12)') '  x(',n,') = ', x(n)
    END DO

  END SUBROUTINE Fit_PrintX

  SUBROUTINE Fit_WriteX(x, nDim)

    USE disloc_elasticity_ani, ONLY : nd, cDislo, bDislo
    USE dDipoleModule, ONLY : nDDipole, c1DDipole, c2DDipole, bDDipole
    IMPLICIT NONE
    INTEGER, intent(in) :: nDim
    REAL(kind(0.d0)), dimension(1:nDim), intent(out) :: x

    INTEGER :: nx, n, i

    nx = 0

    ! Dislocation positions
    DO n=1, nd  ! Loop on dislocation
       DO i=1, 3        ! Loop on coordinates
          IF (fit_cDislo(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_WriteX >'
                  END IF
                  x(nx) = cDislo(i,n)   ! Initialize to dislocation position
          END IF
       END DO
    END DO
    
    ! Dislocation Burgers vectors
    DO n=1, nd  ! Loop on dislocation
       DO i=1, 3        ! Loop on coordinates
          IF (fit_bDislo(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_WriteX >'
                  END IF
                  x(nx) = bDislo(i,n)   ! Initialize to dislocation Burgers vector
          END IF
       END DO
    END DO

    ! Dislocation dipole: position of center 1
    DO n=1, nDDipole  ! Loop on dipoles
       DO i=1, 3        ! Loop on coordinates
          IF (fit_c1DDipole(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_WriteX >'
                  END IF
                  x(nx) = c1DDipole(i,n)   ! Initialize to dislocation 1 position
          END IF
       END DO
    END DO

    ! Dislocation dipole: position of center 2
    DO n=1, nDDipole  ! Loop on dipoles
       DO i=1, 3        ! Loop on coordinates
          IF (fit_c2DDipole(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_WriteX >'
                  END IF
                  x(nx) = c2DDipole(i,n)   ! Initialize to dislocation 1 position
          END IF
       END DO
    END DO

    ! Dislocation dipole: Burgers vector
    DO n=1, nDDipole  ! Loop on dipoles
       DO i=1, 3        ! Loop on coordinates
          IF (fit_bDDipole(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_WriteX >'
                  END IF
                  x(nx) = bDDipole(i,n)   ! Initialize to dislocation 1 position
          END IF
       END DO
    END DO

  END SUBROUTINE Fit_WriteX

  SUBROUTINE Fit_ReadX(x, nDim)

    USE disloc_elasticity_ani, ONLY : nd, cDislo, bDislo
    USE dDipoleModule, ONLY : nDDipole, c1DDipole, c2DDipole, bDDipole
    USE rearrange
    IMPLICIT NONE
    INTEGER, intent(in) :: nDim
    REAL(kind(0.d0)), dimension(1:nDim), intent(in) :: x

    INTEGER :: nx, n, i, n2

    nx = 0

    ! Dislocation positions
    DO n=1, nd  ! Loop on dislocation
       DO i=1, 3        ! Loop on coordinates
          IF (fit_cDislo(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_ReadX >'
                  END IF
                  cDislo(i,n) = x(nx)   ! Initialize dislocation position
          END IF
       END DO
    END DO
    
    ! Dislocation Burgers vectors
    DO n=1, nd  ! Loop on dislocation
       DO i=1, 3        ! Loop on coordinates
          IF (fit_bDislo(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_ReadX >'
                  END IF
                  bDislo(i,n) = x(nx)   ! Initialize dislocation Burgers vector
                  ! Initialize Burgers vector of second dislocation composing the dipole
                  n2 = iDisloDipole(n)
                  IF (n2.GT.0) bDislo(i,n2) = -x(nx)
          END IF
       END DO
    END DO

    ! Dislocation dipole: position of center 1
    DO n=1, nDDipole  ! Loop on dipoles
       DO i=1, 3        ! Loop on coordinates
          IF (fit_c1DDipole(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_ReadX >'
                  END IF
                  c1DDipole(i,n) = x(nx)    ! Initialize dislocation 1 position
          END IF
       END DO
    END DO

    ! Dislocation dipole: position of center 2
    DO n=1, nDDipole  ! Loop on dipoles
       DO i=1, 3        ! Loop on coordinates
          IF (fit_c2DDipole(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_ReadX >'
                  END IF
                  c2DDipole(i,n) =x(nx)    ! Initialize dislocation 1 position
          END IF
       END DO
    END DO

    ! Dislocation dipole: Burgers vector
    DO n=1, nDDipole  ! Loop on dipoles
       DO i=1, 3        ! Loop on coordinates
          IF (fit_bDDipole(i,n)) THEN
                  nx=nx+1
                  IF (nx.GT.nDim) THEN
                          WRITE(0,'(a,i0)') 'Problem with size of X array: nDim = ', nDim
                          STOP '< Fit_ReadX >'
                  END IF
                  bDDipole(i,n) =x(nx)   ! Initialize dislocation 1 position
          END IF
       END DO
    END DO
  END SUBROUTINE Fit_ReadX

END MODULE FitModule
