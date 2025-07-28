MODULE Associate_lecture

  ! Name of the input file defining the dislocation
  CHARACTER(len=100) :: input_file

CONTAINS

  SUBROUTINE Read_Associate(out)

    USE Associate_data
    USE Structure_module
    IMPLICIT NONE
    INTEGER, intent(in) :: out

    ! Format of the input and reference structure files
    LOGICAL :: inpXyz, inpCfg, inpGin, inpLisa, inpSiesta
    LOGICAL :: refXyz, refCfg, refGin, refLisa, refSiesta
    ! Input and refrence structure files
    CHARACTER(len=100) :: refFile, inpFile
    NAMELIST /input/ inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpFile, &
          refXyz, refCfg, refGin, refLisa, refSiesta, refFile, &
          outXyz, outCfg, outGin, outLisa, outSiesta, outFile, &
          periodic, &
          out_sort, out_index, out_alat, nTypes, label, mass, & 
          verbosity

    INTEGER :: input_unit, i
    CHARACTER(len=6) :: fileType


    !===========================================
    ! Initialization
    verbosity=4
    inpXyz=.FALSE. ; inpCfg=.FALSE. ; inpGin=.FALSE. ; inpLisa=.FALSE. ; inpSiesta=.FALSE. 
    refXyz=.FALSE. ; refCfg=.FALSE. ; refGin=.FALSE. ; refLisa=.FALSE. ; refSiesta=.FALSE. 
    outXyz=.FALSE. ; outCfg=.FALSE. ; outGin=.FALSE. ; outLisa=.FALSE. ; outSiesta=.FALSE. 
    periodic=.FALSE.
    out_sort=.TRUE. ; out_index=.TRUE.
    inpFile='-' ; refFile='-' ; outFile='-'
    imm=0
    im=0 ; im0=0
    at(:,:)=0.d0 ; at0(:,:)=0.d0
    CALL InitStructure()       ! Initialize out_alat, nTypes, mass(:), label(:)
   

    !===========================================
    ! Read input data
    IF (input_file.EQ.'-') THEN
            input_unit=5
    ELSE
            input_unit=50
            OPEN(file=input_file,unit=input_unit,action='read',status='old')
    END IF
    READ(input_unit,nml=input)
    IF (input_unit.NE.5) CLOSE(input_unit)


    !===========================================
    ! Read input structure file
    IF ( Count( (/inpXyz,inpCfg,inpGin,inpLisa,inpSiesta/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa, and Siesta formats&
                & for input structure'
            STOP '< Read_Associate >'
    END IF
    IF (inpXyz.OR.inpCfg.OR.inpGin.OR.inpLisa.OR.inpSiesta) THEN
            IF (inpXyz) fileType='xyz'
            IF (inpCfg) fileType='cfg'
            IF (inpGin) fileType='gin'
            IF (inpSiesta) fileType='siesta'
            IF (inpLisa) fileType='lisa'
            CALL ReadStructure(xp, iTyp, im, imm, at, inpFile, fileType)
    END IF


    !===========================================
    ! Read reference structure file
    IF ( Count( (/refXyz,refCfg,refGin,refLisa,refSiesta/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa, and Siesta formats&
                & for refut structure'
            STOP '< Read_Associate >'
    END IF
    IF (refXyz.OR.refCfg.OR.refGin.OR.refLisa.OR.refSiesta) THEN
            IF (refXyz) fileType='xyz'
            IF (refCfg) fileType='cfg'
            IF (refGin) fileType='gin'
            IF (refSiesta) fileType='siesta'
            IF (refLisa) fileType='lisa'
            CALL ReadStructure(xp0, iTyp0, im0, imm, at0, refFile, fileType)
    END IF

    ! Check number of atoms
    IF (im0.NE.im) THEN
            WRITE(0,'(a)') 'Both input and reference structure files do not contain the&
                  & same number of atoms'
            WRITE(0,'(a,i0)') '  im0 = ', im0
            WRITE(0,'(a,i0)') '  im  = ', im 
            STOP "< Read_Associate >"
    END IF

   ! Check that all atoms are of the same type
   DO i=1, im
      IF (iTyp0(i).nE.iTyp(i)) THEN
              WRITE(0,'(a,i0,a)') 'Atom ', i, ' is not of the same type in both&
                & input and reference structure files'
              WRITE(0,'(a,i0)') '  type in input file: ', iTyp(i)
              WRITE(0,'(a,i0)') '  type in reference file: ', iTyp0(i)
              STOP '< Read_Associate >'
      END IF
   END DO

  END SUBROUTINE Read_Associate

END MODULE Associate_lecture
