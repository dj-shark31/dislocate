MODULE SymmetryXmlModule

        USE Symmetry3DModule
        IMPLICIT NONE

        LOGICAL, parameter, private :: debug=.false. 

        character(len=80), private :: tag
        logical, private :: endtag
        character(len=80), dimension(1:2,1:20), private :: attribs
        integer, private :: noattribs
        character(len=200), dimension(1:20), private :: data
        integer, private :: nodata

        PRIVATE :: ReadIRotXml, ReadPerioVectorsXml, mat3inv

CONTAINS

        SUBROUTINE ReadSymmetryFileXml(xmlFile, symGroup)

                USE XmlModule
                IMPLICIT NONE
                CHARACTER(len=*), intent(in) :: xmlFile
                TYPE(sym3D_group_t), intent(out) :: symGroup                

                Type(Xml_Parse) :: parser
                INTEGER :: n, nSymRead
                REAL(kind(0.d0)) :: alat 
                REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_uPBC

                ! Initialization
                CALL InitSym3Dgroup(symGroup)
                nSymRead = 0
                alat = 1.d0
                
                ! Open xml file
                CALL Xml_Open( parser, xmlFile, mustread=.true. )

                ! Read xml file
                DO 

                   CALL xml_get( parser, tag, endtag, attribs, noattribs, data, nodata )

                   IF (.NOT.endTag) THEN
                           IF (debug) CALL XmlPrint(tag, endtag, attribs, noattribs, data, nodata, 6)
                           SELECT CASE(Trim(tag))
                           CASE('symmetry')
                           CASE('nSym')
                                   symGroup%nSym = XmlReadData_int(data, nodata)
                           CASE('irot') 
                                   nSymRead = nSymRead + 1
                                   CALL ReadIRotXml(symGroup%op(nSymRead)%iRot(:,:), symGroup%op(nSymRead)%uReduced(:), parser)
                           CASE('periodicity_vectors')
                                   CALL ReadPerioVectorsXml(symGroup%uPBC(:,:), parser)
                           CASE('distance_scaling')
                                   alat = XmlReadData_real(data, nodata)
                           CASE DEFAULT
                                   IF (.NOT.XmlCommentTag(tag)) THEN
                                           WRITE(6,'(a)') 'WARNING <ReadStructFileXml> : unknown xml tag'
                                           CALL XmlPrint(tag, endtag, attribs, noattribs, data, nodata, 6)
                                   END IF
                           END SELECT

                   ELSE IF (tag.EQ.'symmetry') THEN
                           EXIT
                   ELSE
                           CYCLE
                   END IF

                END DO
                CALL Xml_Close( parser )

                ! Check number of symmetry operations
                IF (nSymRead.NE.symGroup%nSym) THEN
                        IF (symGroup%nSym.LE.0) THEN
                                symGroup%nSym=nSymRead
                        ELSE
                                WRITE(0,'(a,i0)') "  number of symmetry operations read:    ", nSymRead
                                WRITE(0,'(a,i0)') "  number of symmetry operations defined: ", symGroup%nSym
                                STOP "< ReadSymmetryFileXml >"
                        END IF
                END IF

                ! Multiply by lattice parameter
                symGroup%uPBC(:,:) = alat*symGroup%uPBC(:,:)

                ! Translation in cartesian coordinates for each symmetry operations
                DO n=1, symGroup%nSym
                   symGroup%op(n)%u(1:3) = MatMul( symGroup%uPBC(:,:), symGroup%op(n)%uReduced(:) )
                END DO

                ! Rotation matrix in cartesian coordinates
                CALL Mat3Inv( symGroup%uPBC(:,:), inv_uPBC(:,:) )
                DO n=1, symGroup%nSym
                   symGroup%op(n)%rot(:,:) = MatMul( symGroup%uPBC(:,:) , &
                           MatMul( Dble( symGroup%op(n)%iRot(:,:) ),  inv_uPBC(:,:) ) )
                END DO

                ! Check compatibility between rotation matrixes and periodicity vectors
                CALL CheckSym3Dgroup(symGroup, 0)

        END SUBROUTINE ReadSymmetryFileXml

        SUBROUTINE ReadIRotXml(iRot, uReduced, parser)

                USE XmlModule
                IMPLICIT NONE

                ! Rotation matrix and translation vector in reduced coordinates
                INTEGER, dimension(1:3,1:3), intent(out) :: iRot
                REAL(kind(0.d0)), dimension(1:3), intent(out) :: uReduced
                Type(Xml_Parse) :: parser

                ! Initialization
                iRot(:,:) = 0
                uReduced(:) = 0.d0

                DO 
                   CALL xml_get( parser, tag, endtag, attribs, noattribs, data, nodata )
                   IF (.NOT.endTag) THEN
                           IF (debug) CALL XmlPrint(tag, endtag, attribs, noattribs, data, nodata, 6)
                           SELECT CASE(Trim(tag))
                           CASE('irot')
                           CASE('isymop')
                                   iRot(1:3,1:3) = XmlReadData_int_matrix(data, nodata, 3, 3)
                           CASE('gtrans')
                                   uReduced(1:3) = XmlReadData_real_vector(data, nodata, 3)
                           CASE DEFAULT
                                   IF (.NOT.XmlCommentTag(tag)) THEN
                                           WRITE(6,'(a)') 'WARNING <ReadIRotXml> : unknown xml tag'
                                           CALL XmlPrint(tag, endtag, attribs, noattribs, data, nodata, 6)
                                   END IF
                           END SELECT

                   ELSE IF (tag.EQ.'irot') THEN
                           EXIT
                   ELSE
                           CYCLE
                   END IF
                END DO

        END SUBROUTINE ReadIRotXml

        SUBROUTINE ReadPerioVectorsXml(at, parser)

                USE XmlModule
                IMPLICIT NONE

                REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
                Type(Xml_Parse) :: parser

                ! Initialization
                at(:,:) = 0.d0

                DO 
                   CALL xml_get( parser, tag, endtag, attribs, noattribs, data, nodata )
                   IF (.NOT.endTag) THEN
                           IF (debug) CALL XmlPrint(tag, endtag, attribs, noattribs, data, nodata, 6)
                           SELECT CASE(Trim(tag))
                           CASE('periodicity_vectors')
                           CASE('a1')
                                   at(1:3,1) = XmlReadData_real_vector(data, nodata, 3)
                           CASE('a2')
                                   at(1:3,2) = XmlReadData_real_vector(data, nodata, 3)
                           CASE('a3')
                                   at(1:3,3) = XmlReadData_real_vector(data, nodata, 3)
                           CASE DEFAULT
                                   IF (.NOT.XmlCommentTag(tag)) THEN
                                           WRITE(6,'(a)') 'WARNING <ReadPerioVecorsXml> : unknown xml tag'
                                           CALL XmlPrint(tag, endtag, attribs, noattribs, data, nodata, 6)
                                   END IF
                           END SELECT

                   ELSE IF (tag.EQ.'periodicity_vectors') THEN
                           EXIT
                   ELSE
                           CYCLE
                   END IF
                END DO

        END SUBROUTINE ReadPerioVectorsXml

        SUBROUTINE Mat3Inv(A, B)
          ! return matrix B which is the inverse of matrix A
          ! A and B are 3x3 matrices

          implicit none

          REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
          REAL(kind(0.d0)), dimension(3,3), intent(out) :: B

          !REAL(kind(0.d0)), dimension(3,3) :: C ! DEBUG

          REAL(kind(0.d0)) :: invdet

             
          b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
          b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
          b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
             
          b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
          b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
          b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
             
          b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
          b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
          b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

          invdet = 1.d0/( a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1) )
          b(:,:)=b(:,:)*invdet

          !!!!! DEBUG !!!!
          !c = MatMul(a,b)
          !WRITE(6,'(3g14.7)') C(1,1:3)
          !WRITE(6,'(3g14.7)') C(2,1:3)
          !WRITE(6,'(3g14.7)') C(3,1:3)
          !STOP

        END SUBROUTINE mat3inv

END MODULE SymmetryXmlModule
