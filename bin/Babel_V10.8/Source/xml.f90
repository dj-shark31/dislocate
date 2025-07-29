MODULE XmlModule

        USE XmlParse
        
CONTAINS

        FUNCTION XmlCommentTag(tag) RESULT(test)

                IMPLICIT NONE
                character(len=*), intent(in)                 :: tag
                LOGICAL :: test

                IF (tag(1:3).EQ.'!--') THEN
                        test=.true.
                ELSE
                        test=.false.
                END IF

        END FUNCTION XmlCommentTag

        FUNCTION XmlReadData_int(data, nodata) RESULT(iValue)

                IMPLICIT NONE
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                INTEGER :: iValue

                CHARACTER(len=nodata*(len(data(1))+1)) :: full_data
                INTEGER :: n

                full_data = data(1)
                DO n=2, nodata
                   full_data = Trim(full_data) // " " // data(n)
                END DO
                READ(full_data,*) iValue
                
        END FUNCTION XmlReadData_int

        FUNCTION XmlReadData_int_vector(data, nodata, nDim) RESULT(iValue)

                IMPLICIT NONE
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                INTEGER, intent(in) :: nDim
                INTEGER, dimension(1:nDim) :: iValue

                CHARACTER(len=nodata*(len(data(1))+1)) :: full_data
                INTEGER :: n

                full_data = data(1)
                DO n=2, nodata
                   full_data = Trim(full_data) // " " // data(n)
                END DO
                READ(full_data,*) iValue(1:nDim)
                
        END FUNCTION XmlReadData_int_vector

        FUNCTION XmlReadData_int_matrix(data, nodata, nDim, mDim) RESULT(iValue)

                IMPLICIT NONE
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                INTEGER, intent(in) :: nDim, mDim
                INTEGER, dimension(1:nDim, 1:mDim) :: iValue

                CHARACTER(len=nodata*(len(data(1))+1)) :: full_data
                INTEGER :: n

                full_data = data(1)
                DO n=2, nodata
                   full_data = Trim(full_data) // " " // data(n)
                END DO
                READ(full_data,*) ( iValue(n, 1:mDim), n=1, nDim )
                
        END FUNCTION XmlReadData_int_matrix

        FUNCTION XmlReadData_real(data, nodata) RESULT(rValue)

                IMPLICIT NONE
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                REAL(kind(0.d0)) :: rValue

                CHARACTER(len=nodata*(len(data(1))+1)) :: full_data
                INTEGER :: n

                full_data = data(1)
                DO n=2, nodata
                   full_data = Trim(full_data) // " " // data(n)
                END DO
                READ(full_data,*) rValue
                
        END FUNCTION XmlReadData_real

        FUNCTION XmlReadData_real_vector(data, nodata, nDim) RESULT(rValue)

                IMPLICIT NONE
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                INTEGER, intent(in) :: nDim
                REAL(kind(0.d0)), dimension(1:nDim) :: rValue

                CHARACTER(len=nodata*(len(data(1))+1)) :: full_data
                INTEGER :: n

                full_data = data(1)
                DO n=2, nodata
                   full_data = Trim(full_data) // " " // data(n)
                END DO
                READ(full_data,*) rValue(1:nDim)
                
        END FUNCTION XmlReadData_real_vector

        FUNCTION XmlReadData_real_matrix(data, nodata, nDim, mDim) RESULT(rValue)

                IMPLICIT NONE
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                INTEGER, intent(in) :: nDim, mDim
                REAL(kind(0.d0)), dimension(1:nDim, 1:mDim) :: rValue

                CHARACTER(len=nodata*(len(data(1))+1)) :: full_data
                INTEGER :: n

                full_data = data(1)
                DO n=2, nodata
                   full_data = Trim(full_data) // " " // data(n)
                END DO
                READ(full_data,*) ( rValue(n,1:mDim), n=1, nDim )
                
        END FUNCTION XmlReadData_real_matrix

        SUBROUTINE XmlPrint(tag, endtag, attribs, noattribs, data, nodata, out)

                IMPLICIT NONE
                character(len=*), intent(in)                 :: tag
                logical, intent(in)                          :: endtag
                character(len=*), dimension(:,:), intent(in) :: attribs
                integer, intent(in)                          :: noattribs
                character(len=*), dimension(:), intent(in)   :: data
                integer, intent(in)                          :: nodata
                INTEGER, intent(in) :: out

                INTEGER :: n

                IF (endtag) THEN
                        WRITE(out,'(2(a,1x))')   'end tag: ', Trim(tag)
                        WRITE(out,*)
                ELSE
                        WRITE(out,'(2(a,1x))')   'tag: ', Trim(tag)
                        WRITE(out,'(a,i0)') 'noattribs: ', noattribs
                        DO n=1, noattribs
                           WRITE(out,'(a,i0,a,2(a,1x))')   'attribs ', n,': ', Trim(attribs(1,n)), Trim(attribs(2,n))
                        END DO
                        WRITE(out,'(a,i0)') 'nodata: ', nodata
                        DO n=1, nodata
                           WRITE(out,'(a,i0,2a)')   'data', n, ': ', Trim(data(N))
                        END DO
                        WRITE(out,*)
                END IF

        END SUBROUTINE XmlPrint

END MODULE XmlModule

