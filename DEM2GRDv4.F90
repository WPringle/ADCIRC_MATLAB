! COPYRIGHT (C) 2013 MATT BILSKIE / UNIVERSITY OF CENTRAL FLORIDA
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

!	DEVELOPED BY:
!								
!	MATTHEW V. BILSKIE, E.I. (Matt.Bilske@gmail.com)
!       UNIVESITY OF CENTRAL FLORIDA
!	CHAMPS LAB (CHAMPSLAB.COM)
!
!   PLEASE EMAIL WITH QUESTIONS/BUGS

!   PORTIONS OF THIS CODE UTILIZE ROUTINES FOUND IN CRYSTAL FULCHERS SURFACE
!   ROUGHNESS CODES - LINES 298 & 405 
!   (MANNINGS_N_FINDER, SURFACE_ROUGHNESS_CALC, AND SURFACE_CANOPY)
!   ALONG WITH USING FORTRAN TYPES THANKS TO ZACH COBELL

!   CITATION (PLEASE CITE THE FOLLOWING PUBLICATION):
!
!   Matthew V. Bilskie, Scott C. Hagen, Topographic accuracy assessment of bare
!   earth lidar-derived unstructured meshes, Advances in Water Resources,
!   Volume 52, February 2013, Pages 165-177

    MODULE MISC
        CONTAINS
        LOGICAL FUNCTION FileExist(fileName,flag)
            CHARACTER(*),INTENT(IN)       :: fileName
            INTEGER,INTENT(IN)            :: flag
            LOGICAL                       :: fExist
            
            INQUIRE(FILE=TRIM(fileName),EXIST=fExist)
            IF(.NOT.fExist)THEN
                WRITE(*,'(A)')'File '//TRIM(fileName)//' not found.'
                FileExist = .FALSE.
                IF (flag.eq.0) THEN
                   STOP
                ELSE
                   WRITE(*,'(A)') 'Continue setting flag = 1 everywhere'
                ENDIF
            ELSE
                FileExist = .TRUE.
            ENDIF
        END FUNCTION

        REAL(8) FUNCTION Distance(lon1,lat1,lon2,lat2,coord)
            IMPLICIT NONE
            REAL(8),INTENT(IN)          :: lon1
            REAL(8),INTENT(IN)          :: lat1
            REAL(8),INTENT(IN)          :: lon2
            REAL(8),INTENT(IN)          :: lat2
            LOGICAL,INTENT(IN)          :: coord

            REAL(8)                     :: a
            REAL(8)                     :: c
            REAL(8),PARAMETER           :: R = 6371.64D0 !km
            REAL(8)                     :: deg2rad
            REAL(8),PARAMETER           :: PI = 3.14159265359D0
            REAL(8)                     :: dlon
            REAL(8)                     :: dlat

            INTRINSIC                   :: SQRT
            INTRINSIC                   :: DSIN
            INTRINSIC                   :: DCOS
            INTRINSIC                   :: DATAN2

            deg2rad = PI / 180.0D0

            IF(coord)THEN !...If geographic coordinates
                ! Haversine Formula (www.movable-type.co.uk_scripts_latlong.html)
                dlon = (lon2 - lon1) * deg2rad
                dlat = (lat2 - lat1) * deg2rad
                a = DSIN(dlat/2.0D0)**2 + DCOS(lat1*deg2rad)*DCOS(lat2*deg2rad) &
                    * DSIN(dlon/2.0D0)**2
                c = 2.D0 * DATAN2(SQRT(a),SQRT(1.D0-a))
                Distance = R * c * 1000 !...Convert km to m
            ELSE !...Cartesian coordinates
                Distance = SQRT((lon2 - lon1)**2 + (lat2 - lat1)**2)
            ENDIF
                
        END FUNCTION
    END MODULE MISC
    
    MODULE ADCGRID
        USE MISC

        TYPE BoundaryListing
            INTEGER                     :: NumNodes
            INTEGER                     :: Code
            INTEGER,ALLOCATABLE         :: N1(:)
            INTEGER,ALLOCATABLE         :: N2(:)
            REAL(8),ALLOCATABLE         :: Crest(:)
            REAL(8),ALLOCATABLE         :: Supercritical(:)
            REAL(8),ALLOCATABLE         :: Subcritical(:)
        END TYPE
    
        TYPE MESH
            CHARACTER(200)              :: title
            INTEGER                     :: ne
            INTEGER                     :: np
            REAL(8),ALLOCATABLE         :: nodes(:,:)
            INTEGER,ALLOCATABLE         :: elem(:,:)
            INTEGER                     :: NumOpenBoundaries
            INTEGER                     :: TotNumOpenBoundaryNodes
            INTEGER                     :: NumLandBoundaries
            INTEGER                     :: TotNumLandBoundaryNodes
            TYPE(BoundaryListing),ALLOCATABLE   :: OceanBC(:)
            TYPE(BoundaryListing),ALLOCATABLE   :: LandBC(:)
        END TYPE MESH
        
    CONTAINS
    
    SUBROUTINE ReadMesh(meshName,myMesh)
        IMPLICIT NONE
        CHARACTER(200),INTENT(IN)       :: meshName
        TYPE(MESH),INTENT(OUT)          :: myMesh
        INTEGER                         :: I,J
        INTEGER                         :: tempI
        INTEGER                         :: tempI2
        LOGICAL                         :: exists
        
        exists = FileExist(meshName,0)
        
        WRITE(*,'(A)')'Reading mesh...'
        OPEN(UNIT=14,FILE=TRIM(meshName),ACTION='READ')
        READ(14,*)myMesh%title
        READ(14,*)myMesh%ne,myMesh%np
        ALLOCATE(myMesh%nodes(1:myMesh%np,1:3))
        ALLOCATE(myMesh%elem(1:myMesh%ne,1:3))
        DO I = 1, myMesh%np
            READ(14,*)tempI,myMesh%nodes(I,1),myMesh%nodes(I,2), & 
                myMesh%nodes(I,3)
            IF(tempI.NE.I)THEN
                WRITE(*,'(A)')'Mesh needs renumbered...'
                WRITE(*,'(A)')''
                STOP
            ENDIF
        ENDDO
        DO I = 1, myMesh%ne
            READ(14,*)tempI,tempI2,myMesh%elem(I,1),myMesh%elem(I,2), &
                myMesh%elem(I,3)
            IF(tempI.NE.I)THEN
                WRITE(*,'(A)')'Mesh needs renumbered...'
                WRITE(*,'(A)')''
                    STOP
            ENDIF
        ENDDO

        !...Read Open BCs
        READ(14,*,END=100) myMesh%NumOpenBoundaries
        READ(14,*,END=100) myMesh%TotNumOpenBoundaryNodes
        ALLOCATE(myMesh%OceanBC(MyMesh%NumOpenBoundaries))
        DO I = 1, myMesh%NumOpenBoundaries
            READ(14,*) myMesh%OceanBC(I)%NumNodes
            ALLOCATE(myMesh%OceanBC(I)%N1(1:myMesh%OceanBC(I)%NumNodes))
            DO J = 1, myMesh%OceanBC(I)%NumNodes
                READ(14,*), myMesh%OceanBC(I)%N1(J)
            ENDDO
        ENDDO

        !...Read Land BCs
        READ(14,*,END=200) myMesh%NumLandBoundaries
        READ(14,*,END=200) myMesh%TotNumLandBoundaryNodes
        ALLOCATE(myMesh%LandBC(myMesh%NumLandBoundaries))
        DO I = 1, myMesh%NumLandBoundaries
            ReAD(14,*) myMesh%LandBC(I)%NumNodes,&
                myMesh%LandBC(I)%Code
            SELECT CASE(myMesh%LandBC(I)%Code)
                CASE(0,1,10,11,12,20,21,22,52)
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*) myMesh%LandBC(I)%N1(J)
                    ENDDO
                CASE(13,23)
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Crest(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Supercritical(1:myMesh%LandBC(I)%NumNodes))
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*) myMesh%LandBC(I)%N1(J), &
                            myMesh%LandBC(I)%Crest(J), &
                            myMesh%LandBC(I)%Supercritical(J)
                    ENDDO
                CASE(24)
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%N2(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Crest(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Supercritical(1:myMesh%LandBC(I)%NumNodes))
                    ALLOCATE(myMesh%LandBC(I)%Subcritical(1:myMesh%LandBC(I)%NumNodes))
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*)myMesh%LandBC(I)%N1(J),&
                            myMesh%LandBC(I)%N2(J), &
                            myMesh%LandBC(I)%Crest(J), &
                            myMesh%LandBC(I)%Subcritical(J), &
                            myMesh%LandBC(I)%Supercritical(J)
                    ENDDO
                CASE DEFAULT
                    ALLOCATE(myMesh%LandBC(I)%N1(1:myMesh%LandBC(I)%NumNodes))
                    WRITE(*,'(2A,I0)') 'WARNING: Unknown boundary ',&
                        'condition. ADCIRC TYPE = ',&
                        myMesh%LandBC(I)%Code
                    WRITE(*,'(A)') '        READ AS A SINGLE NODE.'
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        READ(14,*) myMesh%LandBC(I)%N1(J)
                    ENDDO
            END SELECT
        ENDDO

        WRITE(*,'(A)')'done!'
        CLOSE(14)

        RETURN
        
100     CONTINUE
        myMesh%NumOpenBoundaries = 0
        myMesh%TotNumOpenBoundaryNodes = 0
200     CONTINUE
        myMesh%NumLandBoundaries = 0
        myMesh%TotNumLandBoundaryNodes = 0
        RETURN

    END SUBROUTINE ReadMesh

    SUBROUTINE WriteMesh(myMesh,fileName,coord)
        IMPLICIT NONE
        TYPE(MESH),INTENT(IN)           :: myMesh
        CHARACTER(200),INTENT(IN)       :: fileName
        LOGICAL,INTENT(IN)              :: coord

        INTEGER                         :: I,J
        
        WRITE(*,'(A)')'Writing new mesh...'
        OPEN(UNIT=14,FILE=TRIM(fileName),ACTION='WRITE')
        WRITE(14,*)'Interpolated from DEM2GRD.F90'
        WRITE(14,'(I9,x,I9)')myMesh%ne,myMesh%np
        DO I = 1, myMesh%np
            IF(.NOT.coord)THEN !...Cartesian
                WRITE(14,'(I9,x,F15.7,5x,F15.7,5x,F15.7)')I,myMesh%nodes(I,1), &
                myMesh%nodes(I,2),myMesh%nodes(I,3)
            ELSE
                WRITE(14,'(I9,x,F14.10,5x,F14.10,5x,F15.7)')I,myMesh%nodes(I,1), &
                myMesh%nodes(I,2),myMesh%nodes(I,3)
            ENDIF
        ENDDO
        DO I = 1, myMesh%ne
            WRITE(14,*)I,'3',myMesh%elem(I,1),myMesh%elem(I,2),myMesh%elem(I,3)  
        ENDDO
        WRITE(14,'(I10)') myMesh%NumOpenBoundaries
        WRITE(14,'(I10)') myMesh%TotNumOpenBoundaryNodes
        DO I = 1, myMesh%NumOpenBoundaries
            WRITE(14,'(I10)') myMesh%OceanBC(I)%NumNodes
            DO J = 1, myMesh%OceanBC(I)%NumNodes
                WRITE(14,*) myMesh%OceanBC(I)%N1(J)
            ENDDO
        ENDDO
        WRITE(14,*) myMesh%NumLandBoundaries
        WRITE(14,*) myMesh%TotNumLandBoundaryNodes
        DO I = 1, myMesh%NumLandBoundaries
            WRITE(14,'(I6,4x,I6,4x,A,I6)') myMesh%LandBC(I)%NumNodes, &
                myMesh%LandBC(I)%Code,"!=seg",I
            SELECT CASE(myMesh%LandBC(I)%Code)
                CASE(0,1,10,11,12,20,21,22,52)
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10)') myMesh%LandBC(I)%N1(J)
                    ENDDO
                CASE(13,23)
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10,2x,F16.3,2x,F16.3)') myMesh%LandBC(I)%N1(J), &
                            myMesh%LandBC(I)%Crest(J),&
                            myMesh%LandBC(I)%Supercritical(J)
                    ENDDO
                CASE(24)
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10,2x,I10,2x,F16.3,2x,F16.3,2x,F16.3,2x,F16.3)') &
                            myMesh%LandBC(I)%N1(J), &
                            myMesh%LandBC(I)%N2(J), &
                            myMesh%LandBC(I)%Crest(J),&
                            myMesh%LandBC(I)%Supercritical(J), &
                            myMesh%LandBC(I)%Subcritical(J)
                    ENDDO
                CASE DEFAULT
                    DO J = 1, myMesh%LandBC(I)%NumNodes
                        WRITE(14,'(I10)') myMesh%LandBC(I)%N1(J)
                    ENDDO
            END SELECT
        ENDDO

        CLOSE(14)
        WRITE(*,'(A)')'done!'

        RETURN

    END SUBROUTINE WriteMesh

    SUBROUTINE ElementSize(myMesh,coord,esize)
        USE MISC
        IMPLICIT NONE
        TYPE(MESH),INTENT(IN)           :: myMesh
        LOGICAL,INTENT(IN)              :: coord
        REAL(8),ALLOCATABLE,INTENT(OUT) :: esize(:)
        
        INTEGER                         :: I
        REAL(8)                         :: DIST(3)

        ALLOCATE(esize(1:myMesh%np))
        WRITE(*,'(A)')'Computing element size...'
        DO I = 1, myMesh%ne
            DIST(1) = Distance(myMesh%nodes(myMesh%elem(I,1),1), &
                myMesh%nodes(myMesh%elem(I,1),2), &
                myMesh%nodes(myMesh%elem(I,2),1), &
                myMesh%nodes(myMesh%elem(I,2),2), coord)
            DIST(2) = Distance(myMesh%nodes(myMesh%elem(I,2),1), &
                myMesh%nodes(myMesh%elem(I,2),2), &
                myMesh%nodes(myMesh%elem(I,3),1), &
                myMesh%nodes(myMesh%elem(I,3),2), coord)
            DIST(3) = Distance(myMesh%nodes(myMesh%elem(I,3),1), &
                myMesh%nodes(myMesh%elem(I,3),2), &
                myMesh%nodes(myMesh%elem(I,1),1), &
                myMesh%nodes(myMesh%elem(I,1),2), coord)
            IF(DIST(1).GT.eSize(MyMesh%elem(I,1)))THEN
                eSize(MyMesh%elem(I,1)) = DIST(1)
            ENDIF
            IF(DIST(2).GT.eSize(MyMesh%elem(I,2)))THEN
                eSize(MyMesh%elem(I,2)) = DIST(2)
            ENDIF
            IF(DIST(3).GT.eSize(MyMesh%elem(I,3)))THEN
                eSize(MyMesh%elem(I,3)) = DIST(3)
            ENDIF

        ENDDO
#ifdef WRITEMESHSIZE
        OPEN(UNIT=14,FILE='ElementSize.grd',ACTION='WRITE')
        WRITE(14,*)myMesh%title
        WRITE(14,'(I9,x,I9)')myMesh%ne,myMesh%np
        DO I = 1, myMesh%np
            WRITE(14,*)I,myMesh%nodes(I,1),myMesh%nodes(I,2),esize(I)
        ENDDO
        DO I = 1, myMesh%ne
            WRITE(14,*)I,'3',myMesh%elem(I,1),myMesh%elem(I,2),myMesh%elem(I,3)  
        ENDDO
        CLOSE(14)
#endif
        WRITE(*,'(A)')'done!'
    END SUBROUTINE ElementSize
    
    END MODULE ADCGrid
    
    MODULE DEM
        USE MISC 
        USE NETCDF 
 
        TYPE DEM_G
            INTEGER                     :: NCOLS
            INTEGER                     :: NROWS
            REAL(8)                     :: x_botlt
            REAL(8)                     :: y_botlt
            REAL(8)                     :: x_y_dist
            REAL(8)                     :: nodata_value
            REAL(4),ALLOCATABLE         :: values(:,:)
        END TYPE DEM_G
        
    CONTAINS
        
        SUBROUTINE ReadDEM(DEMName,myDEM)
            IMPLICIT NONE
            CHARACTER(200),INTENT(IN)   :: DEMName
            TYPE(DEM_G),INTENT(OUT)     :: myDEM
            LOGICAL                     :: exists
            CHARACTER(10)               :: ndum, lat_name, lon_name, depth_name
            CHARACTER(10)               :: row_name, col_name
            INTEGER                     :: I
            INTEGER                     :: J
            INTEGER                     :: L, status, FlipUD
            INTEGER                     :: NC_ID, lat_ID, lon_ID, depth_ID
            REAL(4),ALLOCATABLE         :: DUMMY(:,:)
            REAL(8)                     :: t1, t2
            
            exists = FileExist(DEMName,0)
            L = LEN_TRIM(DEMName)
            IF (DEMName(L-2:L) == 'flt') THEN 
               WRITE(*,'(A)')'Reading raster (*.flt) file...'
               ! Read the header file
               exists = FileExist(DEMName(1:L-3)//'hdr',0)
               OPEN(UNIT=10,FILE=DEMName(1:L-3)//'hdr',ACTION='READ')
               READ(10,*) ndum, myDEM%NCOLS
               READ(10,*) ndum, myDEM%NROWS
               READ(10,*) ndum, myDEM%x_botlt
               READ(10,*) ndum, myDEM%y_botlt
               READ(10,*) ndum, myDEM%x_y_dist
               READ(10,*) ndum, myDEM%nodata_value
               CLOSE(10)
            ELSEIF (DEMName(L-1:L) == 'nc') THEN
               WRITE(*,'(A)')'Reading NETCDF (*.nc) file...'
               exists = FileExist(DEMName(1:L-2)//'hdr',0)
               OPEN(UNIT=10,FILE=DEMName(1:L-2)//'hdr',ACTION='READ')
               READ(10,*) row_name
               READ(10,*) col_name
               READ(10,*) lat_name
               READ(10,*) lon_name
               READ(10,*) depth_name
               READ(10,*) FlipUD
               CLOSE(10)
               !  Open NETCDF and get info about lon and lat dimensions etc.
               CALL Check(NF90_OPEN(TRIM(DEMName),NF90_NOWRITE,NC_ID))
               ! Get dimids of lat and lon
               CALL Check(NF90_INQ_DIMID(NC_ID,row_name,lat_ID))
               CALL Check(NF90_INQ_DIMID(NC_ID,col_name,lon_ID))
               ! Get the lengths of lat and lon
               CALL Check(NF90_INQUIRE_DIMENSION(NC_ID,lat_ID,len=myDEM%NROWS))
               CALL Check(NF90_INQUIRE_DIMENSION(NC_ID,lon_ID,len=myDEM%NCOLS))
               ! Get the varids of lat and lon
               CALL Check(NF90_INQ_VARID(NC_ID,lat_name,lat_ID))
               CALL Check(NF90_INQ_VARID(NC_ID,lon_name,lon_ID))
               ! Get the first entries of lat and lon
               CALL Check(NF90_GET_VAR(NC_ID,lat_ID,t1,start=[1,1]))
               CALL Check(NF90_GET_VAR(NC_ID,lat_ID,t2,start=[myDEM%NROWS,1]))
               myDEM%y_botlt = min(t1,t2)
               CALL Check(NF90_GET_VAR(NC_ID,lon_ID,t1,start=[1,1]))
               CALL Check(NF90_GET_VAR(NC_ID,lon_ID,t2,start=[myDEM%NCOLS,1]))
               myDEM%x_botlt = min(t1,t2)
               ! Get the last entry of lon and then minus off first dividing by len - 1 to get dx
               CALL Check(NF90_GET_VAR(NC_ID,lon_ID,myDEM%x_y_dist,start=[2,1]))
               myDEM%x_y_dist = myDEM%x_y_dist - myDEM%x_botlt
               ! Get varid of depth
               CALL Check(NF90_INQ_VARID(NC_ID,depth_name,depth_ID))
               ! Get the FillValue of the depth
               status = NF90_GET_ATT(NC_ID,depth_ID,'_FillValue',myDEM%nodata_value)
               IF (status /= nf90_noerr) THEN
                   myDEM%nodata_value = -9999d0  
               ENDIF
               write(*,*) myDEM%NROWS, myDEM%NCOLS, myDEM%x_y_dist
               write(*,*) myDEM%x_botlt, myDEM%y_botlt
            ELSE
               STOP 'Unrecognized filetype'
            ENDIF

            ALLOCATE(myDEM%values(1:myDEM%NROWS,1:myDEM%NCOLS))
            
            IF (DEMName(L-2:L) == 'flt') THEN 
               ! Read the main attribute in FLT
               OPEN(UNIT=10,FILE=TRIM(DEMName),ACTION='READ', &
               FORM='UNFORMATTED',ACCESS='STREAM')
               DO I = 1, myDEM%NROWS
                   READ(10) (myDEM%values(I,J),J=1,myDEM%NCOLS)
               ENDDO
               CLOSE(10)
            ELSEIF (DEMName(L-1:L) == 'nc') THEN
               ! Read the main attribute in NETCDF
               status = NF90_GET_VAR(NC_ID,depth_ID,myDEM%values)
               IF (status /= nf90_noerr) THEN
                   ! We need to read in to dummy and then reshape it back
                   ALLOCATE(DUMMY(myDEM%NCOLS,myDEM%NROWS))
                   CALL Check(NF90_GET_VAR(NC_ID,depth_ID,DUMMY))
                   WRITE(*,'(A)') 'Reshaping depth array'
                   myDEM%values = TRANSPOSE(DUMMY)
               ENDIF
               ! Flipud
               IF (FlipUD.eq.1) THEN
                  DO I = 1,myDEM%NCOLS
                     myDEM%values(:,I) = myDEM%values(myDEM%NROWS:1:-1,I) 
                  ENDDO
               ENDIF
               CALL ChecK(NF90_CLOSE(NC_ID))
            ENDIF
            WRITE(*,'(A)')'done!'
        END SUBROUTINE ReadDEM
    
        subroutine check(status)
          integer, intent ( in) :: status
    
          if(status /= nf90_noerr) then 
             print *, trim(nf90_strerror(status))
             stop "Stopped"
          end if
        end subroutine check  

    END MODULE DEM
    
    PROGRAM DEM2GRD
        
        USE ADCGRID
        USE DEM

        IMPLICIT NONE
        
        CHARACTER(200)                  :: CMD
        CHARACTER(100)                  :: inputFile
        CHARACTER(200)                  :: meshFile
        CHARACTER(200)                  :: flagGrid
        CHARACTER(200)                  :: DEMFile
!        CHARACTER(200)                  :: HDRFile
        CHARACTER(200)                  :: outMeshFile
        CHARACTER(60)                   :: VERSION = 'VERSION 4.NC'

        REAL(8)                         :: mult_fac
        REAL(8)                         :: interpval
        
        INTEGER                         :: I
        INTEGER                         :: cs
        INTEGER                         :: cell_dist
        
        LOGICAL                         :: FoundInputFile
        LOGICAL                         :: exists
        LOGICAL                         :: coord

        TYPE(MESH)                      :: myMesh
        TYPE(MESH)                      :: flagMesh
        TYPE(MESH)                      :: zMesh
        TYPE(DEM_G)                     :: myDEM
        
        IF(IARGC().GE.1)THEN
            I = 0
            DO WHILE(I.LT.IARGC())
                I = I + 1
                CALL GETARG(I,CMD)
                SELECT CASE(TRIM(CMD))
                    CASE('-i','-I')
                        I = I + 1
                        CALL GETARG(I,inputFile)
                        FoundInputFile = .TRUE.
                    CASE('-v','-V','-version')
                        WRITE(*,'(A)')VERSION
                        STOP
                    CASE DEFAULT
                        WRITE(*,'(A)')'ERROR: I did not understand argument: '//TRIM(CMD)
                        STOP
                END SELECT
            ENDDO
        ENDIF
        
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'******************************************************'
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'                 DEM2GRD.F90                          '
        WRITE(*,'(A)')'                 VERSION 4.NC                         '
        WRITE(*,'(A)')' PROGRAM TO INTERPOLATE A FLT/GRD RASTER OR           '
        WRITE(*,'(A)')' NETCDF DEM TO ADCIRC MESH NODES.                     '
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'MATTHEW V. BILSKIE, E.I.                              '
        WRITE(*,'(A)')'Matt.Bilskie@gmail.com'
        WRITE(*,'(A)')'Copyright M. Bilskie (2013)'
        WRITE(*,'(A)')'Please cite:'
        WRITE(*,'(A)')'Bilsie & Hagen (2013) Adv. Water Res., 52, 165-177,'
        WRITE(*,'(A)')'http://dx.doi.org_10.1016_j.advwatres.2012.09.003'
        WRITE(*,'(A)')''
        WRITE(*,'(A)')'GO NEW ZEALAND!'
        WRITE(*,'(A)')'******************************************************'
        WRITE(*,'(A)')''

        IF(.NOT.FoundInputFile)THEN
            WRITE(*,'(A,$)')'Name of input file:'
            READ(*,*)inputFile
            exists = fileExist(TRIM(inputFile),0)
        ENDIF

        OPEN(UNIT=10,FILE=TRIM(inputFile),ACTION='READ')
        READ(10,*)
        READ(10,*)
        READ(10,'(A60)') meshFile
        READ(10,*)flagGrid
        READ(10,*)cs !... 0 for cartesian, 1 for lat/lon
        READ(10,*)DEMFile
!        HDRFile = TRIM(DEMFile)//'.hdr'
!        DEMFile = TRIM(DEMFile)//'.flt'
        READ(10,*)mult_fac
        READ(10,*)cell_dist
        READ(10,*)interpval
        READ(10,'(A60)')outMeshFile
        CLOSE(10)

        IF(cs.EQ.0)THEN
            coord = .FALSE. !...cartesian coordinates
        ELSE
            coord = .TRUE.  !...Geographic (lat/lon)
        ENDIF

        !IF(cell_dist.EQ.0)THEN
        !    WRITE(*,'(A)')'Automatic cell averaging will be used...'
        !ELSEIF(cell_dist.GT.0)THEN
        !    WRITE(*,'(A,I2,A)')'Automatic cell averaging * ',cell_dist, &
        !        ' will be used...'
        !ENDIF

        CALL ReadMesh(meshFile,myMesh)
        exists = FileExist(flagGrid,1)
        IF (exists) THEN
           CALL ReadMesh(flagGrid,flagMesh)
        ELSE
           flagMesh = myMesh
           flagMesh%nodes(:,3) = cell_dist  
        ENDIF
        IF (myMesh%np.NE.flagMesh%np) THEN
            WRITE(*,'(A)')'Mesh and flag grid geometry does not match'
            WRITE(*,*)
            STOP
        ENDIF
        CALL ReadDEM(DEMFile,myDEM)
        CALL CAA(myMesh,flagMesh,zMesh,myDEM,cell_dist,mult_fac,coord,interpval)  
        CALL WriteMesh(zMesh,outMeshFile,coord) 

        WRITE(*,'(A)')''
        
    END PROGRAM

    SUBROUTINE CAA(myMesh,flagMesh,newMesh,myRaster,cdist,mf,coord,interpval)
        USE ADCGRID
        USE DEM
 
        IMPLICIT NONE

        TYPE(MESH),INTENT(IN)           :: myMesh
        TYPE(MESH),INTENT(IN)           :: flagMesh
        TYPE(MESH),INTENT(OUT)          :: newMesh
        TYPE(DEM_G),INTENT(IN)          :: myRaster
        INTEGER,INTENT(IN)              :: cdist
        REAL(8),INTENT(IN)              :: mf !...mult. factor
        LOGICAL,INTENT(IN)              :: coord
        REAL(8),INTENT(IN)              :: interpval

        INTEGER                         :: I
        INTEGER                         :: K
        INTEGER                         :: L
        INTEGER                         :: counter
        INTEGER                         :: no_data_count
        INTEGER                         :: in_pts
        INTEGER                         :: CA
        INTEGER                         :: n_col, n_row

        REAL(8),ALLOCATABLE             :: eSize(:)
        REAL(8)                         :: x_inc
        REAL(8)                         :: y_inc
        REAL(8)                         :: x_botrt
        REAL(8)                         :: y_botrt
        REAL(8)                         :: x_toplt
        REAL(8)                         :: y_toplt
        REAL(8)                         :: N
        REAL(8)                         :: x
        REAL(8)                         :: y
        REAL(8)                         :: z
        REAL(8)                         :: dx
        REAL(8)                         :: dy
        REAL(8)                         :: col
        REAL(8)                         :: row
        REAL(8)                         :: decimal
        REAL(8)                         :: rastz
        REAL(8)                         :: avgz
        REAL(8)                         :: fz !...final z
        REAL(8)                         :: z_check
        REAL(8)                         :: deg2rad
        REAL(8),PARAMETER               :: PI = 3.14159265359D0
        REAL(8),PARAMETER               :: R = 6371.64D0 !km

        INTRINSIC                       :: ABS

        deg2rad = PI / 180.0D0

        CALL ElementSize(myMesh,coord,eSize)

        WRITE(*,'(A)')'Interpolating...'

        newMesh = myMesh

        !IF(coord)THEN
        !    x_inc = R * 1000.D0 * PI * myRaster%x_y_dist / 180.D0
        !    y_inc = x_inc
        !ELSE
        !    x_inc = myRaster%x_y_dist
        !    y_inc = myRaster%x_y_dist
        !ENDIF
        x_inc = myRaster%x_y_dist
        y_inc = myRaster%x_y_dist

        x_botrt = myRaster%x_botlt + myRaster%NCOLS * x_inc
        y_botrt = myRaster%y_botlt
        x_toplt = myRaster%x_botlt
        y_toplt = myRaster%y_botlt + myRaster%NROWS * y_inc

        no_data_count = 0
        in_pts = 0

! WJP: We can easily do this in parallel, use OPENMP
!$omp parallel do private(x,y,z,dx,dy,n_col,n_row,rastz,counter,&
!$omp                     avgz,in_pts,CA,fz,K,L,z_check,no_data_count)       
        DO I = 1, myMesh%np
            x = myMesh%nodes(I,1)
            y = myMesh%nodes(I,2)
            z = myMesh%nodes(I,3)
            
            !IF (I.eq.392657) THEN
            !    continue
            !ENDIF
            !IF(z.NE.0) GOTO 2000 !...Only perform interpolation on values not equal to 0
            IF (interpval.eq.-888.and.z.ne.-888) GOTO 2000 !...Only interpolate values equal to -888

            IF((x.GT.x_toplt).AND.(x.LT.x_botrt).AND.(y.GT.y_botrt).AND. &
                    (y.LT.y_toplt))THEN

            
                dx = x - x_toplt
                dy = y - y_botrt
                
                ! Ensure rounded to nearest value
                n_col = NINT(ABS(dx / x_inc)) + 1
                n_row = NINT(ABS(dy / y_inc)) + 1

                IF(n_col.EQ.myRaster%NCOLS+1) n_col = myRaster%NCOLS
                IF(n_row.EQ.myRaster%NROWS+1) n_row = myRaster%NROWS
		!!!!!!Debugging	file!!!!!!!!!!!!!!!!!!!!!
                !IF (mod(I,1000000).eq.0) THEN
                !   write(6,*) x,y,dx,dy,n_col,n_row     
                !ENDIF
	        !    OPEN(UNIT=17,FILE="debug_output.txt",POSITION="append")
		!    WRITE(17,*) 'myRaster%NROWS = ', myRaster%NROWS
		!    WRITE(17,*) 'myRaster%NCOLS = ', myRaster%NCOLS
		!    WRITE(17,*) 'x_toplt = ', x_toplt
		!    WRITE(17,*) 'x_botrt = ', x_botrt
		!    WRITE(17,*) 'y_toplt = ', y_toplt
		!    WRITE(17,*) 'y_botrt = ', y_botrt
		!    WRITE(17,*) 'dx = ', dx
		!    WRITE(17,*) 'dy = ', dy
		!    WRITE(17,*) 'n_row = ', n_row
		!    WRITE(17,*) 'n_col = ', n_col
		!    CLOSE(17)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                rastz = myRaster%values(n_row,n_col)
                
                ! This part checks the nearest neighbour
                !col = ABS(dx / x_inc)
                !decimal = n_col - col
                !IF(decimal.LT.0.5) n_col = n_col - 1
                !IF(decimal.GT.0.5) n_col = n_col + 1

                !row = ABS(dy / y_inc)
                !decimal = n_row - row
                !IF(decimal.LT.0.5) n_row = n_row - 1
                !IF(decimal.GT.0.5) n_row = n_row + 1

                !IF (n_col.EQ.0) n_col = 1
                !IF (n_row.EQ.0) n_row = 1
                !rastz = myRaster%values(n_row,n_col)

                counter = 0
                avgz = 0.D0
                IF(rastz.NE.MyRaster%nodata_value)THEN !...Node has value
                    in_pts = in_pts + 1
                    ! Correction by WJP:
                    ! Had ceiling on CA at beginning meaning that CA < 1 could
                    ! not be obtained. Changed to calc real at start then
                    ! convert to integer before the averaging approach
                    IF(coord)THEN 
                        N = (0.25D0 * eSize(I)) / &
                        (R * 1000.D0 * PI * myRaster%x_y_dist / 180.D0) 
                    ELSE
                        N = 0.25D0 * eSize(I) / x_inc
                    ENDIF

                    !IF(cdist.GT.0) CA = CA * cdist
                    ! Use flagged mesh
                    ! flagmesh%nodes(I,3) < 1, no smoothing
                    ! flagmesh%nodes(I,3) = 1, direct lookup
                    ! flagmesh%nodes(I,3) > 1, smoothing value

                    IF (ABS(flagMesh%nodes(I,3)).GT.1) THEN
                        N = N * ABS(flagMesh%nodes(I,3))
                    ENDIF
                    !IF (MOD(I,100000).eq.0) THEN
                    !    write(*,*) eSize(I), N, ceiling(N)
                    !ENDIF
                    ! Either do convesion to nearest integer here (meaning 0.5 -
                    ! 1 values become CA = 1 or do after IF statement meaning 
                    ! all values < 1 become direct lookup
                    CA = INT(N)
                    !CA = CEILING(N)
                    IF(CA.LT.1)THEN
                        fz = rastz * mf
                        GOTO 2100
                    ENDIF
!                    CA = NINT(N)
                    
                    DO K = n_row - CA, n_row + CA
                        DO L = n_col - CA, n_col + CA
                            !Make sure K is within the DEM rows & L within the columns
                            IF( (K.GE.CA).AND.(K.LE.myRaster%NROWS).AND.(L.GE.CA).AND.(L.LE.myraster%NCOLS) ) THEN
                            !IF((K.GT.CA).AND.(K.LT.myRaster%NROWS))THEN !...If cell window falls outside the DEM
                                z_check = myRaster%values(K,L)
                                IF(z_check.NE.myRaster%nodata_value)THEN
                                    IF(flagMesh%nodes(I,3).LT.0)THEN
                                        IF(mf*myRaster%values(K,L).GT.0)THEN
                                            avgz = avgz + myRaster%values(K,L)
                                            counter = counter + 1
                                        ENDIF
                                    ELSE
                                        avgz = avgz + myRaster%values(K,L)
                                        counter = counter + 1
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDDO
                    ENDDO

                    IF(counter.NE.0)THEN
                        avgz = avgz / counter
                        fz = avgz * mf
                    ELSE
                        fz = -9999
                    ENDIF
2100                CONTINUE
                    avgz = 0.D0
                    counter = 0
                ELSE
                    no_data_count = no_data_count + 1
                    !fz = 0.D0
                    fz = -9999
                ENDIF
            ELSE !...Node is outside the DEM
                no_data_count = no_data_count + 1
                !fz = 0.D0
                fz = -9999
            ENDIF
            IF (fz == -9999) THEN
               IF (z.eq.-888) THEN
                  ! We wanted to interp here but got null data
                  newMesh%nodes(I,3) = fz
               ENDIF
               ! Else we do not update the value
            ELSE
               IF (interpval.eq.-888.or.fz < interpval) THEN
                  newMesh%nodes(I,3) = fz
               ENDIF
            ENDIF
2000    CONTINUE
        ENDDO
!$omp end parallel do
        WRITE(*,'(A)')'done!'
    END SUBROUTINE CAA
