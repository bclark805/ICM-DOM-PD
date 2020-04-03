MODULE MOD_FILENAMES

!INTEGER ifindext
!CHARACTER *(*) fnameprefix
!CHARACTER *(*) fnameext 
  CONTAINS

	!functions:
	!	function ifindext()
	!	function fnameprefix()
	!	function fnameext()
!
!WLong: added a few functions that deals with filenames and extensions etc
!

       INTEGER Function ifindext(FILENAME)

      !
      !Fortran function to find the index of the last dot '.' in a filename string
      !which is useful for detecting whether the filename has an extension
      !
       implicit none
       integer(4)::index_dot
       character(*)::FILENAME
       index_dot=index(TRIM(FILENAME),'.',.TRUE.)
       ifindext=index_dot
       RETURN
       END FUNCTION ifindext

       CHARACTER(1024) Function fnameprefix(FILENAME)
       !
       !Fortran function to find the filename after removing the extension (.dot etc)
       !
       implicit none
       integer(4)::index_dot
       character(*)::FILENAME
       character(LEN=len(FILENAME)):: FILENAME_PREFIX,FILENAME_NEW

       FILENAME_NEW=TRIM(FILENAME)
       index_dot=index(TRIM(FILENAME),'.',.TRUE.)
       IF(index_dot > 0) THEN
         FILENAME_PREFIX=TRIM(FILENAME_NEW(1:index_dot-1))
       ELSE
         FILENAME_PREFIX=TRIM(FILENAME_NEW)
       ENDIF
       fnameprefix=FILENAME_PREFIX
       RETURN
       END FUNCTION fnameprefix

       CHARACTER(1024) Function fnameext(FILENAME)
       !
       !fortran function to find extension name if any
       !
       implicit none
       integer(4)::index_dot,length
       character(*)::FILENAME
       character(LEN=len(FILENAME))::FILENAME_EXT,FILENAME_NEW

       FILENAME_NEW=TRIM(FILENAME)
       index_dot=index(TRIM(FILENAME),'.',.TRUE.)
       length=len(TRIM(FILENAME_NEW))

       IF(index_dot > 0) THEN
          FILENAME_EXT = TRIM(FILENAME_NEW(index_dot:length))
       ELSE
          FILENAME_EXT=''
       ENDIF
       fnameext=FILENAME_EXT
       RETURN
       END FUNCTION fnameext


END MODULE MOD_FILENAMES

