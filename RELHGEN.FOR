      program relhgen
C***********************************************************************
C     GENERATES A RELATIVE HUMIDITY FILE FOR THE NEEDS OF THE 
C     BBRI PAS-CALIBRATION WALL VALIDATION EXERCISE
C     
C     ATHANASSIOS A. ARGIRIOU, 11-4-1994
C***********************************************************************
      open(100,file='relhum')
      rh=40
      do i=1,24
      do j=1,21
      write(100,'(f9.0)')rh
      enddo
      enddo
      close(100)
      end
       