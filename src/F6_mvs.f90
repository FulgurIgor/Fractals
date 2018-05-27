MODULE Pixel

CONTAINS
	FUNCTION ITER_Z(zC, CConst, inpA, inpB, colorR, colorG, colorB, iterations)
		COMPLEX(8), INTENT(IN) :: zC, CConst
		COMPLEX(8) :: c, z
		REAL(8), INTENT(IN) :: inpA, inpB, colorR, colorG, colorB
		REAL(8) :: iter2
		INTEGER :: iter, iterations
		CHARACTER, DIMENSION(0:2) :: ITER_Z

		iter = 0
		z = zC
		c = inpB*z
		do while ((iter .le. iterations) .and. (CDABS(z) .le. 2._8))
			z = inpA*z**2 + c + CConst
			iter = iter + 1
		end do
		iter2 = DFLOAT(iter) - DLOG(DLOG(CDABS(z))) / DLOG(2._8)
		ITER_Z(0) = CHAR(myround(iter2 * colorR))
	        ITER_Z(1) = CHAR(myround(iter2 * colorG))
	        ITER_Z(2) = CHAR(myround(iter2 * colorB))
	        RETURN
	        
	END FUNCTION ITER_Z

	INTEGER function myround(a) result(ret)
		implicit none
		REAL(8), INTENT(IN) :: a
		ret = ABS(INT(a -  255._8 * NINT(a / 255._8)))
	end function myround

END MODULE Pixel

PROGRAM frac
	USE Pixel
	INCLUDE 'mpif.h'
	CHARACTER(len=128) :: ARG128
	CHARACTER(len=32)  :: ARG
	REAL(8) :: xshift, yshift, x0, y0, mag_ratio, mag_ratio_pow_start, mag_ratio_pow_finish
	REAL(8) :: CXmin, CXmax, CXshag, CYmin, CYmax, CYshag, inpA, inpB, colorR, colorG, colorB, iter2
	COMPLEX(8) :: CConst, z, c
	INTEGER :: resox, resoy, iterations, x, y, iter, ompthread, sar
	INTEGER :: i, rank, size, ierr                 !MPI
	INTEGER :: status(MPI_STATUS_SIZE), request    !MPI

	CHARACTER, DIMENSION(:), ALLOCATABLE :: matrix, mpiarray
	
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
	
	IF (rank .eq. 0) THEN
		IF (IARGC() .lt. 18) THEN
			WRITE (*, *) "INPUT ALL PARAMETRES:"
			WRITE (*, *) " 1: INTEGER - Размер изображения по X"
			WRITE (*, *) " 2: INTEGER - Размер изображения по Y"
			WRITE (*, *) " 3: REAL - Смещение изображения по X (-1)"
			WRITE (*, *) " 4: REAL - Смещение изображения по Y (-1)"
			WRITE (*, *) " 5: REAL - Смещение комплексной оси по X"
			WRITE (*, *) " 6: REAL - Смещение комплексной оси по Y"
			WRITE (*, *) "Приближение фрактала, A^B - начальное значение" 
			WRITE (*, *) "A^C - конечное значение"
			WRITE (*, *) " 7: REAL - A"
			WRITE (*, *) " 8: REAL - B"
			WRITE (*, *) " 9: REAL - C"
			
			WRITE (*, *) "10: INTEGER - число итераций"
			WRITE (*, *) "11: INTEGER - OMP-потоки (в MPI версии не используется)"
			WRITE (*, *) "Фрактал A*z**2 + B*c + U"
			WRITE (*, *) "U - комплексное число"
			WRITE (*, *) "A, B - действительные числа"
			WRITE (*, *) "12: REAL - Действительна часть параметра U"
			WRITE (*, *) "13: REAL - Мнимая часть параметра U"
			WRITE (*, *) "14: REAL - A"
			WRITE (*, *) "15: REAL - B"
			WRITE (*, *) "Цветовая гамма"
			WRITE (*, *) "16: REAL - Red"
			WRITE (*, *) "17: REAL - Green"
			WRITE (*, *) "18: REAL - Blue"
			WRITE (*, *) "19: CHARs - название выходного файла"
			CALL EXIT(1)
		END IF
	END IF
	
	!читаем разрешение изображения
	CALL GETARG(1, ARG)
	READ (ARG, *) resox
	CALL GETARG(2, ARG)
	READ (ARG, *) resoy
	
	!отражение и растяжение
	CALL GETARG(3, ARG)
	READ (ARG, *) xshift
	CALL GETARG(4, ARG)
	READ (ARG, *) yshift
	
	!координаты центра в комплексной системе координат
	CALL GETARG(5, ARG128)
	READ (ARG128, *) x0
	CALL GETARG(6, ARG128)
	READ (ARG128, *) y0
	
	!приближение
	CALL GETARG(7, ARG)
	READ (ARG, *) mag_ratio
	CALL GETARG(8, ARG)
	READ (ARG, *) mag_ratio_pow_start
	
	CALL GETARG(9, ARG)
	READ (ARG, *) mag_ratio_pow_finish
	
	!число итераций
	CALL GETARG(10, ARG)
	READ (ARG, *) iterations
	
	!число OMP потоков
	CALL GETARG(11, ARG)
	READ (ARG, *) ompthread
	
	!внешнее комплексное число
	CALL GETARG(12, ARG128)
	READ (ARG128, *) inpA 
	CALL GETARG(13, ARG128)
	READ (ARG128, *) inpB
	 CConst = CMPLX(inpA, inpB)
	
	!коэффициенты при комплексных числах
	CALL GETARG(14, ARG128)
	READ (ARG128, *) inpA 
	CALL GETARG(15, ARG128)
	READ (ARG128, *) inpB
	
	!цветовая гамма
	CALL GETARG(16, ARG)
	READ (ARG, *) colorR
	CALL GETARG(17, ARG)
	READ (ARG, *) colorG
	CALL GETARG(18, ARG)
	READ (ARG, *) colorB
	
	!имя выходного файла
	CALL GETARG(19, ARG128)
	
	!начало цикла
	
	!пре-расчет (настройки изображения и координатной сетки)
	mag_ratio = mag_ratio**mag_ratio_pow_start
	xshift = (xshift / mag_ratio) + x0
	yshift = (yshift / mag_ratio) + y0
	
	CXmin = -1._8 / mag_ratio
	CXmax = 1._8 / mag_ratio
	CXshag = (DABS(CXmin) + DABS(CXmax)) / resox
	CYmin = -1._8 / mag_ratio
	CYmax = 1._8 / mag_ratio
	CYshag = (DABS(CYmin) + DABS(CYmax)) / resoy
	sar = resox*resoy*3/size
	
	ALLOCATE(matrix(0 : resox * resoy * 3))
	ALLOCATE(mpiarray(0 : sar-1))
	
	do y = resoy/size*rank, resoy/size*(rank+1)
		! !$OMP PARALLEL DO PRIVATE(z)
		do x = 0, resox - 1
			z = DCMPLX(x * CXshag + xshift, y * CYshag + yshift)
			matrix((x + y * resox) * 3 + 0 : (x + y * resox) * 3 + 2) = ITER_Z(z, CConst, inpA, inpB, colorR, colorG, colorB, iterations)
		end do
		! !$OMP END PARALLEL DO
	end do
	
	CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
	
	IF (rank .eq. 0) THEN
		DO i = 1, size-1		
			CALL MPI_IRECV(mpiarray, sar, MPI_CHAR, i, 12, MPI_COMM_WORLD, request, ierr)
			CALL MPI_WAIT(request, status, ierr)
			matrix(sar*i : sar*(i+1)-1) = mpiarray
		END DO
		open(unit=7, status = 'replace', file = TRIM(TRIM(ARG128)//".pnm"))
	    	write(7, '(a)') "P6" !(a) это текстовая строка
		write(7, '(a)') ""
		write(7, '(I0, a, I0)') resox, ' ', resoy
		write(7, '(I0)') 255
		write(7, *) matrix
		close(7)
	ELSE
		mpiarray = matrix(sar*rank : sar*(rank+1)-1)
		CALL MPI_ISEND(mpiarray, sar, MPI_CHAR, 0, 12, MPI_COMM_WORLD, request, ierr)	
		CALL MPI_WAIT(request, status, ierr)
	END IF
	
	CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
	DEALLOCATE(matrix)
	DEALLOCATE(mpiarray)
	
	!конец цикла
	
	CALL MPI_FINALIZE(ierr)
END