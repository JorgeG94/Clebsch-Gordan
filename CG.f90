module math

IMPLICIT NONE

PUBLIC FACT, getdum, getl

REAL(KIND=8), PARAMETER :: ONE = 1.0D+00

CONTAINS

! SUBROUTINES
! DECK FACTORIAL


	SUBROUTINE FACT(N, FAC)

! Calculates factorials up to N=100

	IMPLICIT NONE
	INTEGER N
	INTEGER I
	REAL(KIND=8), DIMENSION(:),  ALLOCATABLE :: FAC(:)

	IF (N .gt. 100) then
		STOP
	END IF

	FAC(0) = ONE

	DO I=1, N
		FAC(I) = FAC(I-1)*I
	END DO

	RETURN

	END SUBROUTINE FACT


! DECK CLEBSCH GORDAN COEFFICIENTS
! Required subroutines for the Clebsch Gordan coefficients

	SUBROUTINE getdum(l1, l2, DUM, lmin, lmax)

! Triangle inequality to determine the array of ls

	IMPLICIT NONE
	INTEGER :: l1, l2, lmin, lmax, DUM

	lmin = ABS(l1 - l2)
	lmax = l1 + l2
	DUM = lmax - lmin + 1

	RETURN

	END SUBROUTINE getdum



	SUBROUTINE getl(lmin, lmax, l)

! Gets the array of ls

	IMPLICIT NONE
	INTEGER :: lmin, lmax
	INTEGER :: I
	INTEGER, DIMENSION(:), ALLOCATABLE :: l(:)

	DO I=lmin, lmax, 1
		l(I) = I
	END DO


	RETURN


	END SUBROUTINE getl


	SUBROUTINE special_cases(l1,l2,l,DUM,FAC,lmin,lmax)

		IMPLICIT NONE

! Checks for the special cases of the CG coefficients
! C(l1,l2,l;000) = 0 unless l1 + l2 + l = even (case 1)
! C(l1,l2,l;000) where l = l1 + l2 (case 2)
! case 1: equation A-161 Theory of Molecular Fluids, Gray and Gubbins
! case 2: equation A-162 Theory of Molecular Fluids, Gray and Gubbins

	INTEGER, INTENT(IN) :: l1, l2, DUM, lmin, lmax
	INTEGER :: I, l(DUM), TEST
	INTEGER :: L12
	REAL(KIND=8) :: FAC(0:100)

	L12 = l1 + l2


	DO I=lmin, lmax
		l(I) = I
		TEST = L12 + l(I)
		IF (MOD(TEST,2) .eq. 0 .AND. TEST .eq. 2*L12) then
				CALL case_l1l2(l1,l2,DUM,FAC, I) 
		ELSE IF (MOD(TEST,2) .eq. 0) then
			CALL case_even(l1,l2,DUM,I,FAC)
		ELSE
			CALL case_odd(DUM, I)
		END IF
	END DO

	RETURN


	END SUBROUTINE special_cases

	subroutine case_even(l1,l2,DUM,I,FAC)

		IMPLICIT NONE

	! Calculates the l1 + l2 = even case following the equation
	! found in Gray and Gubbins

	INTEGER, INTENT(IN) :: l1, l2, DUM, I
	INTEGER :: x
	REAL(KIND=8) :: FAC(0:100), taocalc, num, denom
	REAL(kind=8) :: c1, c2, c3, c4,cexp, csqrt
	REAL(KIND=8) :: CG(DUM)

		x = l1 + l2 + I
		write(*,*) x
		
		num = FAC(x/2)

		csqrt = sqrt(REAL((2*I + 1))/REAL((x + 1)))
		
		denom = SQRT(FAC(x))
		
		taocalc = num/denom
		c1 = taocalc

		x = l1 + l2 - I
		write(*,*) x

		num = FAC(x/2)

		cexp = (-1.0)**(x/2)

		denom = SQRT(FAC(x))

		taocalc = num/denom
		c2 = taocalc	

		x= l1 - l2 + I
		write(*,*) x
		num = FAC(x/2)

		denom = SQRT(FAC(x))

		taocalc = num/denom
		c3 = taocalc


		x = -l1 + l2 + I
		write(*,*) x
		num = FAC(x/2)

		denom = SQRT(FAC(x))

		taocalc = num/denom
		c4 = taocalc


		CG(I) = cexp*csqrt*c1/(c2*c3*c4)

		write(*,*) "EVEN", l1,l2,I, CG(I)


		return

	end subroutine case_even

	subroutine case_odd(DUM,I)

		IMPLICIT NONE


		INTEGER, INTENT(IN) ::  DUM, I
		REAL(KIND=8) :: CG(DUM)

		CG(I) = 0.0

		write(*,*) "ODD",I, CG(I)

		return 

	end subroutine case_odd

	subroutine case_l1l2(l1,l2,DUM,FAC, I)

		IMPLICIT NONE

	!Calculates the CG coefficients for the case l1 + l2 = l 

		INTEGER, INTENT(IN) :: l1, l2, DUM, I
		REAL(KIND=8) :: Fl1, Fl2,  F2l1, F2l2, Fl(DUM), F2l(DUM)
		REAL(KIND=8) :: FAC(0:100),  CG(DUM)

		Fl1 = FAC(l1)
		Fl2 = FAC(l2)
		F2l1 = FAC(2*l1)
		F2l2 = FAC(2*l2)

		Fl(I) = FAC(I)
		F2l(I) = FAC(2*I)
		CG(I) = (Fl(I)/(Fl1*Fl2))*sqrt((F2l1*F2l2)/F2l(I))

		write(*,*) "l1 + l2 = l",l1,l2,I, CG(I)


		return 

	end subroutine case_l1l2



	subroutine case_fullCG(l1,l2,DUM,FAC,lmin,I,m1,m2,mm,lmax)

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: l1, l2, DUM, lmin, m1, m2, mm, lmax
		REAL(KIND=8) :: FAC(0:100), CG(DUM)
		REAL(kind=8), dimension(:), allocatable :: sum(:)
		INTEGER :: I, z, nz
		REAL(KIND=8) :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11, csqrt1
		REAL(KIND=8) :: d1,d2,d3,d4,d5,d6,d7, denom, full
		INTEGER :: t1,t2,t3,t4,t5,t6

		nz = 0


		allocate(sum(0:4))


			!write(*,*) l1,l2,I, m1, m2, mm

		! To do: check that the factorials are non negative


		 
		do I=lmin, lmax
			 
		c1 = 2*I + 1
		c2 = FAC(l1+l2-I)
		c3 = FAC(l1-l2+I)
		c4 = FAC(-l1 + l2 + I)
		c5 = FAC(l1 + l2 + I + 1)

		t1 = l1 + m1
		t2 = l1 - m1
		t3 = l2 + m2
		t4 = l2 - m2
		t5 = I + mm
		t6 = I - mm

		if (t1 .lt. 0 .or. t2 .lt. 0 .or. & 
			t3 .lt. 0 .or. t4 .lt. 0 .or. &
			t5 .lt. 0 .or. t6 .lt. 0 ) then 
				goto 100 !non physical factorials
		else 
				continue !physical factorials
		end if

		c6 = FAC(l1 + m1)
		c7 = FAC(l1 - m1)
		c8 = FAC(l2 + m2)
		c9 = FAC(l2 - m2)
		c10 = FAC(I + mm)
		c11 = FAC(I - mm)


		csqrt1 = sqrt((c1*c2*c3*c4*c6*c7*c8*c9*c10*c11)/c5)

			!sum(z-1) = 0.0
			sum = 0.0
			do z=0, 4 !loops over z until the factorials are all physical
				
				d1 = l1 + l2 - I - z
				d2 = l1 - m1 - z
				d3 = l2 + m2 - z
				d4 = I - l2 + m1 + z
				d5 = I - l1 - m2 + z

					if (d1 .lt. 0 .or. d2 .lt. 0 .or. & 
					d3 .lt. 0 .or. d4 .lt. 0 .or. &
					d5 .lt. 0 ) then 
						goto 200 !non physical factorial
					else 
						continue !physical factorial 
					end if

				d1 = FAC(l1 + l2 - I - z)
				d2 = FAC(l1 - m1 - z)
				d3 = FAC(l2 + m2 - z)
				d4 = FAC(I - l2 + m1 + z)
				d5 = FAC(I - l1 - m2 + z)
				d6 = FAC(z)
				d7 = (-1)**z
				
				denom = d1*d2*d3*d4*d5*d6
				full = d7/denom

				sum(z) = sum(z-1) + full
				nz = z
				




				200 continue !non physical factorials, loop again

			end do 

			CG(I) = csqrt1*sum(nz)

			write(*,*) "**************************************"
			write(*,*) l1, l2, I, m1, m2, mm, CG(I)

			sum = 0.0
			
			100 continue
		end do 
	
	end subroutine case_fullCG

! FUNCTIONS

	FUNCTION makem(m1,m2) RESULT(m)

		IMPLICIT NONE

	! Gets the m value of the CG coefficient

	INTEGER, INTENT(IN) :: m1, m2
	!INTEGER, DIMENSION(1) :: m
	INTEGER :: m
	INTEGER :: mcalc

	mcalc = m1 + m2

	m = mcalc

	END FUNCTION makem



end module math

program CG

use math

implicit none

  INTEGER :: N, I
	INTEGER :: l1, l2, lmin, lmax, DUM, m1, m2, mm
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: FAC(:)
	INTEGER, DIMENSION(:), ALLOCATABLE :: l(:)
    ALLOCATE(FAC(0:100))

    N=50
	l1 = 4
	l2 = 4
	m1 = 1
	m2 = 1

	CALL getdum(l1, l2, DUM, lmin, lmax)

	ALLOCATE(l(DUM))

	CALL getl(lmin, lmax, l)

	CALL FACT(N, FAC)

	!write(*,*) FAC(2)

	mm = makem(m1,m2)

	if ( m1 .eq. 0 .AND. m2 .eq. 0) then

		CALL special_cases(l1,l2,l,DUM,FAC,lmin,lmax)

	else
		!write(*,*) "m1 and m2 different from 0, proceeding to full calculation"

		call case_fullCG(l1,l2,DUM,FAC,lmin,I,m1,m2,mm,lmax)

	end if


end program CG
