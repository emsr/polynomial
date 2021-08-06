! gfortran -g -o madsen_reid madsen_reid.f95

double precision function norml1(z)
    double complex z
    norml1 = abs(real(z)) + abs(aimag(z))
end function

! Evaluate polynomial at z, set fz, return squared modulus.
double precision function norml2(z)
    double complex z
    norml2 = real(z)**2 + aimag(z)**2
end function

! Evaluate polynomial at z, set fz, return squared modulus.
double precision function eval(z, fz, n1, a)

    implicit none
    double complex z, fz, a(n1), p
    double precision norml2
    integer n1, n, i

    n = n1 - 1
    p = a(1)
    do i = 1, n
        p = p * z + a(i + 1)
    end do
    fz = p

    eval = real(fz)**2 + aimag(fz)**2
end function

!
! Root search...
!
subroutine solve(a1, m, root, a, mp1)

    implicit none

    double complex a1(m+1), root(m), a(mp1)
    double complex z0, f0z, z, dz, f1z, fz, w, fw, dzk
    double precision f0, ff, f, fa, dz1, fmin, f2
    logical stage1, div2
    double precision BIG, SMALL, BASE, EPS, SSMALL, ALOGB
    double precision norml1, norml2, eval

    integer m, mp1, n, j, i, n1, k
    double precision r0, t, u, u0, u1, u2, r, r1

    data BIG/1.0d+70/ ! Overflow limit
    data SMALL/1.0d-70/ ! Underflow limit.
    !!!data BASE/16/
    data BASE/2/
    data EPS/2.3d-16/

    SSMALL = sqrt(SMALL)
    ALOGB = log(BASE)

    n = m

    ! Store original polynomial in a and in root.
    j = m + 1
    a(j) = a1(1)
    do i = 1, m
        root(i) = a1(i)
        a(i) = a1(j)
        j = j - 1
    end do

    ! Test for zeros at infinity.
 20 continue
    if (norml1(a(1)) .gt. 0.0) then
        goto 40
    end if
    do i = 1, n
        a(i) = a(i + 1)
    end do
    root(n) = BIG
    n = n - 1
    if (n .gt. 0) then
        goto 20
    end if
    goto 310

 40 continue
    if (n .le. 1) then
        goto 260
    end if

    n1 = n + 1

    ! Scale the coefficients.
    u1 = 0.0
    u2 = BIG
    do k = 1, n1
        u = norml1(a(k))
        if (u .le. 0.0) then
            cycle
        end if
        if (u .gt. u1) then
            u1 = u
        end if
        if (u .lt. u2) then
            u2 = u
        end if
    end do
    u = sqrt(u1) * sqrt(u2)
    i = -log(u) / ALOGB ! ilogb?
    u = BASE**i
    do k = 1, n
        a(k) = u * a(k)
        a1(k) = a(k) * dble(n1 - k)
    end do
    a(n1) = u * a(n1)

    ! Test for zeros at (0, 0)
    z = complex(0.0, 0.0)
    if (norml1(a(n1)) .le. SSMALL) then
        goto 290
    end if
    z0 = complex(0.0, 0.0)
    f0 = norml2(a(n1))
    fmin = f0 * (dble(n) * 16.0 * EPS)**2

    ! z is the current point
    ! f = |f(z)|^2
    ! z0 is the last point
    ! f0z = f'(z0)
    ! f0 = |f(z0)|^2
    ! r0 = 3|z - z0|
    ! dz is the last tentative step if the last step was successful or is the required next step.
    ff = f0
    u0 = f0
    t = BIG
    do k = 1, n
        u = norml2(a(k))
        if (u .eq. 0.0) then
            cycle
        end if
        u = log(u0 / u) / dble(2 *(n1 - k))
        if (u .lt. t) then
            t = u
        end if
    end do
    t = exp(t)
    f0z = a(n)
    z = complex(1.0, 0.0)
    if (norml1(f0z) .gt. 0.0) then
        z = -a(n1) / a(n)
    end if
    u = 0.5 * t / norml1(z)
    z = u * z
    dz = z
    f =  eval(z, fz, n + 1, a)
    r0 = 0.5 * t

    ! Calculate tentative step.
120 continue
    u = eval(z, f1z, n, a1)
    if (u .eq. 0.0) then
        goto 140
    end if
    dz = -fz / f1z
    f2 = norml2(f0z - f1z) / norml2(z0 - z)
    stage1 = (f * f2 / u .gt. 0.25 * u) .or. (f .ne. ff)
    r = norml1(dz)
    if (r .le. 3.0 * r0) then
        goto 150
    end if

    dz1 = real(dz)
    dz = complex(1.8 * dz1 - 2.4 * aimag(dz) * r0 / r, &
                 2.4 * dz1 + 1.8 * aimag(dz) * r0 / r)
    goto 150

140 continue
    dz1 = real(dz)
    dz = complex(1.8 * dz1 - 2.4 * aimag(dz), &
                 2.4 * dz1 + 1.8 * aimag(dz))
    stage1 = .true.

150 continue
    f0z = f1z

    ! Find the next point in te iteration.
    ! This is where iteration starts if the previous one was unsuccessful.
160 continue
    z0 = z
    f0 = f
    dzk = dz
    z = z0 + dz
    ! If either part of z is small replace by zero to avoid underflows.
    if (abs(real(z)) .lt. EPS * abs(aimag(z))) then
        z = complex(0.0, aimag(z))
    end if
    if (abs(aimag(z)) .lt. EPS * abs(real(z))) then
        z = complex(real(z), 0.0)
    end if
    w = z
    f = eval(z, fz, n + 1, a)
    ff = f
    if (.not.stage1) then
        goto 240
    end if

    ! Beginning of stage 1 search.
    j = 1
    div2 = f .ge. f0
180 continue
    if (div2) then
        dz = 0.5 * dz
        w = z0 + dz
    else
        w = w + dz
    end if

    fa = eval(w, fw, n + 1, a)
    if (fa .ge. f) then
        goto 240
    end if
    f = fa
    fz = fw
    z = w
    j = j + 1
    if (div2 .and. j .eq. 3) then
        goto 220
    end if
    if (j .le. n) then
        goto 180
    end if
    goto 240
220 continue
    dz1 = real(dz)
    dz = complex(0.6 * dz1 - 0.8 * aimag(dz), &
                 0.8 * dz1 + 0.6 * aimag(dz))
    z = z0 + dz
    f = eval(z, fz, n + 1, a)

    ! End of stage 1 search.

240 continue
    r0 = norml1(z0 - z)

    ! Convergence test.
    if (f .lt. f0) then
        goto 250
    end if
    z = z0
250 continue
    r1 = norml1(z)
    if (r0 .lt. EPS * r1) then
        goto 270
    end if
    if (f .lt. f0) then
        goto 120
    end if
    f = f0
    if (f .lt. fmin) then
        goto 270
    end if
    dz = complex(-0.3 * real(dzk) + 0.4 * aimag(dzk), &
                 -0.4 * real(dzk) - 0.3 * aimag(dzk))
    stage1 = .true.
    goto 160

260 continue
    ! Deal with n .eq. 1 case.
    z = -a(2) / a(1)
    goto 290

    ! Deflate, store root, restore coefficient of original polynomial.
270 continue
    do k = 2, n
        a(k) = a(k - 1) * z + a(k)
    end do
290 continue
    a1(n) = root(n)
    root(n) = z
    n = n - 1
    if (n .lt. 1) then
        return
    else if (n .eq. 1) then
        goto 260
    else
        goto 40
    end if
310 continue
    return
end subroutine



!
! Main
!
program mr

    implicit none

    integer m, k
    double complex a1(21), a(21), root(20)

    write(6,*) 'Enter degree: '
    read(5,*) m
    do k = 1, m + 1
        write(6,*) 'Enter coefficient ', k - 1, ': '
        read(5,*) a1(m+2-k)
    end do
    
    call solve(a1, m, root, a, m+1)

    do k = 1, m
        write(6,*) 'Zero ', k-1, ': ', root(k)
    end do

end program
