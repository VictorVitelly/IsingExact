program main
    use iso_fortran_env, only : dp => real64, i32 => int32, i64 => int64
    implicit none

    ! ── Parameters ─────────────────────────────────────────────────────────────
    integer(i32), parameter :: L   = 4       ! lattice side length  (change me)
    integer(i32), parameter :: Nps = 50      ! number of temperature points
    integer(i64), parameter :: Nstates = 2_i64**int(L*L, i64)
    !   L=2 → 16 states  |  L=4 → 65 536  |  L=6 → ~68 billion (slow)

    ! ── Variables ───────────────────────────────────────────────────────────────
    integer(i32) :: spin(L, L)               ! current spin configuration
    integer(i32), allocatable :: ip(:), im(:)! periodic neighbours
    real(dp),     allocatable :: E_res(:)    ! energy per spin at each T

    real(dp)     :: Tf, Ti, T, vol
    real(dp)     :: Ztot, Etot, E_micro, weight
    integer(i64) :: istate
    integer(i32) :: iT, k, row, col

    ! ── I/O ─────────────────────────────────────────────────────────────────────
    open(10, file='energy.dat',    status='replace')
    open(20, file='partition.dat', status='replace')
    open(30, file='entropy.dat', status='replace')

    vol = real(L*L, dp)
    Tf  = 5.1_dp
    Ti  = 0.1_dp

    call init_vecs()

    ! ── Main temperature loop ────────────────────────────────────────────────────
    do iT = 1, Nps
        T    = Tf + (Ti - Tf) * real(iT-1, dp) / real(Nps-1, dp)
        Ztot = 0._dp
        Etot = 0._dp
        write(*,*) iT
        ! Loop over ALL 2^(L²) microstates encoded as integers
        do istate = 0_i64, Nstates - 1_i64

            ! Decode integer bits → spin array
            !   bit k=0 ↔ spin(1,1), bit k=1 ↔ spin(1,2), ..., bit k=L*L-1 ↔ spin(L,L)
            !   bit = 1 → spin = +1,  bit = 0 → spin = -1
            do k = 0, L*L - 1
                row = k / L + 1
                col = mod(k, L) + 1
                spin(row, col) = merge(1, -1, btest(istate, k))
            end do

            E_micro = Hamilt(spin)
            weight  = exp(-E_micro / T)
            Ztot    = Ztot + weight
            Etot    = Etot + weight * E_micro

        end do

        E_res(iT) = Etot / (Ztot * vol)
        write(10,*)  T, E_res(iT)
        write(20,*)  T, Ztot
        write(30,*)  E_res(iT), E_res(iT)/T+log(Ztot)/vol
    end do

    close(10)
    close(20)
    close(30)
contains

    ! ── Build periodic-boundary neighbour index arrays ───────────────────────────
    subroutine init_vecs()
        integer(i32) :: i
        allocate(ip(L), im(L), E_res(Nps))
        do i = 1, L-1
            ip(i) = i + 1
        end do
        ip(L) = 1          ! wraps last site → first site

        do i = 2, L
            im(i) = i - 1
        end do
        im(1) = L          ! wraps first site → last site
    end subroutine init_vecs

    ! ── Ising Hamiltonian  H = -J Σ s_i s_j   (J=1, sum over bonds, no double count)
    function Hamilt(spin) result(H)
        integer(i32), dimension(:,:), intent(in) :: spin
        real(dp)     :: H
        integer(i32) :: i, j

        H = 0._dp
        do i = 1, size(spin, 1)
            do j = 1, size(spin, 2)
                ! Count each bond once: right neighbour + down neighbour
                H = H - real(spin(i,j), dp) &
                      * real(spin(ip(i),j) + spin(i,ip(j)), dp)
            end do
        end do
    end function Hamilt

end program main
