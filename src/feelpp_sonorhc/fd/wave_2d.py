"""This module implements the solution of the acoustic wave equation in 2D, 
in a homogeneous rectangular medium with Neumann boundary conditions, 
using the finite difference method.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import dok_array
from typing import Callable


def get_fdtd_matrix(c: float, dx: float, dy: float, nx: int, ny: int) -> sp.csr_array:
    """Get the FDTD matrix for a 2D wave equation with a constant velocity.

    Args:
        c: Velocity of the wave.
        dx: Spatial step in x direction.
        dy: Spatial step in y direction.
        nx: Number of spatial steps in x direction.
        ny: Number of spatial steps in y direction.

    Returns:
        FDTD matrix.
    """

    # matrix size
    N = (nx + 1) * (ny + 1)
    # initialize sparse matrix
    A = dok_array((N, N),dtype=np.float32)

    # in the interior of the domain
    for i in range(1, nx): # i = 1 à (nx - 1)
        for j in range(1, ny): # j = 1 à (ny - 1)
            # bijection between (i,j) et n
            n = i + j * (nx + 1)
            # coefficients computation
            A[n, n - nx - 1] = - c**2 / dy**2
            A[n, n - 1] = - c**2 / dx**2
            A[n, n] = 2 * c**2 * ((1 / dx**2) + (1 / dy**2))
            A[n, n + 1] = - c**2 / dx**2
            A[n, n + nx + 1] = - c**2 / dy**2

    # at bottom boundary (j = 0)
    for n in range(1, nx): # n = 1 à (nx - 1)
        # coefficients computation
        A[n, n - 1] = - c**2 / dx**2
        A[n, n] = 2 * c**2 * ((1 / dx**2) + (1 / dy**2))
        A[n, n + 1] = - c**2 / dx**2
        A[n, n + nx + 1] = - 2 * c**2 / dy**2

    # at left boundary (i = 0)
    for j in range(1, ny): # j = 1 à (ny - 1)
        # bijection between (i,j) et n
        n = j * (nx + 1)
        # coefficients computation
        A[n, n - nx - 1] = - c**2 / dy**2
        A[n, n] = 2 * c**2 * ((1 / dx**2) + (1 / dy**2))
        A[n, n + 1] = - 2 * c**2 / dx**2
        A[n, n + nx + 1] = - c**2 / dy**2

    # at bottom left corner (i = 0, j = 0)
    n = 0
    # coefficients computation
    A[n, n + 1] = - np.sqrt(2) / (2 * dx)
    A[n, n] = (np.sqrt(2) / 2) * ((1 / dx) + (1 / dy))
    A[n, n + nx + 1] = - np.sqrt(2) / (2 * dy)

    # at left top corner (i = 0, j = ny)
    n = ny * (nx + 1)
    # coefficients computation
    A[n, n + 1] = - np.sqrt(2) / (2 * dx)
    A[n, n] = (np.sqrt(2) / 2) * ((1 / dx) + (1 / dy))
    A[n, n - nx - 1] = - np.sqrt(2) / (2 * dy)

    # at top boundary (j = ny)
    for i in range(1, nx): # i = 1 à (nx - 1)
        # bijection between (i,j) et n
        n = i + ny * (nx + 1)
        # coefficients computation
        A[n, n - nx - 1] = - 2 * c**2 / dy**2
        A[n, n - 1] = - c**2 / dx**2
        A[n, n] = 2 * c**2 * ((1 / dx**2) + (1 / dy**2))
        A[n, n + 1] = - c**2 / dx**2

    # at right top corner (i = nx, j = ny)
    n = nx + ny * (nx + 1)
    # coefficients computation
    A[n, n] = (np.sqrt(2) / 2) * ((1 / dx) + (1 / dy))
    A[n, n - 1] = - np.sqrt(2) / (2 * dx)
    A[n, n - nx - 1] = - np.sqrt(2) / (2 * dy)

    # at right boundary (i = nx)
    for j in range(1, ny): # j = 1 à (ny - 1)
        # bijection between (i,j) et n
        n = nx + j * (nx + 1)
        # coefficients computation
        A[n, n - nx - 1] = - c**2 / dy**2
        A[n, n - 1] = - 2 * c**2 / dx**2
        A[n, n] = 2 * c**2 * ((1 / dx**2) + (1 / dy**2))
        A[n, n + nx + 1] = - c**2 / dy**2

    # at right bottom corner (i = nx, j = 0)
    n = nx
    # coefficients computation
    A[n, n - 1] = - np.sqrt(2) / (2 * dx)
    A[n, n] = (np.sqrt(2) / 2) * ((1 / dx) + (1 / dy))
    A[n, n + nx + 1] = - np.sqrt(2) / (2 * dy)

    return A.tocsr()


def get_fdtd_second_member(c: float, dx: float, dy: float, nx: int, ny: int, S: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Get the FDTD source vector for a 2D wave equation with a constant velocity.
    
    Args:
        c: Velocity of the wave.
        dx: Spatial step in x direction.
        dy: Spatial step in y direction.
        nx: Number of spatial steps in x direction.
        ny: Number of spatial steps in y direction.
        S: Equation's source term vector.
        K: Neumann second member vector.
        
    Returns:
        FDTD second member.
    """

    # initialize second member
    F = np.zeros((nx + 1) * (ny + 1))

    assert S.shape == K.shape == F.shape, "S, K and F must have the same shape"

    # in the interior of the domain
    for i in range(1, nx): # i = 1 à (nx - 1)
        for j in range(1, ny): # j = 1 à (ny - 1)
            # bijection between (i,j) et n
            n = i + j * (nx + 1)
            # coefficients computation
            F[n] = S[n]

    # at bottom boundary (j = 0)
    for n in range(1, nx): # n = 1 à (nx - 1)
        # coefficients computation
        F[n] = S[n] + (2 * c**2 / dy) * K[n]

    # at left boundary (i = 0)
    for j in range(1, ny): # j = 1 à (ny - 1)
        # bijection between (i,j) et n
        n = j * (nx + 1)
        # coefficients computation
        F[n] = S[n] + (2 * c**2 / dx) * K[n]

    # at bottom left corner (i = 0, j = 0)
    F[0] = K[0]

    # at left top corner (i = 0, j = ny)
    n = ny * (nx + 1)
    F[n] = K[n]

    # at top boundary (j = ny)
    for i in range(1, nx): # i = 1 à (nx - 1)
        # bijection between (i, j) et n
        n = i + ny * (nx + 1)
        # coefficients computation
        F[n] = S[n] + (2 * c**2 / dy) * K[n]

    # at right top corner (i = nx, j = ny)
    n = nx + ny * (nx + 1)
    F[n] = K[n]

    # at right boundary (i = nx)
    for j in range(1, ny): # j = 1 à (ny - 1)
        # bijection between (i, j) et n
        n = nx + j * (nx + 1)
        # coefficients computation
        F[n] = S[n] + (2 * c**2 / dx) * K[n]

    # at right bottom corner (i = nx, j = 0)
    F[nx] = K[nx]

    return F


def solve_wave2D_fdtd(c: float, s: Callable, f: Callable, g: Callable, k: Callable, L: float, nx: float, H: float, ny: float, T: float) -> np.ndarray:
    """Get the solution of the wave equation with a constant velocity and a source term.
    
    Args:
        c: Velocity of the wave.
        s: Source term function.
        f: Fist initial condition function.
        g: Second intial condition function.
        k: Neumann condition function.
        L: Length of the domain.
        nx: Number of spatial steps in x direction.
        H: Height of the domain.
        ny: Number of spatial steps in y direction.
        T: Final time.
        
    Returns:
        Solution of the wave equation per time step. Each column is a time step.
    """

    # spatial step
    dx = L / nx
    dy = H / ny

    # discretization
    x = np.linspace(0, L, nx + 1)
    y = np.linspace(0, H, ny + 1)
    x, y = np.meshgrid(x, y)
    x = x.flatten()
    y = y.flatten()

    # time step (with CFL condition)
    dt = 1. / (c * np.sqrt((1 / dx**2) + (1 / dy**2)))

    # initial conditions
    U0 = f(x, y)
    V0 = g(x, y)
    assert U0.shape == V0.shape, "U0 and V0 must have the same shape"

    # FDTD matrix
    A = get_fdtd_matrix(c, dx, dy, nx, ny)

    # FDTD second member (t=0)
    F = get_fdtd_second_member(c, dx, dy, nx, ny, s(x, y, 0), k(x, y, 0))

    # solve
    ncol = int(T / dt) + 1

    U = np.zeros((U0.shape[0], ncol))
    U[:, 0] = U0
    U[:, 1] = U0 + dt * V0 + 0.5 * dt**2 * (F - A.dot(U0))

    for n in range(1, ncol-1):
        F = get_fdtd_second_member(c, dx, dy, nx, ny, s(x, y, n * dt), k(x, y, n * dt))
        U[:, n + 1] = 2 * U[:, n] - U[:, n - 1] + dt**2 * (F - A.dot(U[:, n]))

        # At the corners, the second-order temporal derivatives are zero
        for i in [0, ny * (nx + 1), nx + ny * (nx + 1), nx]:
            U[i, n + 1] = 2 * U[i, n] - U[i, n - 1]

    return U