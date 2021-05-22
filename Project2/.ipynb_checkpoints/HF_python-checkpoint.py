import sys

sys.path.append("/OneDrive/Skrivebord/quantum-systems/")
import quantum_systems as qs
import numpy as np
np.set_printoptions(precision=2) # 2 decimals when printing arrays
np.set_printoptions(suppress=True) # No scientific notation for small numbers
import matplotlib.pyplot as plt

# from C:\Users\brbre\OneDrive\Skrivebord\quantum_systems import ODQD, GeneralOrbitalSystem # library developed by Øyvind Schøyen and others, https://github.com/Schoyen/quantum-systems

l_0 = 10                # number of basis functions
grid_length = 10        # compute from x = -10 to x = 10 in 1D
num_grid_points = 2001
omega = 0.25            # strength of harmonic oscillator potential
n = 2                   # number of particles

# this sets up the harmonic oscillator basis functions and integrals between them
odho = qs.ODQD(l_0, grid_length, num_grid_points, a = 0.25, alpha = 1, potential = qs.ODQD.HOPotential(omega))

# this makes a spin up and spin down variant of the odho basis function and sets up the integrals between them
system = qs.GeneralOrbitalSystem(n = n, basis_set=odho)
l = system.l

print(f"l = {system.l}")
print(f"grid shape = {system._basis_set.grid.shape}")
print(f"h shape = {system.h.shape}")
print(f"u shape = {system.u.shape}")
print(f"x shape = {system.position.shape}")
print(f"spf shape = {system.spf.shape}")

np.diag(system.h.real)

def getP(C):
    l = C.shape[0]
    P = np.zeros((l, l), dtype=complex)

    for g in range(l):
        for d in range(l):
            for i in range(n):
                P[d, g] += np.conj(C[g, i]) * C[d, i]
    return P

def getF(C):
    l = C.shape[0]
    F = np.zeros((l, l), dtype=complex)
    P = getP(C)

    for b in range(l):
        for a in range(l):
            F[b, a] += system.h[b, a]
            for g in range(l):
                for d in range(l):
                    F[b, a] += P[d, g] * system.u[b, g, a, d]
    return F

C = np.eye(l, dtype=complex)

for i in range(20):
    F = getF(C)
    vals, C = np.linalg.eigh(F)


gs = C[:, 0].real
print(gs)

gs_vals = np.zeros(num_grid_points, dtype = complex)
for i in range(l):
    gs_vals += gs[i] * system.spf[i]

plt.plot(system.grid, np.abs(gs_vals)**2 * 2001, label=r"$\psi_{GS}$",)
plt.xlim(-6, 6)

plt.grid()
plt.legend()
plt.show()

system.change_basis(C)
print(system.compute_reference_energy())
odho = qs.ODQD(l_0, grid_length, num_grid_points, a = 0.25, alpha = 1, potential = ODQD.HOPotential(omega))
system = qs.GeneralOrbitalSystem(n = n, basis_set=odho)
