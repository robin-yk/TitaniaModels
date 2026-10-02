"""Sphere and rutile(110) site capacities, in umol per gram of initial TiO2.

Every O and Ti site of the particle is counted once: explicit near-surface
sites first, then the bulk remainder. Values come from
specification/parameter_registry.csv (a test checks they agree).
"""
import math

NA = 6.02214076e23
KB = 8.617333262e-5          # eV/K
A_NM, C_NM, U = 0.4594, 0.2958, 0.30479
RHO = 4.2487                 # g/cm3
M = 79.866                   # g/mol TiO2

D110 = A_NM / math.sqrt(2)                  # Ti plane spacing along [110], nm
CELL = math.sqrt(2) * A_NM * C_NM           # (110) 1x1 cell area, nm2
H_BRI = D110 - math.sqrt(2) * U * A_NM      # bridging O above its Ti plane, nm

O_TOTAL = 2e6 / M
TI_TOTAL = 1e6 / M

# One O-Ti2O2-O trilayer per 1x1 cell: site class -> (multiplicity, depth below
# the bridging O of the same trilayer).
TRILAYER_O = {'BRI': (1, 0.0), 'IPL': (2, H_BRI), 'SBR': (1, 2 * H_BRI)}
TRILAYER_TI = (2, H_BRI)


def bridging_capacity(d_nm):
    """Bridging O per gram on a smooth sphere of diameter d_nm."""
    area_m2_g = 6.0 / (RHO * d_nm * 1e-9) / 1e6
    return area_m2_g * 1e18 / CELL / NA * 1e6


def layer_sites(d_nm, n_layers):
    """Explicit trilayer sites: list of (layer, site, O capacity, depth nm) and
    per-layer Ti capacity. Area shrinks as (1 - z/R)^2 with depth."""
    R, cb = d_nm / 2, bridging_capacity(d_nm)
    o, ti = [], []
    for k in range(1, n_layers + 1):
        top = (k - 1) * D110
        for site, (m, dz) in TRILAYER_O.items():
            z = top + dz
            o.append((k, site, m * cb * (1 - z / R) ** 2, z))
        m, dz = TRILAYER_TI
        ti.append(m * cb * (1 - (top + dz) / R) ** 2)
    return o, ti


def shell_fraction(R, z_lo, z_hi):
    """Volume fraction of the sphere between depths z_lo and z_hi."""
    return (1 - z_lo / R) ** 3 - (1 - z_hi / R) ** 3


def particle_mass(d_nm):
    """Mass of one smooth sphere, g."""
    return RHO * (math.pi / 6) * d_nm ** 3 * 1e-21


def oxygen_shells(rmax_nm=1.0):
    """O-O neighbour shells of rutile out to rmax: [(r_nm, count per O)].

    Rutile has one oxygen Wyckoff site (4f), so every O sees the same shells.
    """
    import numpy as np
    cell = np.diag([A_NM, A_NM, C_NM])
    basis = np.array([[U, U, 0], [-U, -U, 0], [.5 + U, .5 - U, .5], [.5 - U, .5 + U, .5]])
    n = int(rmax_nm / min(A_NM, C_NM)) + 2
    g = np.arange(-n, n + 1)
    shifts = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
    vec = ((basis[None, :, :] + shifts[:, None, :]).reshape(-1, 3) - basis[0]) @ cell
    r = np.sqrt((vec ** 2).sum(1))
    r = np.round(r[(r > 1e-9) & (r <= rmax_nm + 1e-12)], 9)
    vals, counts = np.unique(r, return_counts=True)
    return [(float(v), int(c)) for v, c in zip(vals, counts)]


# Bulk vacancy hop: Iddir 2007 step C (equatorial to apical O) is the
# octahedron edge between an apical and an equatorial O, the second O-O shell
# (8 neighbours). D = z lambda^2 Gamma / 6 on that network.
def bulk_hop():
    r, z = oxygen_shells(0.3)[1]
    return r, z
