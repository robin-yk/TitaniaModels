"""Constrained free-energy minimum over a finite state table.

A domain d holds C_d umol of cells per gram; each cell is in one state a with
energy E_a, degeneracy g_a, v_a oxygen vacancies and e_a Ti3+ electrons. Its
charge is z_a = 2 v_a - e_a. Populations x_a = C_d p_a minimise

    F = sum x_a E_a + kT sum x_a ln(x_a / (C_d g_a)) + 1/2 q^T H q

subject to sum v_a x_a = N (measured inventory) and one charge closure:

  NEUTRAL  electrons implicit; no charge terms (domains carry no region).
  LOCAL    every region r is neutral, sum_{a in r} z_a x_a = 0.
  GLOBAL   only the particle is neutral, sum e_a x_a = 2 N. Shell charges
           q_l = sum_{a in l} z_a x_a act through H_lm = K / max(r_l, r_m),
           the potential of concentric charged shells in a medium of
           permittivity eps (K = e^2/(4 pi eps0 eps) per umol/g of charge).

The minimiser has Gibbs form, p_a ~ g_a exp(-E_a/kT + (multipliers) . (v, e, z)).
The unknowns are the multipliers: mu (and chi per region) for NEUTRAL/LOCAL,
mu, the electron potential and one electrostatic potential per shell for
GLOBAL. They minimise a convex dual.

NEUTRAL/LOCAL: region charge rises monotonically with its own chi and the
inventory with mu once regions are neutral, so both are bracketed scalar roots
(inner Newton-bisection for all chi at once, outer for mu). They cannot fail
to converge.

GLOBAL: the dual adds 1/2 u^T (beta H)^-1 u for u = phi/kT. H is c*min(s_l,s_m)
with s = 1/r, whose inverse is tridiagonal: the capacitance matrix of nested
spherical capacitors, 1/(1/r_i - 1/r_{i+1}) between neighbours and r_1 to
infinity. The dual Hessian is therefore tridiagonal plus a rank-two border, and
each damped Newton step costs O(shells). It starts from the LOCAL solution on
the same regions, which is the eps -> 0 limit of GLOBAL.

A domain may put its electrons on a second, adjacent shell (shell2): its
vacancy charge 2v sits on shell, its electron charge -e on shell2. The two
shells are neighbours in the radial order, so the dual Hessian stays
tridiagonal plus the border.

Frozen populations (fixed=...) enter as a vacancy count and fixed charges per
region or shell; they are not optimised.
"""
import numpy as np
from scipy.linalg import solve_banded

from .particle import KB, NA

COULOMB = 1.43996454784255     # e^2 / (4 pi eps0), eV nm
BOUND, MAX_ITER, COLLAPSE = 300.0, 400, 1e-14


class Model:
    """domains: list of dicts with C, E, v and optional e, g, region, shell, tag.

    electrostatics: dict(radii={shell id: r_nm}, eps=..., mass_g=...) for GLOBAL.
    """

    def __init__(self, domains, electrostatics=None):
        S = max(len(d['E']) for d in domains)
        D = len(domains)
        self.C = np.array([d['C'] for d in domains], float)
        self.E = np.zeros((D, S)); self.v = np.zeros((D, S)); self.e = np.zeros((D, S))
        self.logg = np.full((D, S), -np.inf)
        for i, d in enumerate(domains):
            n = len(d['E'])
            self.E[i, :n] = d['E']; self.v[i, :n] = d['v']
            self.e[i, :n] = d.get('e', [0] * n)
            self.logg[i, :n] = np.log(np.asarray(d.get('g', [1] * n), float))
        if np.any(self.C <= 0):
            raise ValueError('every domain needs positive capacity')
        self.z = 2 * self.v - self.e
        # Charge on the domain's own shell (zA) and on its electron shell (zB).
        split = np.array([d.get('shell2') is not None for d in domains])
        self.zA = np.where(split[:, None], 2 * self.v, self.z)
        self.zB = np.where(split[:, None], -self.e, 0.0)
        region = [d.get('region', -1) for d in domains]
        # Outside any region the electrons are implicit and carry no charge.
        self.q = np.where(np.array(region)[:, None] < 0, 0.0, self.z)
        self.region_ids = sorted({r for r in region if r >= 0})
        ids = {r: k for k, r in enumerate(self.region_ids)}
        self.region = np.array([ids.get(r, -1) for r in region])
        self.nreg = len(ids)
        self.tags = [d.get('tag') for d in domains]
        self.es = None
        if electrostatics:
            self._setup_electrostatics(domains, electrostatics)

    # ---------------------------------------------------------------- setup
    def _setup_electrostatics(self, domains, es):
        if self.nreg == 0 or np.any(self.region < 0):
            raise ValueError('GLOBAL needs every domain in a region: its start is the LOCAL solution')
        radii = es['radii']
        order = sorted(radii, key=lambda s: -radii[s])          # outermost first
        r = np.array([radii[s] for s in order], float)
        if np.any(np.diff(r) >= 0):
            raise ValueError('electrostatic shells need distinct radii')
        pos = {s: k for k, s in enumerate(order)}
        self.shell = np.array([pos[d['shell']] for d in domains])
        self.shell2 = np.array([pos[d.get('shell2', d['shell'])] for d in domains])
        if np.any(np.abs(self.shell2 - self.shell) > 1):
            raise ValueError('a split domain needs adjacent shells')
        reg_of_shell = {}
        for i, s in enumerate(np.r_[self.shell, self.shell2]):
            if reg_of_shell.setdefault(s, self.region[i % len(self.shell)]) != self.region[i % len(self.shell)]:
                raise ValueError('a shell must lie inside one region')
        s_inv = 1.0 / r
        c = 1.0 / np.diff(np.r_[0.0, s_inv])                     # capacitances
        main = c + np.r_[c[1:], 0.0]
        self.es = dict(r=r, order=order, K=COULOMB / es['eps'] * 1e-6 * NA * es['mass_g'],
                       main=main, off=-c[1:], s=s_inv, L=len(r),
                       region_of_shell=np.array([reg_of_shell[k] for k in range(len(r))]))

    def potentials(self, q):
        """phi_l = sum_m H_lm q_m (eV per elementary charge), O(L) via cumsums."""
        s, K = self.es['s'], self.es['K']
        # min(s_l, s_m): s ascending with l, so inner shells (m > l) give s_l.
        head = np.cumsum(s * q)                 # sum_{m <= l} s_m q_m
        tail = np.cumsum(q[::-1])[::-1]         # sum_{m >= l} q_m
        return K * (head - s * q + s * tail)

    def _Pmul(self, u, kT):
        es = self.es
        out = es['main'] * u
        out[:-1] += es['off'] * u[1:]
        out[1:] += es['off'] * u[:-1]
        return kT / es['K'] * out

    # ------------------------------------------------------------- NEUTRAL / LOCAL
    def _p(self, y_mu, chi, beta):
        L = self.logg - beta * self.E + y_mu * self.v
        if self.nreg:
            L = L + np.r_[0.0, chi][self.region + 1][:, None] * self.q
        w = np.exp(L - L.max(axis=1)[:, None])
        return w / w.sum(axis=1)[:, None]

    def _neutralise(self, y_mu, chi, beta, tol, N, qfix):
        """chi that makes every region neutral at this mu, with d(inventory)/d(mu)."""
        R, inr = self.nreg, self.region >= 0
        ri, C = self.region[inr], self.C[inr]
        lo, hi = np.full(R, -BOUND), np.full(R, BOUND)
        dx_old = np.full(R, 2 * BOUND)
        for _ in range(MAX_ITER):
            p = self._p(y_mu, chi, beta)[inr]
            mq = (p * self.q[inr]).sum(1)
            h = np.bincount(ri, C * mq, R) + qfix
            dh = np.bincount(ri, C * ((p * self.q[inr] ** 2).sum(1) - mq ** 2), R)
            if np.all((np.abs(h) <= tol * N) | (hi - lo <= COLLAPSE * np.maximum(1, np.abs(chi)))):
                break
            hi = np.where(h > 0, chi, hi); lo = np.where(h > 0, lo, chi)
            step = chi - h / np.maximum(dh, 1e-300)
            # Newton only while it lands inside the bracket and at least halves
            # the step before last; otherwise bisect (a sharp sigmoid can make
            # Newton bounce between the two bracket ends).
            newton = (step > lo) & (step < hi) & (np.abs(chi - step) < 0.5 * dx_old)
            nxt = np.where(newton, step, 0.5 * (lo + hi))
            dx_old, chi = np.abs(nxt - chi), nxt
        else:
            raise RuntimeError('local neutrality did not converge')
        mv = (p * self.v[inr]).sum(1)
        b = np.bincount(ri, C * ((p * self.v[inr] * self.q[inr]).sum(1) - mv * mq), R)
        return chi, float(np.sum(b * b / dh))

    def _solve_local(self, N, T, tol, fixed):
        beta = 1.0 / (KB * T)
        Nf = N - fixed.get('v', 0.0)
        qfix = np.zeros(self.nreg)
        for r, qv in fixed.get('q_region', {}).items():
            qfix[self.region_ids.index(r)] += qv
        if not 0 < Nf < float(np.sum(self.C * self.v.max(axis=1))):
            raise ValueError('inventory outside the represented capacity')
        y, chi, lo, hi = 0.0, np.zeros(self.nreg), -BOUND, BOUND
        dx_old = 2 * BOUND
        for it in range(1, MAX_ITER + 1):
            schur = 0.0
            if self.nreg:
                chi, schur = self._neutralise(y, chi, beta, tol, Nf, qfix)
            p = self._p(y, chi, beta)
            mv = (p * self.v).sum(1)
            g = float(self.C @ mv) - Nf
            # Stop at the tolerance, or when the bracket has shrunk to rounding.
            if abs(g) <= tol * Nf or hi - lo <= COLLAPSE * max(1, abs(y)):
                break
            dg = float(self.C @ ((p * self.v ** 2).sum(1) - mv ** 2)) - schur
            hi, lo = (y, lo) if g > 0 else (hi, y)
            step = y - g / max(dg, 1e-300)
            nxt = step if lo < step < hi and abs(y - step) < 0.5 * dx_old else 0.5 * (lo + hi)
            dx_old, y = abs(nxt - y), nxt
        else:
            raise RuntimeError('equilibrium did not converge')
        x = self.C[:, None] * p
        inr = self.region >= 0
        chg = np.bincount(self.region[inr], (x * self.q).sum(1)[inr], self.nreg) + qfix
        return x, y, chi, it, float((x * self.v).sum() - Nf), float(np.abs(chg).max(initial=0.0))

    # ----------------------------------------------------------------- GLOBAL
    def _global_eval(self, y, u, beta, Nf, Ef, qfix, need_hess=True):
        es = self.es
        L = (self.logg - beta * self.E + y[0] * self.v + y[1] * self.e
             - u[self.shell][:, None] * self.zA - u[self.shell2][:, None] * self.zB)
        Lmax = L.max(axis=1)
        w = np.exp(L - Lmax[:, None]); Zs = w.sum(axis=1); p = w / Zs[:, None]
        kT = 1.0 / beta
        Pu = self._Pmul(u, kT)
        psi = float(self.C @ (Lmax + np.log(Zs))) - Nf * y[0] - Ef * y[1] + 0.5 * u @ Pu - u @ qfix
        if not need_hess:
            return psi, p
        C, sh, sh2, Lsh = self.C, self.shell, self.shell2, es['L']
        zA, zB = self.zA, self.zB
        mv = (p * self.v).sum(1); me = (p * self.e).sum(1)
        mA = (p * zA).sum(1); mB = (p * zB).sum(1)
        g = np.empty(2 + Lsh)
        g[0] = C @ mv - Nf
        g[1] = C @ me - Ef
        g[2:] = -np.bincount(sh, C * mA, Lsh) - np.bincount(sh2, C * mB, Lsh) + Pu - qfix
        cov = lambda a, b, ma, mb: C * ((p * a * b).sum(1) - ma * mb)   # noqa: E731
        Hyy = np.array([[np.sum(cov(self.v, self.v, mv, mv)), np.sum(cov(self.v, self.e, mv, me))],
                        [0.0, np.sum(cov(self.e, self.e, me, me))]])
        Hyy[1, 0] = Hyy[0, 1]
        B = np.stack([-np.bincount(sh, cov(self.v, zA, mv, mA), Lsh) - np.bincount(sh2, cov(self.v, zB, mv, mB), Lsh),
                      -np.bincount(sh, cov(self.e, zA, me, mA), Lsh) - np.bincount(sh2, cov(self.e, zB, me, mB), Lsh)],
                     axis=1)
        AB = cov(zA, zB, mA, mB)
        same = sh == sh2
        diag = (np.bincount(sh, cov(zA, zA, mA, mA), Lsh) + np.bincount(sh2, cov(zB, zB, mB, mB), Lsh)
                + 2 * np.bincount(sh[same], AB[same], Lsh) + kT / es['K'] * es['main'])
        off = kT / es['K'] * es['off'] + np.bincount(np.minimum(sh, sh2)[~same], AB[~same], Lsh)[:-1]
        return psi, p, g, Hyy, B, diag, off

    def _solve_global(self, N, T, tol, fixed):
        beta = 1.0 / (KB * T)
        es = self.es
        # Start: the LOCAL solution on the same regions (eps -> 0 limit).
        qreg = {}
        for s_id, qv in fixed.get('q_shell', {}).items():
            k = es['order'].index(s_id)
            r = self.region_ids[es['region_of_shell'][k]]
            qreg[r] = qreg.get(r, 0.0) + qv
        _, yL, chi, _, _, _ = self._solve_local(N, T, tol, dict(v=fixed.get('v', 0.0), q_region=qreg))
        chi_sh = chi[es['region_of_shell']]
        u = chi_sh[0] - chi_sh
        y = np.array([yL + 2 * chi_sh[0], -chi_sh[0]])
        Nf = N - fixed.get('v', 0.0)
        Ef = 2.0 * N
        qfix = np.zeros(es['L'])
        for s_id, qv in fixed.get('q_shell', {}).items():
            qfix[es['order'].index(s_id)] += qv
        for it in range(1, MAX_ITER + 1):
            psi, p, g, Hyy, B, diag, off = self._global_eval(y, u, beta, Nf, Ef, qfix)
            if max(abs(g[0]), abs(g[1]), np.abs(g[2:]).max()) <= tol * Nf:
                break
            ab = np.zeros((3, es['L'])); ab[0, 1:] = off; ab[1] = diag; ab[2, :-1] = off
            sol = solve_banded((1, 1), ab, np.column_stack([g[2:], B]))
            w, X = sol[:, 0], sol[:, 1:]
            dy = np.linalg.solve(Hyy - B.T @ X, -g[:2] + B.T @ w)
            du = -w - X @ dy
            slope = g[:2] @ dy + g[2:] @ du
            # Newton decrement at rounding level: nothing left to gain.
            if -slope <= 1e-26 * max(1.0, abs(psi)):
                break
            t = 1.0
            while True:
                psi_new, _ = self._global_eval(y + t * dy, u + t * du, beta, Nf, Ef, qfix, False)
                if psi_new <= psi + 1e-4 * t * slope + 1e-14 * abs(psi) or t < 1e-12:
                    break
                t *= 0.5
            y, u = y + t * dy, u + t * du
        else:
            raise RuntimeError('GLOBAL equilibrium did not converge')
        x = self.C[:, None] * p
        q = (np.bincount(self.shell, (x * self.zA).sum(1), es['L'])
             + np.bincount(self.shell2, (x * self.zB).sum(1), es['L']) + qfix)
        poisson = float(np.abs(self.potentials(q) - u / beta).max())
        return x, y, u, it, float((x * self.v).sum() - Nf), float((x * self.e).sum() - Ef), q, poisson

    # ------------------------------------------------------------------ solve
    def solve(self, N, T, tol=1e-12, fixed=None):
        fixed = fixed or {}
        if self.es is None:
            x, y, chi, it, inv, chg = self._solve_local(N, T, tol, fixed)
            return Solution(self, x, T, y * KB * T, it, inv, charge_residual=chg, chi=chi)
        x, y, u, it, inv, el, q, poisson = self._solve_global(N, T, tol, fixed)
        return Solution(self, x, T, y[0] * KB * T, it, inv, charge_residual=el,
                        phi=u * KB * T, mu_e=y[1] * KB * T, shell_charge=q,
                        poisson_residual=poisson)


class Solution:
    def __init__(self, model, x, T, mu_eV, iterations, inventory_residual, **extra):
        self.model, self.x, self.T, self.mu_eV = model, x, T, mu_eV
        self.iterations, self.inventory_residual = iterations, inventory_residual
        self.chi = extra.pop('chi', None)
        self.phi = extra.pop('phi', None)
        self.__dict__.update(extra)

    def occupancy(self, i):
        """Mean vacancies per cell of domain i."""
        return float(self.x[i] @ self.model.v[i]) / self.model.C[i]

    def vacancies(self, tag):
        """Vacancy population (umol/g) summed over domains with this tag."""
        m = self.model
        return float(sum((self.x[i] * m.v[i]).sum() for i, t in enumerate(m.tags) if t == tag))
