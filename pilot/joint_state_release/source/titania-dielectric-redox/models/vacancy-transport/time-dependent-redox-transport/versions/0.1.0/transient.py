"""Source-driven, fast-electron transport pilot on the unmodified state tables.

The radial mobility is a declared effective transport closure, not a site-specific
DFT transition network. The free energy, capacities and electrostatics are reused.
"""
import argparse
import hashlib
import json
import time
import warnings

import numpy as np
from scipy.linalg import solve, lstsq, LinAlgWarning
from scipy.special import expit

# pilot accepts a positional checkout; avoid passing this script's flags to it.
import sys
_argv = sys.argv
sys.argv = [sys.argv[0]]
from pilot import HERE, build, quantities, source, T, pt
sys.argv = _argv


def bernoulli(z):
    z = np.asarray(z)
    small = np.abs(z) < 1e-4
    # Clip only evaluation arguments in unused branches to avoid overflow.
    a = np.clip(z, -700, 700)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        b = a / np.expm1(a)
        bp = (np.expm1(a) - a * np.exp(a)) / np.expm1(a)**2
    b = np.where(small, 1 - z/2 + z*z/12 - z**4/720, b)
    bp = np.where(small, -.5 + z/6 - z**3/180, bp)
    return b, bp


class Transport:
    def __init__(self, n_bulk=80, diffusivity_nm2_s=None):
        self.m, self.lay = build(n_bulk)
        m = self.m
        self.L = L = m.es['L']
        self.kT = pt.KB*T
        self.size = 2*L + 1
        q0 = quantities(m, self.lay, np.zeros_like(m.E), np.zeros(L))
        self.cap = np.array(q0['oxygen_capacity_umol_g'])
        self.P = self.kT / m.es['K'] * (np.diag(m.es['main']) +
                 np.diag(m.es['off'], 1) + np.diag(m.es['off'], -1))
        self.features = [m.v, m.e, -m.zA, -m.zB]
        self.indices = [m.shell, np.full(len(m.C), L), L+1+m.shell, L+1+m.shell2]
        self.D = (1e13*.3**2/6*np.exp(-1.10/self.kT)
                  if diffusivity_nm2_s is None else diffusivity_nm2_s)
        r = m.es['r']
        face = (r[:-1]+r[1:])/2
        density = pt.O_TOTAL/(4*np.pi/3*450**3)
        self.G = self.D * density*4*np.pi*face**2 / (-np.diff(r))
        # The atomic part represents 75% (110) area; the bulk is a full sphere.
        self.G[:12] *= .75

    def state(self, Y, hessian=True):
        m, L = self.m, self.L
        alpha, ye, u = Y[:L], Y[L], Y[L+1:]
        logw = (m.logg - m.E/self.kT + alpha[m.shell,None]*m.v + ye*m.e
                - u[m.shell,None]*m.zA-u[m.shell2,None]*m.zB)
        w = np.exp(logw-logw.max(1)[:,None])
        p = w/w.sum(1)[:,None]
        x = m.C[:,None]*p
        means = [(p*f).sum(1) for f in self.features]
        n = np.bincount(m.shell, m.C*means[0], minlength=L)
        ne = m.C@means[1]
        charge = (-np.bincount(m.shell, m.C*means[2], minlength=L)
                  -np.bincount(m.shell2, m.C*means[3], minlength=L))
        g = np.r_[n, ne, self.P@u-charge]
        if not hessian:
            return g, x, charge
        H = np.zeros((self.size,self.size))
        for a in range(4):
            for b in range(a,4):
                cv = m.C*((p*self.features[a]*self.features[b]).sum(1)-means[a]*means[b])
                np.add.at(H,(self.indices[a],self.indices[b]),cv)
                if a != b:
                    np.add.at(H,(self.indices[b],self.indices[a]),cv)
        H[L+1:,L+1:] += self.P
        return g, x, charge, H

    def flux(self, Y, n, H=None):
        L = self.L
        th = np.clip(n/self.cap, 1e-280, 1-1e-14)
        eta = Y[:L]-np.log(th/(1-th))
        delta = np.diff(eta)
        bp, dbp = bernoulli(delta)
        bm, dbm = bernoulli(-delta)
        a = th[:-1]*(1-th[1:])
        c = th[1:]*(1-th[:-1])
        dalpha = Y[:L-1]-Y[1:L]
        # Stable flux near equilibrium (ordinary forward-minus-reverse cancels).
        pos = dalpha >= 0
        J = np.empty(L-1)
        J[pos] = self.G[pos]*bp[pos]*a[pos]*(-np.expm1(-dalpha[pos]))
        J[~pos] = self.G[~pos]*bm[~pos]*c[~pos]*np.expm1(dalpha[~pos])
        if H is None:
            return J
        jd = self.G*(a*dbp+c*dbm)
        jti = self.G*((1-th[1:])*bp+th[1:]*bm)+jd/(th[:-1]*(1-th[:-1]))
        jtj = -self.G*(th[:-1]*bp+(1-th[:-1])*bm)-jd/(th[1:]*(1-th[1:]))
        jac = jti[:,None]*(H[:L-1]/self.cap[:-1,None]) + jtj[:,None]*(H[1:L]/self.cap[1:,None])
        idx = np.arange(L-1)
        jac[idx,idx] -= jd
        jac[idx,idx+1] += jd
        return J, jac

    def residual(self, Y, old, added, dt, target, jacobian=True):
        out = self.state(Y,jacobian)
        g,x,q = out[:3]
        L = self.L
        diff = g[:L]-old
        if jacobian:
            H = out[3]
            J,dJ = self.flux(Y,g[:L],H)
        else:
            J = self.flux(Y,g[:L])
        # Integrate inward from the closed centre. Surface-side cumulative sums
        # cancel the full supplied inventory to resolve an almost empty tail.
        # The independent total row supplies the surface source balance.
        residual = np.r_[-np.cumsum(diff[:0:-1])[::-1]+dt*J, g[:L].sum()-target,
                         g[L]-2*target, g[L+1:]]
        if not jacobian:
            return residual,x,q,J
        jac = np.vstack((-np.cumsum(H[L-1:0:-1],axis=0)[::-1]+dt*dJ,
                         H[:L].sum(0),H[L:]))
        return residual,x,q,J,jac

    def guess(self, inventory):
        s = self.m.solve(inventory,T)
        return np.r_[np.full(self.L,s.mu_eV/self.kT),s.mu_e/self.kT,s.phi/self.kT]

    def step(self,Y,old,added,dt,target):
        scale = max(target,1.)
        tolerance = np.full(self.size, 2e-8*scale)
        # Tight totals prevent mass error masquerading as a free-energy change.
        # Flux balances have a larger floating-point cancellation floor.
        tolerance[self.L-1:self.L+1] = 1e-11
        if getattr(self, 'null_tail', False):
            tolerance[self.L-1:self.L+1] = 1e-9
        tolerance[self.L+1:] = 1e-10*scale
        for it in range(55):
            f,x,q,J,A = self.residual(Y,old,added,dt,target)
            err = np.max(np.abs(f))
            merit = np.max(np.abs(f)/tolerance)
            if merit < 1:
                return Y,x,q,J,it+1,err
            # Row equilibration improves the linear algebra for stiff diffusion.
            rows = np.maximum(np.max(np.abs(A),axis=1),1e-16)
            with warnings.catch_warnings():
                warnings.simplefilter('error',LinAlgWarning)
                try:
                    d = solve(A/rows[:,None],-f/rows,assume_a='gen',check_finite=False)
                except (LinAlgWarning, np.linalg.LinAlgError):
                    # Near-empty deep shells have indeterminate log activities.
                    # A minimum-norm Newton step leaves null directions unchanged;
                    # acceptance still uses the original full residual tolerances.
                    if not getattr(self, 'null_tail', False):
                        raise
                    d = lstsq(A/rows[:,None],-f/rows,cond=1e-12,
                              check_finite=False,lapack_driver='gelsd')[0]
            length = 1.
            # Prevent trial exponent differences overflowing at low occupancy.
            length = min(length, 20/max(np.max(np.abs(d)),1e-30))
            while length > 1e-8:
                Yn = Y+length*d
                fn = self.residual(Yn,old,added,dt,target,False)[0]
                if np.max(np.abs(fn)/tolerance) < merit*(1-1e-4*length):
                    Y = Yn
                    break
                length *= .5
            else:
                raise RuntimeError(f'line search failed dt={dt} error={err} '
                                   f'row={np.argmax(np.abs(f)/tolerance)} merit={merit}')
        raise RuntimeError(f'Newton failed dt={dt} error={err} '
                           f'row={np.argmax(np.abs(f)/tolerance)} merit={merit}')

    def free_energy(self,x,q):
        m = self.m
        use = x>0
        p = x/m.C[:,None]
        term = np.zeros_like(x)
        term[use] = x[use]*(np.log(p[use])-m.logg[use])
        return float((x*m.E).sum()+self.kT*term.sum()+.5*q@self.m.potentials(q))

    def row(self,Y,x,q,time_s,source_total,dt,it,err):
        row = quantities(self.m,self.lay,x,Y[self.L+1:]*self.kT)
        n = np.array(row['vacancy_umol_g'])
        J = self.flux(Y,n)
        row.update(time_s=time_s,source_total_umol_g=source_total,dt_s=dt,iterations=it,
                   equation_max_abs_umol_g=err,
                   vacancy_chemical_potential_eV=(self.kT*Y[:self.L]).tolist(),
                   free_energy_umol_eV_g=self.free_energy(x,q),
                   transport_dissipation_umol_eV_g_s=float(self.kT*J@(Y[:self.L-1]-Y[1:self.L])),
                   inward_flux_umol_g_s=J.tolist(),
                   poisson_error_eV=float(np.max(np.abs(self.m.potentials(q)-Y[self.L+1:]*self.kT))))
        return row


def jacobian_check(tr,Y):
    g=tr.state(Y)[0]
    args=(g[:tr.L]*.9,g[:tr.L].sum()*.1,1.,g[:tr.L].sum())
    f,*_,A=tr.residual(Y,*args)
    rng=np.random.default_rng(19)
    d=rng.normal(size=Y.size)
    d/=np.linalg.norm(d)
    h=1e-5
    fd=(tr.residual(Y+h*d,*args,False)[0]-tr.residual(Y-h*d,*args,False)[0])/(2*h)
    rel=float(np.linalg.norm(fd-A@d)/np.linalg.norm(A@d))
    if rel>1e-6:
        raise AssertionError(f'Jacobian relative error {rel}')
    return rel


def run(args):
    tr=Transport(args.bulk,args.D)
    tr.null_tail=getattr(args,'null_tail',False)
    amount,src=source()
    old=np.zeros(tr.L)  # Exact zero initial inventory; no numerical seed stock.
    t=0.
    dt=.25
    Y=tr.guess(amount(dt))
    checks={'analytic_jacobian_relative_error':jacobian_check(tr,Y)}
    rows=[]
    log=[]
    audit={'max_local_balance_umol_g': 0., 'max_vacancy_balance_umol_g': 0.,
           'max_electron_balance_umol_g': 0., 'min_transport_dissipation': float('inf')}
    checkpoints=[1.,10.,60.,180.,600.,1800.]
    start=time.monotonic()
    while t<args.end-1e-8:
        stop=min([c for c in checkpoints if c>t+1e-8]+[args.end])
        dt=min(dt,args.max_step,args.end-t,stop-t)
        added=amount(t+dt)-amount(t)
        try:
            Yn,x,q,J,it,err=tr.step(Y,old,added,dt,amount(t+dt))
        except (RuntimeError,LinAlgWarning) as exc:
            dt*=.5
            print(f'Retrying t={t:g}, dt={dt:g}: {exc}',flush=True)
            if dt<1e-6:
                raise
            continue
        Y=Yn
        t+=dt
        state=tr.state(Y,False)[0]
        local=state[:tr.L]-old
        local[0]-=added
        local[:-1]+=dt*J
        local[1:]-=dt*J
        audit['max_local_balance_umol_g']=max(audit['max_local_balance_umol_g'],float(np.max(np.abs(local))))
        audit['max_vacancy_balance_umol_g']=max(audit['max_vacancy_balance_umol_g'],float(abs(state[:tr.L].sum()-amount(t))))
        audit['max_electron_balance_umol_g']=max(audit['max_electron_balance_umol_g'],float(abs(state[tr.L]-2*amount(t))))
        audit['min_transport_dissipation']=min(audit['min_transport_dissipation'],float(tr.kT*J@(Y[:tr.L-1]-Y[1:tr.L])))
        old=state[:tr.L]
        log.append([t,dt,it,err])
        if abs(t-stop)<1e-7:
            row=tr.row(Y,x,q,t,amount(t),dt,it,err)
            rows.append(row)
            print(json.dumps(dict(t=t,BRI_pct=row['theta_BRI_pct'],bulk=row['pools_umol_g']['bulk'],
                                  surface=row['pools_umol_g']['bridging'],N=old.sum(),
                                  error=err,elapsed_s=time.monotonic()-start)),flush=True)
        dt=min(args.max_step,dt*1.7)
    checks.update(audit)
    checks['final_analytic_jacobian_relative_error']=jacobian_check(tr,Y)
    # Test relaxation with source switched off at the final inventory.
    if args.relax:
        Ne=amount(t)
        eq=tr.m.solve(Ne,T)
        eqrow=quantities(tr.m,tr.lay,eq.x,eq.phi)
        checks['equilibrium_stationary_flux']=float(np.max(np.abs(tr.flux(tr.guess(Ne),np.array(eqrow['vacancy_umol_g'])))))
        relax=[]
        # Source removal can cross a local state transition. Apply the same
        # step-rejection policy used for the source-driven trajectory.
        h=.01
        closed_end=t+111.11
        closed_increments=[]
        Fprevious=rows[-1]['free_energy_umol_eV_g']
        while t<closed_end-1e-8:
            h=min(h,args.max_step,closed_end-t)
            try:
                Yn,x,q,J,it,err=tr.step(Y,old,0.,h,Ne)
            except (RuntimeError,LinAlgWarning):
                h*=.5
                if h<1e-6:
                    raise
                continue
            Y=Yn
            old=tr.state(Y,False)[0][:tr.L]
            t+=h
            r=tr.row(Y,x,q,t,Ne,h,it,err)
            closed_increments.append(r['free_energy_umol_eV_g']-Fprevious)
            Fprevious=r['free_energy_umol_eV_g']
            relax.append(r)
            h*=1.7
        checks['relaxed_shell_population_L1_error_umol_g']=float(np.sum(np.abs(old-np.array(eqrow['vacancy_umol_g']))))
        checks['relaxed_BRI_pct_error']=relax[-1]['theta_BRI_pct']-eqrow['theta_BRI_pct']
        checks['closed_free_energy_increments'] = closed_increments
    else:
        relax=[]
    out=dict(model='original LI_SBR1 free energy with effective radial SG transport',
             balance_form='centre-side cumulative finite-volume balances plus total inventory',
             numerical_tolerances={'transport': '2e-8 * max(N,1) umol/g',
                                   'global_V_and_e': ('1e-9 umol/g' if tr.null_tail else '1e-11 umol/g'),
                                   'Poisson': '1e-10 * max(N,1) umol/g'},
             code_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in ('pilot.py','transient.py')},
             T_K=T,diameter_nm=900.,initial_vacancy_umol_g=0.,initial_electron_umol_g=0.,
             n_bulk=args.bulk,D_nm2_s=tr.D,max_step_s=args.max_step,source=src,
             null_tail_minimum_norm=tr.null_tail,
             checks=checks,checkpoints=rows,closed_relaxation=relax,step_log=log,
             elapsed_s=time.monotonic()-start)
    path=HERE/args.output
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(f'Saved {path}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--bulk',type=int,default=80)
    parser.add_argument('--max-step',type=float,default=30.)
    parser.add_argument('--end',type=float,default=1800.)
    parser.add_argument('--D',type=float,default=None)
    parser.add_argument('--relax',action='store_true')
    parser.add_argument('--null-tail',action='store_true')
    parser.add_argument('--output',default='transient.json')
    run(parser.parse_args())
