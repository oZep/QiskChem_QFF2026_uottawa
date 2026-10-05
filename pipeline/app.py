"""Molecule -> trial circuit -> measured configurations -> subspace -> diagonalization -> energy (SQD or VQE)."""
import os, time, uuid, threading, traceback
from collections import Counter
from math import comb
import numpy as np
from flask import Flask, jsonify, request, send_from_directory
from pyscf import gto, scf, mcscf, ao2mo, fci, cc
from scipy.optimize import minimize
from scipy.sparse.linalg import expm_multiply
from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_addon_sqd.fermion import solve_fermion

MAX_QUBITS = 14          # simulator limit (statevector). Hardware runs use the same cap.
MHA = 1e3
app = Flask(__name__, static_folder="static")
JOBS = {}


# ---------- 1. chemistry problem ----------
def setup(p):
    mol = gto.M(atom=p["atoms"], basis=p.get("basis", "sto-3g"), charge=int(p.get("charge", 0)),
                spin=0, unit="Angstrom", verbose=0)           # closed-shell only
    mf = scf.RHF(mol).run()
    if not mf.converged: raise ValueError("Hartree-Fock did not converge.")
    nmo, ne = mf.mo_coeff.shape[1], mol.nelectron
    ncore = int(p.get("frozen_core") or 0)
    ncas = int(p.get("n_active") or (nmo - ncore))
    na = (ne - 2 * ncore) // 2
    if ne % 2: raise ValueError("Odd electron count: only closed-shell molecules (spin 0) are supported.")
    if na < 1 or ncas <= na or ncore + ncas > nmo: raise ValueError("Invalid active space for this molecule.")
    if 2 * ncas > MAX_QUBITS: raise ValueError(f"{2*ncas} qubits requested; limit is {MAX_QUBITS}. Freeze more orbitals or shrink the active space.")
    return mol, mf, ncore, ncas, na, nmo

def hamiltonian(mf, ncore, ncas, na):
    mc = mcscf.CASCI(mf, ncas, (na, na)); mc.ncore = ncore
    h1, ecore = mc.get_h1eff()                     # ecore = nuclear repulsion + frozen core (counted once)
    h2 = ao2mo.restore(1, mc.get_h2eff(), ncas)
    e_exact, _ = fci.direct_spin1.FCI().kernel(h1, h2, ncas, (na, na))
    return h1, h2, float(ecore), float(e_exact + ecore)

@app.post("/api/inspect")
def inspect():
    try:
        mol, mf, ncore, ncas, na, nmo = setup(request.json)
        return jsonify(ok=True, electrons=mol.nelectron, orbitals=nmo, qubits=2 * ncas, active_electrons=2 * na,
                       determinants=comb(ncas, na) ** 2, e_hf=float(mf.e_tot))
    except Exception as e:
        return jsonify(ok=False, error=str(e))


# ---------- 2. trial circuit ----------
def ccsd_params(mf, ansatz, ncore, ncas, na, nmo):
    """Classical CCSD amplitudes -> UCCSD angles (no optimisation). Same frozen/active orbitals as the Hamiltonian."""
    frozen = list(range(ncore)) + list(range(ncore + ncas, nmo))
    c = cc.CCSD(mf, frozen=frozen or 0).run()
    x = []
    for occ, virt in ansatz.excitation_list:
        if len(occ) == 1:
            x.append(c.t1[occ[0] % ncas, virt[0] % ncas - na])
        else:
            (i, j), (a, b) = occ, virt
            i, j, a, b = i % ncas, j % ncas, a % ncas - na, b % ncas - na
            same_spin = (occ[0] < ncas) == (occ[1] < ncas)
            x.append(c.t2[i, j, a, b] - c.t2[i, j, b, a] if same_spin else c.t2[i, j, a, b])
    return np.array(x), float(c.e_tot)


# ---------- 3. sampling (simulator or IBM hardware) ----------
def ibm_backend(p, nq):
    from qiskit_ibm_runtime import QiskitRuntimeService
    token = p.get("ibm_token") or os.environ.get("IBM_TOKEN")
    if not token: raise ValueError("IBM hardware selected but no token (paste one or set IBM_TOKEN).")
    kw = dict(channel="ibm_cloud", token=token)
    if p.get("ibm_instance") or os.environ.get("IBM_INSTANCE"): kw["instance"] = p.get("ibm_instance") or os.environ["IBM_INSTANCE"]
    return QiskitRuntimeService(**kw).least_busy(operational=True, simulator=False, min_num_qubits=nq)

def sample(circ, shots, seed, p, hw_backend=None):
    qc = circ.copy(); qc.measure_all()
    if hw_backend is not None:
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
        from qiskit_ibm_runtime import SamplerV2
        isa = generate_preset_pass_manager(backend=hw_backend, optimization_level=3).run(qc)
        counts = SamplerV2(mode=hw_backend).run([isa], shots=shots).result()[0].data.meas.get_counts()
    else:
        isa = transpile(qc, basis_gates=["cx", "rz", "sx", "x"], optimization_level=1, seed_transpiler=0)
        counts = AerSimulator().run(isa, shots=shots, seed_simulator=seed).result().get_counts()
    return counts, dict(depth=isa.depth(), two_qubit_gates=isa.num_nonlocal_gates())


# ---------- 4-5. SQD: configurations -> subspace -> diagonalize ----------
def sqd(counts, h1, h2, ecore, ncas, na, max_strings, include_hf):
    a_cnt, b_cnt, full, valid = Counter(), Counter(), Counter(), 0
    for s, c in counts.items():
        s = s.replace(" ", ""); b, a = s[:ncas], s[ncas:]       # bitstring = [beta][alpha], qubit 0 rightmost
        if a.count("1") != na or b.count("1") != na: continue   # wrong electron count -> discarded
        valid += c; a_cnt[int(a, 2)] += c; b_cnt[int(b, 2)] += c; full[s] += c
    if valid == 0 and not include_hf: raise ValueError("No valid configurations were sampled (all discarded).")
    ranked = [s for s, _ in (a_cnt + b_cnt).most_common()]
    if max_strings: ranked = ranked[:int(max_strings)]
    hf = (1 << na) - 1
    if include_hf and hf not in ranked: ranked.append(hf)
    S = sorted(ranked)                                           # alpha/beta strings (closed-shell symmetrised)
    t0 = time.perf_counter()
    e, *_ = solve_fermion((S, S), h1, h2, open_shell=True)       # full projected H incl. off-diagonals
    return dict(energy=float(e) + ecore, dim=len(S) ** 2, strings=len(S), valid=valid, distinct=len(full),
                solve_s=time.perf_counter() - t0, in_space=set(S),
                top=[dict(bits=s, count=c, kept=(int(s[ncas:], 2) in S and int(s[:ncas], 2) in S))
                     for s, c in full.most_common(12)])


# ---------- job runner ----------
def run_job(jid, p):
    job = JOBS[jid]; t_last = [time.perf_counter()]
    def step(name, summary, **d):
        now = time.perf_counter(); job["steps"].append(dict(name=name, summary=summary, seconds=round(now - t_last[0], 2), **d)); t_last[0] = now
    try:
        method, shots, seed = p.get("method", "sqd"), int(p.get("shots", 1000)), int(p.get("seed", 0))
        mol, mf, ncore, ncas, na, nmo = setup(p)
        h1, h2, ecore, e_exact = hamiltonian(mf, ncore, ncas, na)
        step("Molecule", f"{ncas*2} qubits · {2*na} active electrons · {comb(ncas, na)**2} determinants",
             e_hf=float(mf.e_tot), e_exact=e_exact, frozen_core=ncore, active_orbitals=ncas)

        mapper = JordanWignerMapper()
        ans = UCCSD(ncas, (na, na), mapper, initial_state=HartreeFock(ncas, (na, na), mapper))
        hw = ibm_backend(p, 2 * ncas) if p.get("hardware") else None
        hw_name = hw.name if hw else "Aer simulator (noiseless)"
        res = dict(method=method, e_hf=float(mf.e_tot), e_exact=e_exact, backend=hw_name, qubits=2 * ncas)

        if method == "sqd":
            x, e_ccsd = ccsd_params(mf, ans, ncore, ncas, na, nmo)
            step("Trial circuit", f"HF state + UCCSD ({ans.num_parameters} angles from classical CCSD, no optimisation)",
                 note="CCSD preprocessing is classical cost; the exact ground state is never used.")
            counts, cs = sample(ans.assign_parameters(x), shots, seed, p, hw)
            step("Measured configurations", f"{shots} shots on {hw_name}", shots=shots, **cs)
            r = sqd(counts, h1, h2, ecore, ncas, na, p.get("max_strings"), bool(p.get("include_hf")))
            step("Selected subspace", f"{r['distinct']} distinct configurations → {r['strings']} spin-strings → solve dimension {r['dim']}",
                 invalid_discarded=shots - r["valid"], note="Invalid electron counts are discarded; alpha/beta strings combined into a product space.")
            step("Diagonalization", f"Lowest eigenvalue of the projected Hamiltonian ({r['solve_s']*1e3:.1f} ms)")
            res.update(energy=r["energy"], configs=r["top"], costs=dict(shots=shots, **cs, distinct_configs=r["distinct"],
                       solve_dim=r["dim"], full_sector=comb(ncas, na) ** 2, solve_seconds=round(r["solve_s"], 4)))
        else:
            GEN = [o.to_matrix(sparse=True) for o in ans.operators]
            Hm = mapper.map(ElectronicEnergy.from_raw_integrals(h1, h2).second_q_op())
            hf0 = Statevector(HartreeFock(ncas, (na, na), mapper)).data.astype(complex)
            trace, nfev = [], [0]
            if hw is None:
                M = Hm.to_matrix(sparse=True)
                def f(x):
                    psi = hf0.copy()
                    for k, G in enumerate(GEN): psi = expm_multiply(-1j * x[k] * G, psi)
                    nfev[0] += 1; e = float(np.real(psi.conj() @ (M @ psi))) + ecore; trace.append(e); return e
            else:
                from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
                from qiskit_ibm_runtime import EstimatorV2
                pm = generate_preset_pass_manager(backend=hw, optimization_level=3)
                isa = pm.run(ans); obs = Hm.apply_layout(isa.layout); est = EstimatorV2(mode=hw)
                def f(x):
                    nfev[0] += 1; e = float(est.run([(isa, obs, x)]).result()[0].data.evs) + ecore; trace.append(e); return e
            maxit = int(p.get("max_iter", 30 if hw else 2000))
            step("Trial circuit", f"HF state + UCCSD ({ans.num_parameters} variational angles, start at 0 = Hartree-Fock)")
            r = minimize(f, np.zeros(ans.num_parameters), method="COBYLA", options=dict(maxiter=maxit, rhobeg=0.1, tol=1e-9))
            step("Optimize (VQE loop)", f"{nfev[0]} energy evaluations" + ("" if hw else " (exact statevector expectation, no shot noise)"))
            tq = transpile(ans.assign_parameters(r.x), basis_gates=["cx", "rz", "sx", "x"], optimization_level=1)
            step("Energy", f"Minimum found: {r.fun:.6f} Ha")
            res.update(energy=float(r.fun), trace=trace[::max(1, len(trace)//200)],
                       costs=dict(energy_evaluations=nfev[0], depth=tq.depth(), two_qubit_gates=tq.num_nonlocal_gates(),
                                  parameters=ans.num_parameters, shots="exact (simulated)" if hw is None else "see IBM job usage"))
        res["error_mha"] = (res["energy"] - e_exact) * MHA
        res["hf_error_mha"] = (res["e_hf"] - e_exact) * MHA
        job.update(status="done", result=res)
    except Exception as e:
        traceback.print_exc(); job.update(status="error", error=str(e))

@app.post("/api/run")
def run():
    jid = uuid.uuid4().hex[:8]; JOBS[jid] = dict(status="running", steps=[])
    threading.Thread(target=run_job, args=(jid, request.json), daemon=True).start()
    return jsonify(id=jid)

@app.get("/api/job/<jid>")
def job(jid): return jsonify(JOBS.get(jid, dict(status="error", error="Unknown job")))

@app.get("/")
def index(): return send_from_directory("static", "index.html")

if __name__ == "__main__":
    app.run(debug=False, port=5000)
