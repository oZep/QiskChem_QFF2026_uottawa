# Estimate molecular energies with sample-based quantum diagonalization
https://github.com/user-attachments/assets/d8957230-1f38-42c6-a9b2-db33aeda6028

## Common questions (from the challenge brief), answered with our results

Model: LiH, STO-3G, frozen core (1 orbital), 2 electrons in 5 spatial orbitals = 10 qubits, 25 determinants. Energies are total energies in hartree (Ha); errors are in millihartree (mHa) against the exact energy of the same finite Hamiltonian and (1α, 1β) sector.
Numbers below are from the notebook's reference run (seeds 0–4). Re-run the notebook to regenerate them; small differences from library versions are possible.

### 1. Is SQD just VQE with a different name?

**No.** The two methods use the quantum circuit for different jobs.

- **SQD:** the circuit only proposes configurations (which determinants are occupied). A classical solver builds the projected Hamiltonian in that subspace, including off-diagonal elements, and takes its lowest eigenvalue. No parameters are optimised.
- **VQE:** the circuit's expectation value *is* the energy, and a classical optimiser tunes the angles.

Evidence: at 2.4 Å, SQD on the unoptimised CCSD-seeded circuit gives 0.41 mHa error. The same circuit's own ⟨H⟩ is 14 mHa off. SQD only needs the right *configurations*, not the right *amplitudes*. VQE needs about 1,800 energy evaluations to reach exactness (notebook Step 6).

### 2. Must we use hardware?
**No.** We use finite-shot sampling on the Qiskit Aer simulator (noiseless), with fixed seeds, as the brief allows.

- Every SQD sample comes from actually measuring the compiled circuit (`{cx, rz, sx, x}`).
- We never enumerate statevector amplitudes to build the SQD subspace.
- The fast statevector emulator is used only for VQE optimisation and sanity checks.
- The web app has an IBM hardware option that follows the same pipeline. It is **untested on a real device**, so no results in this README come from hardware.

### 3. Do shots equal subspace dimension?
**No.** Three different numbers are reported and they differ:

| Quantity | Meaning | Example (2.4 Å, CCSD-seeded, 300 shots) |
|---|---|---|
| Shots | circuit measurements | 300 |
| Distinct full configurations | unique valid bitstrings | 7 |
| Solve dimension | size of the α×β product space actually diagonalised | 9 |

Reasons: repeated shots hit the same configuration (HF alone carries about 88% of the weight), and the solver works in a product space of unique α strings × β strings. With closed-shell symmetrisation, dimension = (number of spin-strings)², so it can be larger than the count of sampled bitstrings. The full sector is 25 (5² strings).

### 4. Is configuration recovery required?
**No.** Our core scope is clean-sample SQD plus one controlled design comparison (the trial-circuit parameters, notebook Step 7).

- Invalid samples (wrong α/β electron count) are **discarded** and counted. On the noiseless simulator there are none.
- Recovery and noise models are listed as future work. If added, compare against an explicitly recovery-free workflow, not a one-iteration recovery run.

### 5. Should more shots always improve the energy?
**No.** Independent runs at different shot counts sample different subspaces, so they are not nested.

- Only a genuinely nested enlargement of a subspace is guaranteed not to raise the energy.
- Observed at 2.4 Å, CCSD-seeded circuit: 30 shots gave 9.3 ± 20 mHa across seeds. Some seeds found only the HF string (46.96 mHa, the HF error) and some found a useful subspace. At 100 shots the error was 0.41 mHa.
- Errors are bimodal at low shot counts. This is why we report spread over at least 3 seeds (we use 5).

### 6. What if SQD appears better than the exact reference?
**It indicates a bug.** A correctly solved projected Hamiltonian is a variational upper bound on the exact energy of the matching sector. The checks we apply:

- **Consistency:** PySCF FCI, the Jordan–Wigner qubit Hamiltonian restricted to the (1,1) sector, and CASCI all agree to 1e-8 Ha.
- **Solver check:** `solve_fermion` matches an independent hand-built projection of the qubit Hamiltonian to 1e-8 Ha.
- **Offsets:** nuclear repulsion and frozen-core energy enter once, through CASCI's `ecore`.
- **Bound:** every SQD energy is asserted to be ≥ exact (up to numerical tolerance).

If you see SQD below exact, check, in this order: the particle sector, the offset (`ecore` double counted?), the bit ordering, and eigensolver convergence.

### 7. What if the random baseline performs just as well?
**That is a valid result, and here it is more than "just as well".** The baseline samples uniformly from valid (1α, 1β) configurations and uses the identical subspace rules.

- **Equal shots:** the baseline wins. It covers all 5 spin-strings within about 10 shots and returns the exact energy. The CCSD-seeded circuit concentrates on HF and re-samples the same strings. It needs about 100 shots for a 0.41 mHa error at 2.4 Å and does not reach exactness within that budget.
- **Equal solve dimension (subspace capped):** the circuit samples win clearly. At dimension 9 (3 spin-strings): 0.41 mHa vs about 269 mHa for the baseline. At dimension 4: 18.4 mHa vs about 404 mHa.
- **Interpretation:** the circuit is a good *selector* of important configurations but a poor *explorer* of the whole space. In a 25-determinant sector, exploring is cheap and the full solve takes milliseconds, so there is no sampling advantage to claim.
- **Where it could matter:** when the sector is large and the solve dimension is the bottleneck. We checked this directionally on H₂O (12 qubits, 225 determinants): SQD reached 1.04 mHa error at solve dimension 81 (HF error 51 mHa). We have not yet run the baseline comparison there.

**Takeaway:** the most useful workflow in our experiment is the CCSD-seeded circuit + SQD with a capped subspace. It matches VQE-quality energies with no optimisation loop. We do not claim quantum advantage.

To run the code, first create a venv:
```
cd pipeline
pip install -r requirements.txt
python3 app.py
```
