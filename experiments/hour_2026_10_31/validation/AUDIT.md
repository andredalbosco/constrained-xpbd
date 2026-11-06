# Independent review and validation boundaries

Two independent reviews examined the rigid-chain derivations and the work-graph reduction against source code and recorded measurements.

The rigid review verified the quarter/reaction equivalence for linear potential, canonical coordinate scaling, symmetry-aware work correction, stable Hamiltonian differences, and the modal-transfer denominator. It identified 24 failed zero-step records whose RMS value incorrectly appeared as zero; these are now NaN. The tiny-denominator branch is explicitly disclosed as skipping the exact-work correction.

The work review verified elimination signs and factors, the need for component energy closure, the gamma/N normalization dependence, tree forward/transpose operators, weighted graph QPs, and the disconnected-objective experiment. The older AD reducer now rejects disconnected graphs rather than silently keeping only one energy row. Follow-up measurements expose raw energy and reconstructed full balance residuals beside scaled stopping tests. Failure reasons distinguish evaluation-budget exhaustion from unsuccessful line search in the audited nonlinear controls.

The root CMake suite covers analytic/AD derivatives, exact incremental differences, reference refinement and curvature, AL iteration-matrix differentiation, root-preserving numerical coordinates, cached-factor controls, rigid mechanics/derivatives, and work-reduction derivatives. Graph algebra and large tree-operator checks have separate saved CSV/JSON records.

Viewer schema, JavaScript syntax and saved geometric bounds are checked. Static scientific figures were inspected. Interactive browser playback is not verified: supported browser discovery returned no browser, and automatic approval review rejected the isolated Chrome fallback. No fallback browser was launched after that rejection. A separate native MP4 export was verified as 1,920×1,080, 30 fps and 60 seconds, fully decoded, and checked through representative frames; it does not validate the HTML controls.

Physical conclusions require complete trajectories, converged steps and certified references. Some experiment files intentionally retain finite inexact trajectories or failed prefixes as controls; those are explicitly labeled and are not promoted to successful physical integrations.
