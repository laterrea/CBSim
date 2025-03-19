# CBSim - Carnot Battery Simulator

The Carnot Battery Simulator (CBSim) provides thermodynamic models for various Carnot battery architectures. It can be used to carry out energy and exergy analyses of thermally integrated pumped thermal energy storage (TI-PTES).

<p align="center">
  <img src="figs/CBSim_logo.svg" width="600">
</p>

Documentation
=============

The full documentation is not available yet.
* ```main_hp.py```: main script to simulate the heat pumps (vapour compression heat pumps, VCHP);
* ```main_he.py```: main script to simulate the heat engines (organic Rankine cycles, ORC);
* ```main_cb.py```: main script to simulate the Carnot batteries (heat pump and heat engine connected via thermal energy storage); 

The model can represent basic (B) and recuperated (R) cycles. The regime, either subcritical (S) or transcritical (T), can also chosen.
<p align="center">
  <img src="figs/cb_architecture.svg" width="450">
</p>
|                  |                  |                  | \multicolumn{4}{c|}{Heat Pump (HP)} |
|------------------|------------------|------------------|--------------------------------------|
|                  |                  |                  | \multicolumn{2}{c|}{Subcritical (S)} | \multicolumn{2}{c|}{Transcritical (T)} |
|                  |                  |                  | Basic (B) | Recuperated (R) | Basic (B) | Recuperated (R) |
| **Organic Rankine Cycle (ORC)** | **Subcrit. (S)** | Basic (B)  | #1: SBHP + SBORC  | #3: SRHP + SBORC  | #6: TBHP + SBORC  | #11: TRHP + SBORC |
|                  |                  | Recup. (R)  | #2: SBHP + SRORC  | #4: SRHP + SRORC  | #8: TBHP + SRORC  | #13: TRHP + SRORC |
|                  | **Transcrit. (T)** | Basic (B)  | #5: SBHP + TBORC  | #7: SRHP + TBORC  | #9: TBHP + TBORC  | #15: TRHP + TBORC |
|                  |                  | Recup. (R)  | #10: SBHP + TRORC | #12: SRHP + TRORC | #14: TBHP + TRORC | #16: TRHP + TRORC |


References
==========

The model has been used in the following publications:
* R. Tassenoy, A. Laterre, V. Lemort, F. Contino, M. De Paepe, and S. Lecompte, "Assessing the influence of compressor inertia on the dynamic performance of large-scale vapor compression heat pumps for Carnot batteries", *Journal of Energy Storage*, vol. 101, p. 113948, Nov. 2024.<br>
  https://doi.org/10.1016/j.est.2024.113948
* A. Laterre, O. Dumont, V. Lemort and F. Contino, "Extended mapping and systematic optimisation of the Carnot battery trilemma for sub-critical cycles with thermal integration", *Energy*, vol. 304, p. 132006, Sep. 2024.<br>
  https://doi.org/10.1016/j.energy.2024.132006
* A. Laterre, O. Dumont, V. Lemort and F. Contino, "Is waste heat recovery a promising avenue for the Carnot battery? Techno-economic optimisation of an electric booster-assisted Carnot battery integrated into different data centres", *Energy Conversion and Management*, vol. 301, p. 118030, Feb. 2024.<br>
  https://doi.org/10.1016/j.enconman.2023.118030

Cite CBSim
==========

Please use the following reference to cite CBSim when you use it for any publication:
> A. Laterre, O. Dumont, V. Lemort and F. Contino, "Extended mapping and systematic optimisation of the Carnot battery trilemma for sub-critical cycles with thermal integration", *Energy*, vol. 304, p. 132006, Sep. 2024.<br>
  https://doi.org/10.1016/j.energy.2024.132006
