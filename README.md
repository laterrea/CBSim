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
\[
	\begin{array}{ccc|cccc|}
		\cline{4-7}
		&                                                                                                &                                                      & \multicolumn{4}{c|}{Heat Pump (HP)}                                                                                                                                                                                                                                                                                      \\ \cline{4-7} 
		&                                                                                                &                                                      & \multicolumn{2}{c|}{\begin{tabular}[c]{@{}c@{}}Subcritical\\ (S)\end{tabular}}                                                                                        & \multicolumn{2}{c|}{\begin{tabular}[c]{@{}c@{}}Transcritical\\ (T)\end{tabular}}                                                                 \\ \cline{4-7} 
		&                                                                                                &                                                      & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}Basic\\ (B)\end{tabular}}          & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}Recuperated\\ (R)\end{tabular}}    & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}Basic\\ (B)\end{tabular}}          & \begin{tabular}[c]{@{}c@{}}Recuperated\\ (R)\end{tabular}    \\ \hline
		\multicolumn{1}{|c|}{\multirow{4}{*}{\begin{tabular}[c]{@{}c@{}}Organic\\ Rankine\\ Cycle\\ (ORC)\end{tabular}}} & \multicolumn{1}{c|}{\multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}Subcrit.\\ (S)\end{tabular}}}   & \begin{tabular}[c]{@{}c@{}}Basic\\ (B)\end{tabular}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#1: SBHP\\ + SBORC\end{tabular}}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#3: SRHP\\ + SBORC\end{tabular}}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#6: TBHP\\ + SBORC\end{tabular}}  & \begin{tabular}[c]{@{}c@{}}\#11: TRHP\\ + SBORC\end{tabular} \\ \cline{3-7} 
		\multicolumn{1}{|c|}{}                                                                                           & \multicolumn{1}{c|}{}                                                                          & \begin{tabular}[c]{@{}c@{}}Recup.\\ (R)\end{tabular} & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#2: SBHP\\ + SRORC\end{tabular}}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#4: SRHP\\ + SRORC\end{tabular}}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#8: TBHP\\ + SRORC\end{tabular}}  & \begin{tabular}[c]{@{}c@{}}\#13: TRHP\\ + SRORC\end{tabular} \\ \cline{2-7} 
		\multicolumn{1}{|c|}{}                                                                                           & \multicolumn{1}{c|}{\multirow{2}{*}{\begin{tabular}[c]{@{}c@{}}Transcrit.\\ (T)\end{tabular}}} & \begin{tabular}[c]{@{}c@{}}Basic\\ (B)\end{tabular}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#5: SBHP\\ + TBORC\end{tabular}}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#7: SRHP\\ + TBORC\end{tabular}}  & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#9: TBHP\\ + TBORC\end{tabular}}  & \begin{tabular}[c]{@{}c@{}}\#15: TRHP\\ + TBORC\end{tabular} \\ \cline{3-7} 
		\multicolumn{1}{|c|}{}                                                                                           & \multicolumn{1}{c|}{}                                                                          & \begin{tabular}[c]{@{}c@{}}Recup.\\ (R)\end{tabular} & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#10: SBHP\\ + TRORC\end{tabular}} & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#12: SRHP\\ + TRORC\end{tabular}} & \multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}}\#14: TBHP\\ + TRORC\end{tabular}} & \begin{tabular}[c]{@{}c@{}}\#16: TRHP\\ + TRORC\end{tabular} \\ \hline
	\end{array}
\]

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
