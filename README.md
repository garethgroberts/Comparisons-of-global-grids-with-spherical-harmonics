\section{Introduction}
This document describes a minimum working example to generate statistics used to assess similarities and dissimilarities of dynamic topographies and geoids from different global grids. It assumes that we are comparing `observed' dynamic topography, generated using residual ocean age-depth measurements \citep{Hoggard2016,Holdt2022}, or the observed geoid \cite[{\sf EIGEN5c}][]{Foerste2008}, with predictions from global mantle convection simulations \citep[e.g. {\sf TERRA;}][]{Davies2025}. The following is largely adapted from \citet{OMalley2024}. \\

To aid comparison at specific scales, first the two global grids are expanded into spherical harmonic form. Any real, square-integrable function over the surface of the Earth can be described as a function of longitude $\theta$ and latitude $\phi$ by a linear combination of spherical harmonics of degree $l$ and order $m$,
\begin{equation}\label{eq:sph_harm}
    f(\theta,\phi) = \sum_{l=1}^{L}\sum_{m=-l}^{l}f_{lm}Y_{lm}(\theta,\phi),
\end{equation}


\noindent where the spherical harmonic functions $Y_{lm}$ are the natural orthogonal set of basis functions on the sphere, and $f_{lm}$ are the spherical harmonic coefficients \citep[see e.g.][for discussion of spherical harmonics and application to geophysical problems]{Forte2007, Wieczorek2018}. The spherical harmonic coefficients, $f_{lm}$, for irregularly sampled surfaces (e.g. dynamic topography) can be calculated using a least-squares methodology \cite[see e.g.][]{Hoggard2016,Wieczorek2018}. For regularly sampled grids, the expansion can be achieved by by performing Fast Fourier Transforms in longitude and integrating over latitude, see e.g. the {\sf pyshtools} {\tt expand} modules from \citep{Wieczorek2018}. The power at each degree, $l$, is given by 

\begin{equation}
    P_l = \sum^l_{m=-l} f^2_{lm}.
\end{equation}


Several different normalizations are commonly used in the geodynamic community, which can make the comparisons of spherical harmonic coefficients (and derived quantities, e.g. power) less straightforward than might be first thought. \citet{Wieczorek2018} provide a helpful overview of many of the `standard' normalizations used. It is obviously important to ensure the same normalizations are used to generate the grids being compared. Often this requires conversion and consideration of whether the Condon-Shortly phase factor of $(-1)^m$ is included or not. The following statistics are calculated to aid comparison. \\

\subsection{Euclidean Comparisons of Amplitudes}
\label{euclid}

First, we calculate root-mean-squared difference, $\chi$, between predicted surface deflections in the spatial domain, 

\begin{equation}\label{eq:rms}
\chi = \sqrt{\frac{1}{N}\sum_{n=1}^{N}w_{\phi}\left( h^{a}_{n}-h^{b}_{n}\right)^{2}},
\end{equation}

\noindent where $h^{a}_{n}$ and $h^{b}_{n}$ are predicted surface deflections from the two models being compared. $N=$ number of points in the $1\times 1^{\circ}$ gridded maps being compared. The prefactor $w_{\phi}$ is proportional to $\cos{\phi}$, where $\phi$ is latitude, and is included to correct biases in cell size with latitude; mean $w_{\phi} = 1$. This metric is closely associated with the mean vertical distance ($L^{2}$-norm distance) between predicted and reference surfaces, i.e., $\Delta \bar{h} = 1/N\sum_{n=1}^{N}w_{\phi}|h^{a}_{n} - h^{b}_{n}|$. These metrics are sensitive to differences in amplitudes and locations of surface deflections.

\subsection{Spectral Correlation Coefficients}
Secondly, to aid comparisons of surface deflections as a function of scale they are converted into the frequency domain using spherical harmonics.  The degree-correlation spectrum, $r_l$, is calculated using {\sf pyshtools} \citep{Wieczorek2018}, such that 
%
\begin{equation}
    r_l = \frac{S f_1 f_2}{\sqrt{Sf_1 f_1 \cdot Sf_2 f_2}} 
\end{equation}
where  $f_1$ and $f_2$  are the spherical harmonic coefficients of the two estimates of surface deflection being compared. They vary as a function of order, $m$, and degree, $l$; $f=f^m_l$.  $S f_a f_b$ is the cross spectrum of the two functions $f_a$ and $f_b$. We note that $-1\le r_l\le1$, and we calculate the mean value, $\overline{r_l} = 1/L\sum^{L}_{l=1} r_l$, where $L$ is  total number of degrees. Thirdly, the correlation of the entirety of both functions can be estimated following \cite{Forte2015a}, such that
%
\begin{equation} \label{eq:correlation}
    r = \frac{\sum f_{1}^{*} f_{2}}{\sqrt{\sum f_{1}^{*} f_{1}} \sqrt{\sum f_{2}^{*} f_{2}}}, \quad {\rm where} \quad \sum = \sum_{m=-l}^{+l}, 
\end{equation}
%
where $*$ indicates complex conjugation \citep[see also][]{Becker2002, OConnell1971}.  This metric is not sensitive to the amplitudes of surface deflections. \\

\subsection{Comparing Calculated Power Spectra}
Finally, differences in power spectra between between predicted and independent surface deflections are calculated such that 

\begin{equation}\label{eq:chi_p}
\chi_{p} = \sqrt {\frac{1}{L}\sum_{l=1}^{L}\left( \textnormal{log}_{10} P_{l}-\textnormal{log}_{10} P_{l}^{K}\right)^{2}}
+ \sqrt{\frac{1}{L}\sum_{l=1}^{L}\left( \textnormal{log}_{10} P_{l}-\textnormal{log}_{10} P_{l}^{H}\right)^{2}}
,
\end{equation}
where $L$ is the number of spherical harmonic degrees being considered. $P_l^{K}$ and $P_l^{H}$ are total power per degree estimated independently from  Kaula's law or residual oceanic age-depth measurements, respectively \citep[Equation \ref{eq:kaula};][]{Hoggard2016a, Holdt2022}.  Once power spectra are calculated it is straightforward to compare their spectral slopes, which can be used to assess whether broad patterns of surface deflections are similar even if their amplitudes do not match. \\

Using the total power per degree convention, \citet{Hoggard2016}  derived a rule-of-thumb for estimating the power spectrum of dynamic topography (see their Supporting Information), $P_{l}^{K}$, using \citet{Kaula1963}'s approximation for the long-wavelength gravity field of Earth as a function of $l$:

\begin{equation}\label{eq:kaula}
    P_{l}^{K} \approx \left(\frac{GM}{ZR^{2}}\right)^{2} \left(\frac{2}{l}-\frac{3}{l^{2}} + \frac{1}{l^{4}} \right),
\end{equation}

\noindent where $G$ is the gravitational constant, $M=5.97 \times 10^{24}$ kg is the mass of the Earth, $R \approx 6370$~km is Earth's radius. The value of admittance, $Z$, between gravity and topography varies as a function of  viscosity, as well as the depth and wavelength of  internal density anomalies because of the depth- and degree-dependence of their respective sensitivity kernels \citep[see e.g.][and references therein]{Colli2016}.   However, in the upper mantle, which contributes most to surface deflections, the topography and gravity kernels are approximately proportional to one another across all but the lowest spherical harmonic degrees, even when this layer is assumed to be of relatively low viscosity  \citep[see e.g. ][their Figure 2]{Colli2016}. This behavior can explain why \citet{Hoggard2016} found that assuming an average value of $Z = 12$ mGal km$^{-1}$ provides a reasonable approximation of observed residual topographic trends, thus we make use of that value.  Finally, it is useful to note that \citet{Jeans1923} related spherical harmonic degree to wavelength $\lambda$, which at Earth's surface can be approximated via $\lambda \approx 2\pi R / \sqrt{l(l+1)}$.\\

% \subsubsection{Calculating dynamic topography at Earth's surface}
% Surface deformation can be estimated from numeric simulations of mantle convection by making use of the requirement that normal stress is continuous across the upper boundary of the solid Earth \citep[see e.g][]{Parsons1983, McKenzie1977, Ricard2015}.  If the pressure from the overlying column is hydrostatic, the resultant condition is

% \begin{equation} \label{eq:h_balance}
%     \rho_w g_{s} h = \rho_m g_{s} h + \sigma_{rr},
% \end{equation}

% \noindent where $\sigma_{rr}$  incorporates deviatoric viscous stresses generated by mantle convection and dynamic pressure. In practice, since values for this term are obtained by subtracting radial lithostatic stress from the total stress, values of $\sigma_{rr}$ integrate to zero globally. $g_{s}$ is gravitational acceleration at Earth's surface, $\rho_{m}$ is the mean density for the surficial layer, and $\rho_{w}$ is the density of the overlying fluid. Surface deflection arising in response to predicted convective flow, $h$, is approximated by rearranging Equation~\ref{eq:h_balance},
% \begin{equation}  \label{eq:stress_to_dt}
% h \approx - \frac{\sigma_{rr}}{(\rho_{m}-\rho_{w})g_{s}}.
% \end{equation}\\

% These values could then be transformed into the spherical harmonic domain (e.g. solving for the spherical harmonic coefficients in Equation 1). Nonetheless, we make use of an alternative methodology that calculates surface deflection in response to mantle convection using the analytic propagator matrix technique \citep[e.g.,][]{Craig1987, Gantmacher1959, Ghelichkhan2021, Parsons1983, Richards1984}. In this study, following \citet{Ghelichkhan2021} and references therein, surface deflection for each spherical harmonic coefficient, $h_{lm}$, is calculated in the spectral domain such that
 
% \begin{equation} \label{eq:sph_dt}
%         h_{lm} = \frac{1}{(\rho_{m}-\rho_{w})} \int_{R_{\text{CMB}}}^{R} A_{l}\delta \rho_{lm}(r)\cdot\textnormal{d}r.
% \end{equation}

% \noindent Products of the sensitivity kernel, $A_l$, and density anomalies, $\delta \rho_{lm}$, of spherical harmonic degree, $l$, and order, $m$, are integrated with respect to radius, $r$, between the core-mantle boundary and Earth's surface radii, $R_{\textnormal{CMB}}$ and $R$, respectively. The sensitivity kernel is given by

% \begin{equation}\label{eq:surkernel}
% A_{l} = -\left( \frac{\eta_{0}}{Rg_{R}}\right) \left(u_{1}+\frac{\rho_{w}}{\rho_{0}}u_{3} \right),
% \end{equation}
% where $u_{n}(r)$ represents a set of poloidal variables, which are posed for solution of the set of simultaneous equations by matrix manipulation, such that

% \begin{equation}\label{eq:u_matrix}
% u(r) = \begin{bmatrix} y_{1}\eta_{0} & y_{2}\eta_{0}\Lambda & (y_{3} + {\bar{\rho}(r)}y_{5})r & y_{4}r\Lambda & y_{5}r\rho_{0}\Lambda & y_{6}r^{2}\rho_{0} \end{bmatrix}^{T},
% \end{equation}

% \noindent where $\Lambda = \sqrt{l(l+1)}$, and $y_{1}$ to $y_{6}$ represent the spherical harmonic coefficients of radial velocity $v_{r}$, lateral velocity $v_{\theta,\phi}$, radial stress $\sigma_{rr}$, lateral stress $\sigma_{r\theta,\phi}$, gravitational potential $V$, and gravitational potential gradient $\partial V / \partial r$, respectively \citep{Hager1989, Panasyuk1996}. $\bar{\rho}$ is the layer mean ($l=0$) density. The kernel $A_{l}$ includes both $u_{1}$ and $u_{3}$, two terms in the matrix solution to the governing equations that affect surface topography.  They directly exert stress on the surface boundary ($u_{1}$), and change the gravitational potential at the surface ($u_{3}$). The functional forms of calculated sensitivity kernels depend on chosen radial viscosity profiles and boundary conditions  \citep[e.g., free-slip or no-slip;][]{Parsons1983}. 


% \subsection{Spherical harmonic representation of {\sf TERRA} output}
% The spherical harmonic representations of output (density) from {\sf TERRA} that we make use of do not include the Condon-Shortley phase factor. The real-values spherical harmonic functions (`trig' being sin if $m<0$ or otherwise cos) are normalized so that 

% \begin{equation}
%     Y_l^m(
% \theta,\phi) = \overline{P_l^m}(\cos(\theta)) \,\mbox{trig}\,(m \phi)
% \end{equation}

% \noindent satisfying

% \begin{equation}
%     \int_S Y_l^m(\theta,\phi)\cdot Y_l^m(\theta,\phi) = 4 \pi,
% \end{equation}

% \noindent where $(\theta, \phi)$ represent colatitude and longitude. The normalization factor is the preferred one used in geodesy, 
% \begin{equation}
%     \overline{P^{m}_{l}} = \sqrt{2-\delta_{0,m}}\sqrt{(2l+1)\frac{(l-m)!}{(l+m)!}}P^{m}_{l}.
% \end{equation}
% where $\delta_{ij}$ is the Kronecker delta function, and $P^{m}_{l}$ are the associated Legendre polynomials. \\

\section{Implementation}

The associated density files produced from the {\sf TERRA} runs are usually have the prefix {\tt density$\_$sph.037}, where the highest number (e.g. 37) indicates the output from the final time step, i.e. 0 Ma. The associated radial viscosity file must also be downloaded. If the {\sf TERRA} model has non-radial viscosities they must first be converted into a reasonable one-dimensional function (e.g. from the radial averages). Once these two files are available the following {\tt bash} script can be executed\\

{\tt > ./run.sh}\\

\noindent which inserts the spherical harmonic density and radial viscosity files and into {\tt ./GEOID}, which is generated following the instructions at {\tt https://doi.org/10.5281/zenodo.12696774}. Note that the executable {\tt GEOID} is generated and executed when {\tt run.sh} is run. \\ 

{\tt run.sh} assumes that a file listing the values of constants, {\tt const.dat}, exists. That file is expected to contain the values of constants used to produce the {\sf TERRA} run. They are, ordered by column number: 
 \begin{enumerate}
     \item {\tt Output name}: Name of output directory to store results created during this analyses.       
     \item {\tt Density}: Location and name of density file (using ending {\tt *density$\_$sph.0*}).
     \item {\tt Visc}: Location and name of the radial viscosity file. 
     \item {\tt Max degree (DT calc)}: Maximum spherical harmonic degree, $L$, to perform model comparisons, usually 50. 
     \item {\tt Min depth (for DT calc. km)}: Cut-off depth for dynamic topography calculation, densities above this value are excluded from the calculation (usually set to 0, i.e. to include all values up to the surface of the model). 
     \item {\tt Kappa (1=compressible, else incompressible)}:  Confirm whether {\sf TERRA} model was compressible or not.   
     \item {\tt rho$\_$w}:  Density of overlying fluid, usually assumed to be water, with density $\rho_w = 1.05\times10^3$ kg/m$^3$.
     \item {\tt visc$\_$0}: Reference viscosity, $\eta_\circ$, with value determined by what is used in the {\sf TERRA} run, usually $\eta_\circ = 10^{21}$ Pa s.
     \item {\tt l$\_$min}: Minimum degree, $l$,  used to perform model comparisons, usually 1.      
     \item {\tt bgrho}: Background density, now usually a redundant parameter. 
     \item {\tt grav10}: 1 = assume gravity is constant value ($g = 10$ m /s$^{2}$), 0 = self-gravitation. 
     \item {\tt bdry (surface: 1 = free, 0 = no-slip)}: Surface boundary condition, i.e. free- or no-slip. 
     \item {\tt bdry$\_$bot (CMB: 1 = free, 0 = no-slip)}: Core-mantle boundary condition (free- or no-slip). \\
\end{enumerate}

Output to the terminal from running {\tt run.sh}, for the single model in the provided working example, should include the following statistics, which are produced by the {\tt GEOID} code:

\begin{figure*}[h!]
    \centering
    \includegraphics[width=0.5\linewidth]{figures/example.png}
%    \caption{Output to terminal for single comparison run in the provided working example. }
    \label{fig:example_terminal}
\end{figure*}


\noindent Dynamic topography and the geoid (predictions and reference values) and associate statistics can be plotting by running the {\tt bash script}:\\

{\tt > plot$\_$mc2$\_$paper$\_$fulltopo$\_$geoidstats.gmt}. \\

\noindent That and all following plotting scripts work with {\sf GMT} \cite[v6.4;][]{Wessel2019}. Various additional ({\tt bash} and {\sf python} v3.x) scripts that support the plotting are included in the minimum working directory provided, e.g. {\tt my$\_$histo$\_$global$\_$sph$\_$adj.py}, {\tt spectrum.py}, {\tt write$\_$constants$\_$h.sh}, {\tt cross$\_$spectrum.py}, {\tt forte$\_$calc$\_$corr$\_$ggr.py}, some require {\sf pyshtools} to be installed (see {\tt https://shtools.github.io/SHTOOLS/python-installing.html}). This particular plotting script will calculate the statistics again (producing values that should be identical to those produced by {\tt calc$\_$stats$\_$all$\_$mc2$\_$geoid$\_$dt.sh}; see below). It should produce a figure that looks like the following one. The script {\tt cleanup.sh} can be executed to clean up the directory if you wish. \\


\begin{figure}[h!]
    \centering
    \includegraphics[width=\linewidth]{figures/mc2_m_cc_066_u_full_topo_geoid.jpg}
    \caption{Example output from running {\tt > plot$\_$mc2$\_$paper$\_$fulltopo$\_$geoidstats.gmt}. Comparison of modern surface deflections and the geoid predicted by mantle convection model (MCM; {\tt mc2$\_$m$\_$cc$\_$066$\_$u}) with independent observations up to $l=50$ (see body text for details). (a) Water-loaded surface deflections predicted by MCM. (b) Calculated residual topography from \cite{Hoggard2016a}. (c) Solid black = power spectrum of topography shown in panel a. Thin grey curve and band = expected dynamic topography from Kaula's rule using admittance $Z=12\pm3$ mGal km$^{-1}$. Thick grey = power spectra of residual topography shown in panel b. Orange dashed = expected power spectra for water-loaded residual topography from \cite{Holdt2022}. (d) Black/grey = histograms of amplitudes shown in panels a/b. (e) Spectral correlation coefficients, $r_l$, for panels a and b. (f) Black = power spectrum of geoid calculated using {\sf TERRA}. Grey = {\sf Eigen5c} \cite{Foerste2008}. (g) Black/grey = histograms of geoid amplitudes in MCM/{\sf Eigen5c} models. (h) Correlation coefficients for MCM/{\sf Eigen5c}. Note annotated values of $\chi_p$, $\chi$, $\overline{r_l}$ and $r$ are discussed in the body text. }
    \label{fig:example_fig}
\end{figure}

Finally, a benchmark can be produced by running: \\

{\tt > plot$\_$like$\_$for$\_$like$\_$paper$\_$fulltopo$\_$geoidstats.gmt } \\

\noindent which is a simple modification of {\tt plot$\_$mc2$\_$paper$\_$fulltopo$\_$geoidstats.gmt}, calculating and plotting the statistics for reference vs. reference models, as a check that the way statistics are calculated is sensible. \\

Once {\tt GEOID} has been run for each model, statistics for those (singular or many) model comparisons can be produced by running the following script, which also writes the statistics to output files (separate files for dynamic topography and the geoid):\\

{\tt > calc$\_$stats$\_$all$\_$mc2$\_$geoid$\_$dt.sh}.\\

Example output from that script, stored in {\tt dt$\_$scores$\_$mc2$\_$4dp.dat} is \\

\begin{figure*}[h!]
    \centering
    \includegraphics[width=0.75\linewidth]{figures/example_stats.png}
%    \caption{Example output from running {\tt > calc$\_$stats$\_$all$\_$mc2$\_$geoid$\_$dt.sh}. }
    \label{fig:example_stats}
\end{figure*}

\noindent where the column names refer to the model number (see {\tt const.dat}) and the statistics calculated in Sections 1.1--1.3: {\tt CHI$\_$P} = $\chi_p$, {\tt CHI} = $\chi$, {\tt OVERLINE$\_$RL} = $\overline{r_l}$, {\tt R} = $r$. {\tt MC2NAME} stores the name of the model run archived in the {\sf TERRA} simulation log (called e.g. {\tt *MC2$\_$Simulation$\_$Log*}).

\section*{Acknowledgments}
{\sf TERRA} simulations and spherical harmonic representations of model predictions  were produced by J. Panton, A. Plimmer and J. H. Davies. We thank S. Ghelichkhan for his help and for making the code archived at {\tt https://doi.org/10.5281/zenodo.12696774} available to us. We also thank F. Richards and V. Fernandes for their help. 

\bibliography{bibtexpaper_ggr}


\end{document}
