# Transverse-Cavity-Gross-Pitaevskii-Equation
 a Python repository for simulations for the Bloch Oscillations in quasi-3D Gross-Pitaevskii Equation for cylindrically symmetric clouds coupled to a Cavity. 

 This repository aims to 



 From equation (\ref{3D GPE}) and adding contribution from Harmonic Trap, we have:
$$i\hbar\frac{\partial {\Psi(\bm{r}, t)}}{\partial t} = \left[\frac{-\hbar^2}{2m}\nabla^2 + Fz + \frac{g^2_0|\alpha|^2}{\Delta_a}e^{-\frac{r^2}{w^2}}\cos^2{kz} + \frac{1}{2}m\omega^2_\perp(x^2 + y^2) + \frac{1}{2}m\omega^2z^2 + NU|\Psi|^2 \right]\Psi(\bm{r}, t)$$

Above equation gives us the opportunity to exploit cylindrical symmetry of the system by using cylindrical coordinates $(r, \theta, z)$.

$$\Psi(\bm{r}, t) = \exp(i m\theta)\psi(r, z, t)$$
But for ground state, $m = 0 \Rightarrow e^{im\theta} = 1 $


$$i\hbar\frac{\partial {\psi(r, z, t)}}{\partial t} = \left[\frac{-\hbar^2}{2m}\nabla_{rz}^2 + Fz + \frac{g^2_0|\alpha|^2}{\Delta_a}e^{-\frac{r^2}{w^2}}\cos^2{kz} + \frac{1}{2}m\omega^2_\perp r^2 + \frac{1}{2}m\omega^2z^2 + NU|\psi|^2 \right]\psi(r, z, t)$$

So, methods given in \cite{3952136, PhysRevA.74.013623} suggests how to deal with this numerically.

$$\psi(r, z, t) \rightarrow \boxed{\times \hat{U}_r(\frac{dt}{2})} \rightarrow FFT_z \  and \ HT_r \rightarrow \tilde{\psi}(k_r, k_z, t) \rightarrow \boxed{\times \hat{U}_k(dt)} \rightarrow$$ $$\rightarrow iFFT_z \ and \ iHT_r \rightarrow \boxed{\times \hat{U}_r(\frac{dt}{2})} \rightarrow \psi(r, z, t+dt)$$

Where, $(i)FFT_z$ and $(i)HT_r$ denotes (inverse) Fast Fourier Transform along $z$ and (inverse) Hankel Transform along $r$.

\textit{Note: When referring to Hankel Transform, we mean Zeroth-Order Hankel Transform only.}

\textbf{For Normalization, we use:}
$$ N = 2\pi \int_{-\infty}^{\infty}\int_0^\infty|\Psi(r, z)|^2 r dr dz$$ (In Cylindrical Coordinates)

\textbf{As is with Uncertainties:} $$\Delta z = \sqrt{\langle z^2\rangle - \langle z\rangle^2} $$

$$\langle z \rangle = \int_{-\infty}^{\infty} \int_0^{\infty} z|\Psi(r,z )|^2 rdrdz$$
$$\langle z^2 \rangle = \int_{-\infty}^{\infty} \int_0^{\infty} z^2|\Psi(r,z )|^2 rdrdz$$

Likewise,
$$\Delta r = \sqrt{\langle r^2\rangle - \langle r\rangle^2} $$
$$\langle r \rangle = \int_{-\infty}^{\infty} \int_0^{\infty} |\Psi(r,z )|^2 r^2drdz$$
$$\langle r^2 \rangle = \int_{-\infty}^{\infty} \int_0^{\infty}|\Psi(r,z )|^2 r^3drdz$$

And with, $$\Delta k_z = \sqrt{\langle k_z^2\rangle - \langle k_z\rangle^2} $$

$$\langle k_z \rangle = \int_{-\infty}^{\infty} \int_0^{\infty} k_z|\Phi(k_r,k_z )|^2 k_rdk_rdk_z$$
$$\langle k_z^2 \rangle = \int_{-\infty}^{\infty} \int_0^{\infty} k_z^2|\Phi(k_r,k_z )|^2 k_rdk_rdk_z$$

Likewise,
$$\Delta k_r = \sqrt{\langle k_r^2\rangle - \langle k_r\rangle^2} $$

$$\langle k_r \rangle = \int_{-\infty}^{\infty} \int_0^{\infty} |\Phi(k_r,k_z )|^2 k_r^2dk_rdk_z$$
$$\langle k_r^2 \rangle = \int_{-\infty}^{\infty} \int_0^{\infty}|\Phi(k_r,k_z )|^2 k_r^3dk_rdk_z$$

\textbf{Separate probability density for r and z} \\
\textbf{for z,}
$$|\phi(z)|^2 = 2\pi\int_0^{\infty} |\Psi(r,z )|^2 rdr$$
\textbf{for r,}
$$|\varphi(r)|^2 = \int_{-\infty}^{\infty} |\Psi(r,z )|^2 dz$$

\textbf{Sliced 3D-wavefunction}
$$|\psi(z)|^2 = |\Psi(r = 0, z)|^2$$
$$|\chi(r)|^2 = |\Psi(r, z = 0)|^2$$


\subsubsection{Benchmarks:}
\subsubsection{Schrodinger Equation: Cylindrical}
\subsubsubsection{Harmonic Oscillator}
\textbf{A. Harmonic Oscillator} \\
3D Isotropic Harmonic Oscillator in Cylindrical Coordinates

Isotropic:: $\omega_r = \omega_z = \omega$

$$i\hbar\frac{\partial \Psi(r,z, t)}{\partial t} = -\frac{\hbar^2}{2m}\nabla^2\Psi + \frac{1}{2}m\omega^2r^2\Psi + \frac{1}{2}m\omega^2z^2\Psi$$

with, $$E = \hbar \omega_z + \frac{1}{2}\hbar \omega_r$$

Ground state energy: $$E = \frac{3}{2}\hbar\omega $$
with:
$$\Psi(r,z ) = Ne^{-\frac{\beta}{2}(r^2 + z^2)}$$
$\beta = mw/\hbar$ \\


Extracting scale:
$$i\hbar w\frac{\partial \Psi(r,z, t)}{\partial \tilde{t}} = -\frac{\hbar^2}{2m}\nabla^2\Psi + \frac{1}{2}m\omega_r^2r^2\Psi + \frac{1}{2}m\omega_z^2z^2\Psi$$


non-dimensionalzing:
$$i\frac{\partial \Psi}{\partial \tilde{t}} =\left[ \frac{\tilde{k}_r^2}{2} + \frac{\tilde{k}_z^2}{2} + \alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} \right]\Psi$$

Where, $$\alpha_{r, z} = \frac{\omega_{r,z}^2}{w^2}$$ 
$w = \sqrt{\omega_r \omega_z}$ (can choose as $\omega_{r,z}$ too.)\\

with 
$$\Psi(\tilde{r},\tilde{z} ) = Ne^{-\frac{1}{2}(\sqrt{\alpha_r}\tilde{r}^2 + \sqrt{\alpha_z} \tilde{z}^2)}$$
$$\frac{E}{\hbar w} = \sqrt{\alpha_r} + \frac{1}{2}\sqrt{\alpha_z}$$
Ground state in this unit, when $$\alpha_r = \alpha_z = 1$$ $\Rightarrow$ $$\tilde{E} = 3/2$$

\textbf{A.1 Ground state} \\
\\
%\begin{center}
\begin{tabular}{ |c|c|c|c|c| } 
 \hline
 $\alpha_r$ & $\alpha_z$ & $\tilde{E}_{exact}$ & $\tilde{E}_{num}$ & $\% Error$ \\ 
 \hline
 1.0 & 1.0 & 1.5 & 1.499370537110958 & 0.042  \\ 
 1.0 & 25.0 & 3.5 & 3.4742432518202526 & 0.735\\ 
 9.0 & 0.01 & 3.05 & 3.0388126321508295 & 0.366  \\ 
 \hline
\end{tabular}
%\end{center}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha 1.png}
    \includegraphics[width = 7cm]{figures/HO exact gnd state alpha 1.png}
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha 1 z slice.png}
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha 1 r slice.png}
    %\includegraphics[width = 7cm]{Images/ananum E - P.png}
    \caption{Plots for Probability densities for Isotropic Harmonic Trap with $\alpha_r = 1.0$ and $\alpha_z = 1.0$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha z 25.png}
    \includegraphics[width = 7cm]{figures/HO exact gnd state alpha z 25.png}
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha z 25 z slice.png}
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha z 25 r slice.png}
    %\includegraphics[width = 7cm]{Images/ananum E - P.png}
    \caption{Plots for Probability densities for Anisotropic Harmonic Trap with $\alpha_r = 1.0$ and $\alpha_z = 25.0$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha r 9.png}
    \includegraphics[width = 7cm]{figures/HO exact gnd state alpha r 9.png}
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha r 9 z slice.png}
    \includegraphics[width = 7cm]{figures/HO ITE gnd state alpha r 9 r slice.png}
    %\includegraphics[width = 7cm]{Images/ananum E - P.png}
    \caption{Plots for Probability densities for Anisotropic Harmonic Trap with $\alpha_r = 9.0$ and $\alpha_z = 0.01$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\newpage
\subsubsubsection{Bloch Oscillations}
\textbf{B. Bloch Oscillations}

$$i\hbar\frac{\partial \Psi(r,z, t)}{\partial t} = -\frac{\hbar^2}{2m}\nabla^2\Psi + \frac{1}{2}m\omega^2r^2\Psi + \frac{1}{2}m\omega^2z^2\Psi + Fz\Psi + V_pcos^2(kz)\Psi$$

\textit{Note: Capital '\textbf{R}' for Recoil and small '\textbf{r}' for r-axis.}

taking $$kz = \tilde{z}, E_R = \frac{\hbar^2 k^2}{m}, \omega_R = \frac{\hbar k^2 }{m}, \omega_B = \frac{Fd}{\hbar}, k = \frac{\pi}{d}$$

\textit{Note: Here, $E_R$ is actually twice of the recoil energy.}

non-dim::
$$i\frac{\partial \Psi(r,z, t)}{\partial t} = \left[ \frac{\tilde{k}_r^2}{2} + \frac{\tilde{k}_z^2}{2} + \alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2}  + \tilde{F}\tilde{z} + V_pcos^2(\tilde{z})\right]\Psi$$

$$\alpha_{r, z} = \frac{m^2 \omega^2_{r, z}}{\hbar^2 k^4} = \frac{\omega_{r,z}^2}{\omega_R^2}; \tilde{F} = \frac{Fm}{\hbar^2 k^3} = \frac{\omega_B}{\pi\omega_R}; \tilde{V}_p = \frac{V_p}{E_R}$$

Also, Bloch period $(T_B)$ is scaled as $\omega_R T_B$. So, $$\tilde{T}_B = \frac{2}{|\tilde{F}|}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 0b.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 05 b.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 1b.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for Probability densities for one full oscillation. The Bloch period $(\tilde{T}_B)$ is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 0b z slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 05 b z slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 1b z slice.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for z-slice Probability densities for one full oscillation. The Bloch period $(\tilde{T}_B)$ is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
\newpage
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 0b kz slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 05 b kz slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO gnd t 1b kz slice.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for (momentum ) kz-slice of k-space Probability densities for one full oscillation. The Bloch period $(\tilde{T}_B)$ is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
    
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO exp z.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO exp z sq.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO exp kz.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO exp kz sq.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO exp r.png}
    \includegraphics[width = 7cm]{figures/Schro BO/SE BO exp r sq.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for Expectation values $\langle z, z^2, k_z, k_z^2, r, r^2 \rangle$. Each dotted vertical line represent a Bloch period.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
% \newpage
% blank
\newpage

\textbf{Case:} When Radial Trap frequency ($\omega_r$) is very high.
Changing this frequency does not affect Bloch Oscillations at all.
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 0b.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 05 b.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 1 b.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for Probability densities for one full oscillation when $\omega_r$ is very high. The Bloch period $(\tilde{T}_B)$ is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 0b z slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 05 b z slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 1 b z slice.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for z-slice Probability densities for one full oscillation when $\omega_r$ is very high. The Bloch period $(\tilde{T}_B)$ is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 0b kz slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 05 b kz slice.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO gnd t 1 b kz slice.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for (momentum ) kz-slice of k-space Probability densities for one full oscillation when $\omega_r$ is very high. The Bloch period $(\tilde{T}_B)$ is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
    
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO exp z.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO exp z sq.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO exp kz.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO exp kz sq.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO exp r.png}
    \includegraphics[width = 7cm]{figures/Schro BO high omega r/SE high omegr BO exp r sq.png}
    \caption{Schrodinger's Equation Bloch Oscillations: Plots for Expectation values $\langle z, z^2, k_z, k_z^2, r, r^2 \rangle$. Each dotted vertical line represent a Bloch period.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}


\newpage
blank
\newpage

\subsubsection{Gross-Pitaevskii Equation}
\textbf{A. Harmonic Oscillator}

\begin{equation}\label{GPE HO g low}
    i\frac{\partial \Psi(r,z, t)}{\partial t} = \left[ \frac{\tilde{k}_r^2}{2} + \frac{\tilde{k}_z^2}{2} + \alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} +  \tilde{g}|\Psi|^2 \right]\Psi
\end{equation}


$$\alpha_{r, z} = \frac{\omega_{r,z}^2}{w^2}; \tilde{g} = 4\pi\frac{a_S}{L_{HO}}; L_{HO} = \sqrt{\frac{\hbar}{mw}} $$

\textbf{A.1 For Low values of $\tilde{g}$} \\
In the limit of when $\tilde{g}$ is very low, In equation (\ref{GPE HO g low}) essentially Kinetic Energy and Potential Energy dominates and thus our wavefunction looks similar to as that of Schrodinger Equation's Harmonic Trap (SEHO) wavefunction. However, because of inter-particle interactions, nor wavefunction neither the energy are equal to SEHO's. However, we can still calculate the wavefunction and energy using a variational calculation.

In the absence of interparticle interactions the lowest single-particle state has the familiar wave function:
$$\Psi(r, z) = N_c e^{-\alpha r^2} e^{-\beta z^2}$$

And to calculate the Ground state $\Psi(r,z)$ and the ground state energy for eq () with inter-particle interaction, we use Variational Principle and use the $\Psi(r,z)$ of eq () as Trial wavefunction and vary the parameters $\alpha$ and $\beta$.

$$E_{True \ gnd } \leq \min\frac{\bra{\Psi_T(\alpha, \beta)}H\ket{\Psi_T(\alpha, \beta)}}{\langle\Psi_T(\alpha, \beta)\ket{\Psi_T(\alpha, \beta)}}$$

\begin{figure}[h]
    \centering
    \includegraphics[width = 7cm]{figures/Variational low g/GPE variational gnd state low g.png}
    \includegraphics[width = 7cm]{figures/Variational low g/GPE ITE gnd state low g.png}
    \includegraphics[width = 7cm]{figures/Variational low g/GPE ITE var low g z slice.png}
    \includegraphics[width = 7cm]{figures/Variational low g/GPE ITE var low g r slice.png}
    %\includegraphics[width = 7cm]{Images/ananum E - P.png}
    \caption{Comparison of Probability densities obtained using Variational Method and Imaginary Time Method for low value of $\tilde{g}$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}


Also, We compare this energy with Analytical energy at different value of $\tilde{g}$  

$$E_{Analytical} = E_{har} + E_{int}$$
$$\frac{E_{har}}{\hbar w} = \sqrt{\alpha_r} + \frac{1}{2}\sqrt{\alpha_z}$$
$$\frac{E_{int}}{\hbar w} = 2\pi \tilde{g}\int |\Psi_T|^4 rdrdz$$

However, this $E_{Analytical}$ is valid only for small value of $\tilde{g}$ (compared to HO energy).

\begin{tabular}{ |c|c|c|c|c|c| } 
 \hline
 $\tilde{g}$ & $\alpha$ & $\beta$ & $\tilde{E}_{Numerical}$ & $\tilde{E}_{Analytical}$ & $\% Error$ \\ 
 \hline
 $10^{-9}$ & 0.50084685 & 0.49999993 & 1.4999938660913728 & 1.50000000195879 & 0.0004090578277  \\ 
$10^{-8}$ & 0.50084685 & 0.49999994 & 1.4999938692680177 & 1.50000001958797 & 0.00041002132485\\ 
$10^{-7}$ & 0.50084685 & 0.49999994 & 1.4999939010344503 & 1.500000195879705 & 0.00041965629556 \\ 
$10^{-6}$ & 0.50084687 & 0.49999995 & 1.4999942186987696 & 1.50000195879694 & 0.000516005871392 \\
$10^{-5}$ & 0.50084686 & 0.49999992 & 1.4999973953420307 & 1.50001958797083 & 0.0014794892670 \\
$10^{-4}$ & 0.50084689 & 0.49999993 & 1.5000291617746084 & 1.500195879691797 & 0.0111130765952 \\
$10^{-3}$ & 0.50084715 & 0.49999995 & 1.500346826100257 & 1.50195879585799 & 0.10732449932575 \\
$10^{-2}$ & 0.50084985 & 0.49999989 & 1.503523469340814 & 1.51958785408322 & 1.057154063139 \\
$10^{-1}$ & 0.5008768 & 0.49999993 & 1.5352899001593951 & 1.69586795448250 & 9.468782867125 \\
 \hline
\end{tabular}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Variational low g/GPE comparison var low g sec.png}
    \caption{Ground state energy vs. $\tilde{g}$ for Analytical and Numerical variational Ground state energy}
    %\label{fig:Phaseplot of MBwD}
\end{figure}


\subsubsubsection{Dynamics - When Harmonic Trap is shifted to some offset for low $\tilde{g}$}
$$\alpha_z\frac{\tilde{z}^2}{2} \rightarrow \alpha_z\frac{(\tilde{z} - \tilde{z}_T)^2}{2}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/Psi gnd zt 5 g 1e1 t 0.png}
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/Psi gnd zt 5 g 1e1 t 150.png}
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/Psi gnd zt 5 g 1e1 t 300.png}
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/Psi gnd zt 5 g 1e1 t 450.png}
    
    \caption{Dynamics of GPE - When $\tilde{g}$ is low and Trap is shifted to some $\tilde{z}_T$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/TF gnd zt 5 g 1e1 exp z.png}
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/TF gnd zt 5 g 1e1 exp z sq.png}
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/TF gnd zt 5 g 1e1 exp r.png}
    \includegraphics[width = 7cm]{figures/Low g shifted trap dyn/TF gnd zt 5 g 1e1 exp r sq.png}
    \caption{Dynamics of Expectation value - When Trap is shifted to $\tilde{z_T}$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
\newpage
\\
\subsubsubsection{Thomas-Fermi}
\textbf{A.2 For High values of $\tilde{g}$} \\
In the limit of when $\tilde{g}$ is very high, we use the Thomas-Fermi approximation. Which means, We can neglect the kinetic energy from GPE and thus we have:

$$\mu\Psi(r,z) = \left[ \frac{\tilde{k}_r^2}{2} + \frac{\tilde{k}_z^2}{2} + \alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} +  \tilde{g}|\Psi|^2 \right]\Psi$$

In limit $\tilde{g} \sim $ very high, We have

$$\mu\Psi(r,z) \simeq \left[\alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} +  \tilde{g}|\Psi|^2 \right]\Psi$$
$\Rightarrow$
$$\left[ \mu - \alpha_r\frac{\tilde{r}^2}{2} - \alpha_z\frac{\tilde{z}^2}{2} -  \tilde{g}n(\bm{r}) \right]\Psi = 0$$
$\Rightarrow$
$$\tilde{g}n(\bm{\tilde{r}}) =  \mu - \alpha_r\frac{\tilde{r}^2}{2} - \alpha_z\frac{\tilde{z}^2}{2}$$
$\Rightarrow$
\begin{equation}\label{Num density TF}
    n_{TF}(\bm{\tilde{r}}) = \frac{\mu}{\tilde{g}}\left[ 1 - \frac{\tilde{r}^2}{R^2_{cut}} - \frac{\tilde{z}^2}{Z^2_{cut}} \right]
\end{equation}
Where, $$R^2_{cut} = \frac{2\mu}{\alpha_r}; Z^2_{cut} = \frac{2\mu}{\alpha_z}$$

Equation (\ref{Num density TF}) is Number density in Thomas-Fermi approximation.
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 6.5cm]{figures/Thomas Fermi gnd state/TF ITE g 1300.png}
    \includegraphics[width = 6.5cm]{figures/Thomas Fermi gnd state/TF ITE g 1300 z-slice.png}
    \includegraphics[width = 6.5cm]{figures/Thomas Fermi gnd state/TF ITE g 1300 r-slice.png}
    %\includegraphics[width = 7cm]{Images/ananum E - P.png}
    \caption{Thomas-Fermi Ground State for $\tilde{g} = 1300$. Numerical Peak of Thomas-Fermi: 0.0073224023426329355 and Analytical Peak of Thomas-Fermi: 0.00738371712287171. $\% Error:$  0.8304053259143223}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\subsubsubsection{Thomas-Fermi Dynamics}
$$\alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} \rightarrow \alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2}$$
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF gnd alpha rz eq rz t 0.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF gnd alpha rz eq rz t 20.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF gnd alpha rz eq rz t 40.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF gnd alpha rz eq rz t 80.png}
    \caption{Dynamics of Thomas-Fermi ground state - When no changes are made to Trap}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF dyn alpha rz eq rz exp z.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF dyn alpha rz eq rz exp z sq.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF dyn alpha rz eq rz exp r.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha rz eq rz/TF dyn alpha rz eq rz exp r sq.png}
    \caption{Expectation value of TF ground state}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
\newpage
$$\alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} \rightarrow \lambda\alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF gnd alpha r eq lam r t 0.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF gnd alpha r eq lam r t 20.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF gnd alpha r eq lam r t 40.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF gnd alpha r eq lam r t 60.png}
    %\includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF gnd alpha r eq lam r t 80.png}
    \caption{Breathing dynamics of TF ground state - when the Radial Trap is stretched}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF dyn alpha r eq lam r exp z.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF dyn alpha r eq lam r exp z sq.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF dyn alpha r eq lam r exp r.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r eq lam r/TF dyn alpha r eq lam r exp r sq.png}
    \caption{Expectation value of TF ground state - when the Radial Trap is stretched}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

$$\alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} \rightarrow \alpha_r\frac{\tilde{r}^2}{2} + \lambda \alpha_z\frac{\tilde{z}^2}{2}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF gnd alpha z eq lam z t 0.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF gnd alpha z eq lam z t 10.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF gnd alpha z eq lam z t 30.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF gnd alpha z eq lam z t 50.png}
    %\includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF gnd alpha z eq lam z t 70.png}
    \caption{Breathing dynamics of TF ground state - when the Axial Trap is stretched}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF dyn alpha z eq lam z exp z.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF dyn alpha z eq lam z exp z sq.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF dyn alpha z eq lam z exp r.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha z eq lam z/TF dyn alpha z eq lam z exp r sq.png}
    \caption{Expectation value of TF ground state - when the Axial Trap is stretched}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

$$\alpha_r\frac{\tilde{r}^2}{2} + \alpha_z\frac{\tilde{z}^2}{2} \rightarrow \lambda\alpha_r\frac{\tilde{r}^2}{2} + \lambda \alpha_z\frac{\tilde{z}^2}{2}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF gnd alpha r z eq lam r z t 0.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF gnd alpha r z eq lam r z t 20.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF gnd alpha r z eq lam r z t 40.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF gnd alpha r z eq lam r z t 60.png}
    %\includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF gnd alpha r z eq lam r z t 80.png}
    \caption{Breathing dynamics of TF ground state - when the Axial and Radial Trap is stretched}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF dyn alpha rz eq lam rz exp z.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF dyn alpha rz eq lam rz exp z sq.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF dyn alpha rz eq lam rz exp r.png}
    \includegraphics[width = 7cm]{figures/TF breath alpha r z eq lam r z/TF dyn alpha rz eq lam rz exp r sq.png}
    \caption{Expectation value of TF ground state - when the Axial and Radial Trap is stretched}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\newpage

$$\alpha_z\frac{\tilde{z}^2}{2} \rightarrow \alpha_z\frac{(\tilde{z} - \tilde{z}_T)^2}{2}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 t 0.png}
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 t 150.png}
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 t 300.png}
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 t 450.png}
    \caption{Dynamics of TF - Trap is shifted to some $\tilde{z}_T$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 exp z.png}
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 exp z sq.png}
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 exp r.png}
    \includegraphics[width = 7cm]{figures/TF shifted trap dynamics/TF gnd zt 5 g 10 exp r sq.png}
    \caption{Expectation values when Trap is shifted to some $\tilde{z}_T$}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\newpage
blank
\newpage
\subsubsubsection{Bloch Oscillations}
\textbf{B. Bloch Oscillations}
\begin{equation}
i\hbar\frac{\partial {\psi(r, z, t)}}{\partial t} = \left[\frac{-\hbar^2}{2m}\nabla_{rz}^2 + Fz + V_p\cos^2{kz} + \frac{1}{2}m\omega^2_\perp r^2 + \frac{1}{2}m\omega^2z^2 + NU|\psi|^2 \right]\psi(r, z, t)
\end{equation}

$$kz = \tilde{z}, E_R = \frac{\hbar^2 k^2}{m}, \omega_R = \frac{\hbar k^2 }{m}, \omega_B = \frac{Fd}{\hbar}, k = \frac{\pi}{d}$$

Non-dim::
\begin{equation}
i\frac{\partial \Psi}{\partial t} = \left[-\frac{1}{2}(\tilde{k}_r^2 + \tilde{k}_z^2) + \frac{1}{2}\alpha_r \tilde{r}^2 + \frac{1}{2}\alpha_z \tilde{z}^2 + \tilde{F}\tilde{z} + \tilde{V}_p\cos^2{(\tilde{z})}  + \tilde{g}|\Psi|^2  \right]\Psi
\end{equation}

$$\alpha_{r, z} = \frac{m^2 \omega^2_{r, z}}{\hbar^2 k^4} = \frac{\omega_{r,z}^2}{\omega_R^2}; \tilde{F} = \frac{Fm}{\hbar^2 k^3} = \frac{\omega_B}{\pi\omega_R}; \tilde{V}_p = \frac{V_p}{E_R}; \tilde{g} = \frac{U k^3}{E_R} = 4\pi k a_s  $$

Bloch period $(T_B)$ is scaled as $\omega_R T_B$. So, $$\tilde{T}_B = \frac{2}{|\tilde{F}|}$$

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 0b g 1e1.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 05 b g 1e1.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 1 b g 1e1.png}
    \caption{GPE Bloch Oscillations: Plots for Probability densities for one full oscillation. The Bloch period ($T_B$) is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}

\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 0b g 1e1 z slice.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 05 b g 1e1 z slice.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 1 b g 1e1 z slice.png}
    \caption{GPE Bloch Oscillations: Plots for z-slice Probability densities for one full oscillation. The Bloch period ($T_B$) is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
\newpage
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 0b g 1e1 kz slice.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 05 b g 1e1 kz slice.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO gnd t 1 b g 1e1 kz slice.png}
    \caption{GPE Bloch Oscillations: Plots for (momentum) kz-slice k-space Probability densities for one full oscillation. The Bloch period ($T_B$) is 400.}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
    
\begin{figure}[!htb]
    \centering
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO g 1e1 exp z.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO g 1e1 exp z sq.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO g 1e1 exp kz.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO g 1e1 exp kz sq.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO g 1e1 exp r.png}
    \includegraphics[width = 7cm]{figures/GPE BO g 0.1/GPE BO g 1e1 exp r sq.png}

    \caption{GPE Bloch Oscillations - Plots for Expectation values $\langle z, z^2, k_z, k^2_z, r, r^2 \rangle$. Each dotted vertical line represents a Bloch period}
    %\label{fig:Phaseplot of MBwD}
\end{figure}
\newpage
\subsubsection{Cavity-GPE with Transverse Coupling }

\begin{equation}
i\hbar\frac{\partial {\psi(r, z, t)}}{\partial t} = \left[\frac{-\hbar^2}{2m}\nabla_{rz}^2 + Fz + \frac{g^2_0|\alpha|^2}{\Delta_a}e^{-\frac{r^2}{w^2}}\cos^2{kz} + \frac{1}{2}m\omega^2_r r^2 + \frac{1}{2}m\omega_z^2z^2 + NU|\psi|^2 \right]\psi(r, z, t)
\end{equation}

\begin{equation}
\frac{\partial \alpha}{\partial t} = i\left[\Delta_c - \frac{N}{\Delta_a}\int 2\pi r dr dz g^2_0e^{-\frac{r^2}{w^2}}\cos^2{kz}|\psi(r,z, t)|^2\right]\alpha + \eta - \kappa \alpha
\end{equation}

% *Non-Dimensionalized $\Rightarrow$ ($r$ and $z$ are non-dimensionalized)
Non-dim:
\begin{equation}
i\frac{\partial \Psi}{\partial t} = \left[-\frac{1}{2}(\tilde{k}_r^2 + \tilde{k}_z^2) + \frac{1}{2}\alpha_r \tilde{r}^2 + \frac{1}{2}\alpha_z \tilde{z}^2 + \tilde{F}\tilde{z} + \Gamma|\alpha|^2e^{-\frac{\tilde{r}^2}{\tilde{w}^2}}\cos^2{\tilde{z}}  + g|\Psi|^2  \right]\Psi
\end{equation}

\begin{equation}
\frac{\partial \alpha}{\partial t} = i\left[\tilde{\Delta_c} - 2\pi\Gamma\int e^{-\frac{\tilde{r}^2}{\tilde{w}^2}}\cos^2{\tilde{z}}|\Psi(\tilde{r},\tilde{z}, t)|^2 \tilde{r}d\tilde{r}d\tilde{z}\right]\alpha + \tilde{\eta} - \tilde{\kappa} \alpha
\end{equation}

$$\alpha_{r, z} = \frac{m^2 \omega^2_{r, z}}{\hbar^2 k^4} = \frac{\omega_{r,z}^2}{\omega_R^2}; \tilde{F} = \frac{Fm}{\hbar^2 k^3} = \frac{\omega_B}{2\pi\omega_R}; \Gamma = \frac{g_0^2}{\Delta_a\omega_R}; \tilde{g} = \frac{Uk^3}{E_R} = 4\pi k a_s$$ $$\tilde{\Delta}_c = \frac{\Delta_c}{\omega_R} ;\tilde{\eta} = \frac{\eta}{\omega_R} ;\tilde{\kappa} = \frac{\kappa}{\omega_R}  $$


For Ground state:

\begin{equation}
i\frac{\partial \Psi}{\partial t} = \left[-\frac{1}{2}(\tilde{k}_r^2 + \tilde{k}_z^2) + \frac{1}{2}\alpha_r \tilde{r}^2 + \frac{1}{2}\alpha_z \tilde{z}^2 + \tilde{F}\tilde{z} + \Gamma|\alpha|^2e^{-\frac{\tilde{r}^2}{\tilde{w}^2}}\cos^2{\tilde{z}}  + g|\Psi|^2  \right]\Psi
\end{equation}

\begin{equation}
\frac{\partial \alpha}{\partial t} = i\left[\tilde{\Delta_c} - 2\pi\Gamma\int e^{-\frac{\tilde{r}^2}{\tilde{w}^2}}\cos^2{\tilde{z}}|\Psi(\tilde{r},\tilde{z}, t)|^2 \tilde{r}d\tilde{r}d\tilde{z}\right]\alpha + \tilde{\eta} - \tilde{\kappa} \alpha
\end{equation}

$$\alpha_{r, z} = \frac{m^2 \omega^2_{r, z}}{\hbar^2 k^4} = \frac{\omega_{r,z}^2}{\omega_R^2}; \tilde{F} = \frac{Fm}{\hbar^2 k^3} = \frac{\omega_B}{2\pi\omega_R}; \Gamma = \frac{g_0^2}{\Delta_a\omega_R}; \tilde{g} = \frac{Uk^3}{E_R} = 4\pi k a_s$$ $$\tilde{\Delta}_c = \frac{\Delta_c}{\omega_R} ;\tilde{\eta} = \frac{\eta}{\omega_R} ;\tilde{\kappa} = \frac{\kappa}{\omega_R}  $$

