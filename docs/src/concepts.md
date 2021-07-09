# # Reachability concepts

In this section we consider set-based integration in the context of linear differential
equations and introduce two essential concepts: reach-set and flowpipe.
Then we define the representations used in our approach: hyperrectangles, zonotopes and
support functions. The definitions are illustrated by means of a textbook example:
the harmonic oscillator. More complex and larger problems are considered in other sections.

# ## Trajectories, reach-sets and flowpipes

Let us consider an initial-value problem

```math
\begin{equation}
\dot{\bfx}(t) = \bfA \bfx(t),\qquad \bfx(0) \in \mcX_0,~t \in [0, T],
\end{equation}
```
with state matrix $\bfA \in \mathbb{R}^{n\times n}$ and where $\bfx(t) \in \mathbb{R}^n$
is the state vector at time $t$. The initial condition is $\bfx(0) \in \mcX_0\subset \mathbb{R}^n$
with $\mcX_0$ a closed and convex set. In this article we will refer to any state belonging to $\mcX_0$ as
an \textit{admissible} state.

The guiding principle in set-based integration is to represent reachable states
using sets in the Euclidean space that are *propagated* according to the system's dynamics.
The obtained sequence of sets, by construction, covers the exact trajectories for given
time intervals, for any admissible initial condition or external actions.
For instance, the result for the harmonic oscillator presented in Section XXX
using a single initial condition is shown in Fig.~\ref{fig:sdof_singleton} (green circle),
and in Fig.~\ref{fig:sdof_distributed} (green box) the results for a box of initial conditions are shown.

Given the system of ODEs in Eq.~\eqref{eq:linearODE}, a \emph{trajectory}
from an initial state $\bfx_0$ is the solution for a given initial condition
$\bfx(0) = \bfx_0$, denoted $\varphi(\bfx_0, t) : \R^n \times [0, T] \rightarrow \R^n$,
which we know is $\varphi(\bfx_0, t) = e^{\bfA t}\bfx_0$.

The \textit{reach-set} at a \textit{time point} $t \in [0, T]$ is
\begin{equation}\label{eq:reachset_point}
\mcR^e(\mcX_0, t) := \bigl\{ \varphi(\bfx_0, t): \bfx_0 \in \mcX_0 \bigr\}.
\end{equation}
This is the set of states which are reachable at time $t$ starting from any admissible initial condition.
%
Generalizing Eq.~\eqref{eq:reachset_point} to time intervals is straightforward: the reach-set over $[t_0, t] \subseteq [0, T]$ is
\begin{equation}\label{eq:reachset_interval}
\mcR^e(\mcX_0, [t_0, t]) := \bigl\{ \varphi(\bfx_0, s): \bfx_0 \in \mcX_0, \ \forall~s \in [t_0, t] \bigr\}.
\end{equation}
The difference between reach-sets for time points and for time intervals is that the former are evaluated at a particular time $t$, while the latter include \textit{all} states reachable for \textit{any} time between $t_0$ and $t$, hence $\mcR^e(\mcX_0, t) \subseteq \mcR^e(\mcX_0, [t_0, t])$.
Finally, given a time step $\delta > 0$ and a collection of reach-sets $\left\{\mcR^e(\mcX_0, [\delta(k-1), \delta k)]\right\}_{k=1}^N$ with $N = T/\delta$, the \textit{flowpipe} of Eq.~\eqref{eq:linearODE} is the set union,
\begin{equation}\label{eq:flowpipe}
\mcF^e(\mcX_0, [0, T]) := \bigcup_{k=1}^N \mcR^e(\mcX_0, [\delta(k-1), \delta k)]).
\end{equation}

In Eqs.~\eqref{eq:reachset_point}-\eqref{eq:flowpipe}, the superindex $e$ is used to remark that the definitions refer to the \textit{exact} or true sets. In practice, reachable sets can only be obtained approximately and several algorithms to numerically approximate reachable sets are known in the literature. Different methods crucially depend on how the sets are represented and the cost of operating with those representations.


# ## Set representations

Three common set representations used in reachability analysis of linear systems are zonotopes, hyperrectangles and support functions. Let us briefly recall their definition.

# ### Zonotopes

Given a vector $\bfc \in \mathbb{R}^n$ that we call \emph{center} and a list of $p \ge n$ vectors that we call \emph{generators}, $\bfg_1, \ldots, \bfg_p$, $\bfg_i\in \mathbb{R}^n$, the associated zonotope is the set
\begin{equation} \label{eq:zonotope}
Z = \left\{\bfx \in \mathbb{R}^n : \bfx = \bfc + \sum_{i=1}^p \xi_i \bfg_i,\quad \xi_i \in [-1, 1] \right\}.
\end{equation}
In compact notation we write $Z = \langle \bfc, \bfG\rangle_Z$, where each generator is stored as a column of the generators matrix $\bfG \in \mathbb{R}^{n \times p}$. Zonotopes can be interpreted in different ways, such as the Minkowski sum of line segments, or as the affine map of a unit ball in the infinity norm. The zonotope representation has been successfully used in reachability analysis because it offers a very compact representation relative to the number of vertices or faces, so it is well suited to represent high dimensional sets \citep{GLGM06}.

# ### Hyperrectangles

Given vectors $\bfc \in \mathbb{R}^n$ and $\bfr \in \mathbb{R}^n$ that we call \emph{center} and \emph{radius} respectively, the associated hyperrectangle is the set
\begin{equation} \label{def:hyperrectangle}
H = \left\{\bfx \in \mathbb{R}^n : \vert x_i - c_i \vert \leq r_i, \quad i = 1, \ldots, n\right\}.
\end{equation}
In compact notation we write $H = \langle \bfc, \bfr \rangle_H$. Every hyperrectangle is a zonotope but the converse is not true. In fact, the linear map of a hyperrectangle can be (exactly) represented as a zonotope.
%
While computations with hyperrectangles can be performed very efficiently, using hyperrectangles for reachability can result in coarse overapproximations unless they are used in concurrence with other representations at intermediate parts of the algorithm.

# ### Support functions

The support function of a set $X \subset \mathbb{R}^n$ along direction $\bfd \in \mathbb{R}^n$ is the scalar
\begin{equation}\label{eq:support_function}
\rho(\bfd, X) = \max_{\bfx \in X} \bfd^T \bfx,
\end{equation}
where the superscript $T$ denotes transposition.
The set of points in $X$ that are the maximizers of Eq.~\eqref{eq:support_function} are called the \textit{support vectors}. The intuition behind Eq.~\eqref{eq:support_function} is that support functions represent the farthest (signed) distance to the origin of the set $X$ along direction $\bfd$. Given that the direction can be chosen at will, support functions can be used to obtain the solution corresponding to a linear combination of the state variables at once and at a reduced cost~\citep{LeGuernic09}.

%  it is actually not required to compute the flowpipe for all variables~

# ## Set propagation

# ### Propagation using hyperrectangles

# ### Propagation using zonotopes

# ### Propagation using support functions
