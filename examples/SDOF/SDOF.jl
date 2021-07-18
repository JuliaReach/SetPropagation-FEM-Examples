
# We consider as an illustrative example a simple harmonic oscillator system.
# The set propagation definitions introduced so far are illustrated using a single
# degree of freedom second order problem given by the following equation:

# ```math
# u''(t) + ω^2 u(t) = 0,
# ```
# where $ω$ is a scalar and $u(t)$ is the unknown.
#
# This problem can be associated with a spring-mass system, where $u(t)$ is the
# elongation of the spring at time $t$ and the mass and stiffness set a natural
# frequency $ω$. In this case we consider $ω = 4π$.

# Let us introduce the state variable $v(t) := u'(t)$ and define the vector
# $x(t) = [u(t), v(t)]^T$. Then we can rewrite the system in the following first order form
#
#
# ```math
# x'(t) = \begin{pmatrix} 0 & 1 \\ -\omega^2 & 0 \end{pmatrix} x(t)
#
# Two cases are considered for the initial conditions, first a single initial condition
# given by $u(0)=1$ and $v(0)=0$, and a second case where the initial conditions belong to a set.

# This problem is solved using set propagation and the results obtained for the
# single initial condition are shown in the following figures.
#
# The time step-size used is $\delta = 0.025$.
#
# The set $\Omega_0$ includes the analytic trajectory within $[0, \delta]$, and
# such set is propagated to cover the solution for further time intervals as shown below.
#
# It is worth noting that even when the initial condition is a singleton,
# the set propagation technique returns a sequence of sets, and such sequence guarantees
# an enclosure of the analytic solution at \emph{any} time $t \in [0, T] \subset \mathbb{R}$.
#
# For comparison, we also plot the analytic solution (magenta solid line) and a
# numerical solution obtained using the Newmark's method with time step $\delta/5$
# (red triangle markers).

# The result obtained when a set of initial conditions, together with few dozen
# trajectories with random initial conditions drawn from the initial set
# $\mathcal{X}_0 = \mathcal{U}_0 \times \mathcal{V}_0 = [0.9, 1.1] \times [-0.1, 0.1]$,
# where $\times$ denotes the Cartesian product.

# In both cases, it is shown that a single set-based integration covers infinitely
# many trajectories \emph{in dense time}, i.e. for all intermediate times between
# $0$ and $T$, where $T > 0$ is the time horizon, there is a reach-set that covers
# all the exact solutions.
