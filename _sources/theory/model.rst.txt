Model
=====

We consider a system of non-interacting quantum particles hopping on a lattice in the xy-plane and coupled to a gauge field. For the sake of concreteness, we discuss electrons of charge :math:`-e` and mass :math:`m_\mathrm{e}` in this section, governed by a 2D tight-binding model and coupled to a perpendicular magnetic field :math:`\mathbf{B}=B\hat{\mathbf{e}}_z`. However, the theory is readily generalizable.

The **generalized Hofstadter model** describing the electron motion is given by

.. math::
	H = -\sum_{\kappa} \sum_{\langle ij \rangle_\kappa} t_\kappa e^{\mathrm{i}\theta_{ij}} c^\dagger_i c_j + \mathrm{H.c.},

where :math:`\langle \dots \rangle_\kappa` denotes :math:`\kappa`-th nearest neighbors on some regular Euclidean lattice (defined below) with corresponding hopping amplitudes :math:`t_\kappa`, :math:`e^{\mathrm{i}\theta_{ij}}` is the Peierls phase acquired by hopping from site :math:`i` to site :math:`j`, and :math:`c^{(\dagger)}` is the electron (creation)annihilation operator :cite:`Harper55, Azbel64, Hofstadter76`. We note that a set of nearest neighbors in HofstadterTools is defined as a set of sites that are *equidistant* from an origin site, i.e. a group of points that lie on a circle.

The **Peierls phase** is the lattice-analogue of the Aharonov-Bohm phase and is given as

.. math::
	\theta_{ij} = \frac{2\pi}{\phi_0} \int_i^j \mathbf{A}\cdot \mathrm{d}\mathbf{l},

where :math:`\phi_0` is the flux quantum, :math:`\mathbf{A}` is the vector potential, and :math:`\mathrm{d}\mathbf{l}=(\mathrm{d}x,\mathrm{d}y)` is an infinitesimal line element from site :math:`i=(X_i, Y_i)` to site :math:`j=(X_j, Y_j)` :cite:`Peierls33`. There is freedom in choosing the vector potential, such that :math:`\mathbf{B} = \nabla \times \mathbf{A} = B\hat{\mathbf{e}}_z`.

.. note::
	Although there is freedom in choosing the vector potential, there may be strategy involved. Several investigations define an artificial Cartesian coordinate system to set up the Harper equation, instead of constructing separate equations for the basis sites :cite:`Oh00`. For example, for the honeycomb lattice with lattice constant :math:`a`, physicists often define the coordinates :math:`(x,y)=(mb,nc)`, where :math:`b=a/2` and :math:`c=\sqrt{3}a/6` :cite:`Andrews20`. The advantage of this approach is that it allows you to address every site on the lattice using a single discrete coordinate system. In order to construct the Harper equation, we then eliminate one of these coordinates using the plane wave ansatz, and are left with a difference equation in the remaining direction. If it is not possible to construct a discrete coordinate system in one direction because, for example, the ratio of site spacings is irrational, then it is strategic to choose the vector potential in the orthogonal direction. Moreover, when the Harper equation is defined using such artificial coordinates, then the gauge may have an impact on the periodicity of the spectrum. The simplest example of this is for the nearest-neighbor Hofstadter model on the triangular lattice with :math:`b=a/2` and :math:`c=\sqrt{3}a/2`. A vector potential with a magentic unit cell extended in the x-direction yields an aperiodic spectrum, when :math:`n_\phi=BA_\mathrm{UC}/\phi_0`, whereas a vector potential with a magnetic unit cell extended in the y-direction yields a periodic spectrum (try it!). We note that when addressing the lattice in terms of unit cells and basis sites, as in HofstadterTools, the vector potential is a genuine gauge choice with no unintended side-effects. Nevertheless, inspired by this observation, since our :math:`\mathbf{a}_1=(a, 0)` vector is fixed, our angle :math:`\theta` is defined in the positive direction, and hence a discrete coordinate system in the x-direction is generally impossible to construct for the underlying Bravais lattice, we work with the Landau gauge in the y-direction, with conserved :math:`k_x` momentum.

In HofstadterTools, we choose **Landau gauge** in the y-direction, with a conserved :math:`k_x` momentum, such that :math:`\mathbf{A}=-By\hat{\mathbf{e}}_x`. We note that this is sometimes referred to as Landau gauge in the x-direction, depending on whether you are referring to the dimensions of the magnetic unit cell or the components of the vector potential. Parameterizing the expression for the Peierls phase allows us to write

.. math::
	\theta_{ij} = -\frac{2\pi B}{\phi_0} \int_0^1
	\begin{pmatrix} Y_i + (Y_j - Y_i)\tau \\ 0 \end{pmatrix} \cdot \begin{pmatrix} X_j - X_i \\ Y_j - Y_i \end{pmatrix} \mathrm{d}\tau,

which reduces to

.. math::
	\theta_{ij} = -\frac{2\pi B}{\phi_0} \Delta X \left( Y_i + \frac{\Delta Y}{2} \right),

where :math:`\Delta X = X_j - X_i` and :math:`\Delta Y = Y_j - Y_i`. Notice that this expression depends on *absolute* y-coordinates but only *relative* x-coordinates. The effect of the Peierls transformation is that it enlarges the unit cell (UC) of dimensions :math:`1\times 1` to a **magnetic unit cell (MUC)** of dimensions :math:`1\times q`. The lattice translation operators then form the **magnetic translation group** :cite:`Zak64`. We note that although the dimensions of the MUC is a gauge choice, its area, :math:`q`, is not.

In the Hofstadter model, the magnetic field strength is often expressed in terms of a **flux density** :math:`n_\phi`. Typically, this density is defined with respect to the UC area :math:`A_\mathrm{UC}`, such that :math:`n_\phi=B A_{\mathrm{UC}}/\phi_0`. However, this is not always the case. In cases where the minimal area enclosed by electron hopping :math:`A_\mathrm{min}` is smaller than :math:`A_\mathrm{UC}` or an effective UC area spanned by the hopping terms, then it is necessary to define the flux density with respect to :math:`A_\mathrm{min}` in order to guarantee the correct periodicity and reveal the entire spectrum. For example, for the nearest-neighbor Hofstadter model on the triangular lattice we have :math:`A_\mathrm{UC}=a^2` but :math:`A_\mathrm{min}=a^2/2`, and so we need to define :math:`n_\phi=B A_{\mathrm{min}}/\phi_0` in order to reveal the entire spectrum. This equates to dividing the flux density by a **periodicity factor** of 2.

.. note::
	Alternatively, the periodicity of the spectrum may be recovered by performing a unitary transformation of the wavefunctions in the Schrödinger equation. This is sometimes referred to a Rammal transformation :cite:`Rammal85`. For the sake of simplicity, we have made the decision not to enforce periodicity in this way in HofstadterTools.

Following a suitable definition of :math:`n_\phi`, we can substitute this into our expression for the Peierls phase. Since the Peierls factor is a complex phase factor, we consider rational flux densities :math:`n_\phi=p/q`, where :math:`p` and :math:`q` are coprime integers. Moreover, since the denominator of the flux density :math:`q` is the MUC area in units of UCs, there will ultimately be :math:`N_\mathrm{b}q` bands in the spectrum, where :math:`N_\mathrm{b}` is the number of sites in the basis.

.. image:: ../images/theory/lattice.png
	:align: center
	:width: 60%

In the figure above, we show an example lattice annotated with the relevant unit cells. This figure also serves to define the variables in HofstadterTools. We construct a lattice by repeating some basis in multiples of the Bravais vectors :math:`\mathbf{a}_1=(a,0)` and :math:`\mathbf{a}_2=\alpha a (\cos\theta, \sin\theta)`, where :math:`\alpha` and :math:`\theta` (not to be confused with the Peierls phase) are measures of the Bravais lattice **anisotropy** and **obliqueness**, respectively. In addition, the Bravais lattice has a collection of :math:`N_\mathrm{b}` basis sites at positions :math:`\{\mathbf{a}_\mathrm{b}\}` relative to the UC origin. The UC is defined as the span of :math:`\{\mathbf{a}_1, \mathbf{a}_2\}` and the MUC is defined as the span of :math:`\{\mathbf{a}_{\mathrm{MUC},1}, \mathbf{a}_{\mathrm{MUC},2}\}`, where :math:`\mathbf{a}_{\mathrm{MUC},1}=\mathbf{a}_1` and :math:`\mathbf{a}_{\mathrm{MUC},2}=q\mathbf{a}_2`. If the appropriate hoppings exist, there may be a distinct minimal hopping plaquette area :math:`A_\mathrm{min}`, as shown in the example. By choosing the lattice constant, anisotropy, obliqueness, and a set of basis vectors, we can construct any regular Euclidean lattice.
