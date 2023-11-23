---
title: 'HofstadterTools: A Python package for analyzing the Hofstadter model'
tags:
  - Python
  - theoretical physics
  - electromagnetism
  - fractals
  - quantum physics
  - topological materials
  - mathematical physics
  - topological order
  - condensed matter physics
  - quantum Hall effect
  - band theory
authors:
  - name: Bartholomew Andrews
    orcid: 0000-0002-9079-7433
    affiliation: 1
affiliations:
 - name: Department of Physics, University of California, Berkeley, USA
   index: 1
date: 22 November 2023
bibliography: paper.bib
---

# Summary

*A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.*

The Hofstadter model successfully describes the behavior of non-interacting quantum particles hopping on a lattice coupled to a perpendicular gauge field, and hence is ubiquitous in many fields of research, including condensed matter, optical, and atomic physics. Motivated by this, we introduce [HofstadterTools](https://hofstadter.tools), a Python package that can be used to analyze the energy spectrum of a generalized Hofstadter model on any regular Euclidean lattice. The package can be applied to compute key properties of the band structure, such as quantum geometry and topology, as well as plot colored Hofstadter butterflies and Wannier diagrams.

# Statement of need

*A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.*

The Hofstadter model is one of many tight-binding models that is of interest to physicists and hence, it is often treated as an add-on to other numerical physics packages, such as pyqula, DiagHam, and TeNPy, or simply included as supplementary code together with research articles. However, the Hofstadter model's generalizability, interdisciplinary appeal, and recent experimental realization, motivates us to create a dedicated package that can provide a detailed analysis of its band structure, in the general case.

1) **Generalizability.** Although the Hofstadter model is typically thought of in its original context of electrons in a periodic potential coupled to a perpendicular magnetic field, this set-up is not required. For example, the vector potential can arise from an artificial gauge field in optical flux lattices with bosons, or from time-modulated Floquet Hamiltonians. Moreover, the full scope of the Hofstadter model is a topic of current research, with papers on its manifestation in hyperbolic lattices, higher-dimensional lattices, and synthesized materials, all released within the last couple of years.    

2) **Interdisciplinary appeal.** Aside from the Hofstadter model's well-known connection to condensed matter physics and the quantum Hall effect, it is also intensively studied in other fields. In mathematics, the Harper equation is a special case of an almost Mathieu operator, which is one of the best-understood types of ergodic Schrödinger operators. Moreover, the butterfly spectrum has a clear connection to the fields of chaos and fractals. In physics, the Hofstadter model can also be realized using laser-assisted hoppings in atomic, molecular, and optical physics, or via Floquet pumping.

3) **Recent experimental realization.** Despite the fact that the Hofstadter model was first studied 70 years ago, it has only been experimentally realized within the last decade. Signatures of the Hofstadter butterfly were first observed in moiré materials, and then later found in optical flux lattice and diffraction experiments. Not only does this spur recent theoretical interest, but it also increases the likelihood of newcomers entering the field, with the need for a well-documented, consolidated code repository that they can use for benchmarking experimental data or related computations.

The most likely use-case of HofstadterTools is for computational condensed matter theorists to use as a reference point alongside many-body computations. The Hofstadter model is an infinitely configurable topological flat band model and hence is a popular choice for theorists studying strongly-correlated phenomena, such as the fractional quantum Hall effect and superconductivity. Since there is a relationship between properties of the single-particle band structures and the stability of exotic strongly-correlated states, HofstadterTools may be used to guide theorists who are analyzing quantum many-body systems. 

# Acknowledgements

We thank Gunnar Möller, Titus Neupert, Rahul Roy, Alexey Soluyanov, Mike Zaletel, Johannes Mitscherling, and Mathi Raja, for useful discussions. This project was funded by the Swiss National Science Foundation under Grant No. P500PT_203168.

# References