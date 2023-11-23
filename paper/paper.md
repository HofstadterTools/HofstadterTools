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

The Hofstadter model successfully describes the behavior of non-interacting quantum particles hopping on a lattice coupled to a perpendicular gauge field, and hence is ubiquitous in many fields of research, including condensed matter, optical, and atomic physics. Motivated by this, we introduce HofstadterTools ([https://hofstadter.tools](https://hofstadter.tools)), a Python package that can be used to analyze the energy spectrum of a generalized Hofstadter model on any regular Euclidean lattice. The package can be applied to compute key properties of the band structure, such as quantum geometry and topology, as well as plot colored Hofstadter butterflies and Wannier diagrams.

# Statement of need

The purpose of HofstadterTools is to consolidate the fragmented theory and code relevant to the Hofstadter model into one well-documented Python package, which can be used easily by non-specialists as a benchmark and springboard for their own research projects. The Hofstadter model [@Harper55; @Azbel64; @Hofstadter76] is one of many tight-binding models that is of interest to physicists and hence, it is often treated as an add-on to larger numerical physics packages, such as pyqula [@pyqula], DiagHam [@DiagHam], and TeNPy [@TeNPy], or simply included as supplementary code together with research articles [@Bodesheim23]. However, the Hofstadter model's generalizability, interdisciplinary appeal, and recent experimental realization, motivates us to create a dedicated package that can provide a detailed analysis of its band structure, in the general case.

1) **Generalizability.** The Hofstadter model is typically thought of in its original context of electrons in a periodic potential coupled to a perpendicular magnetic field, however this set-up is not necessary. For example, the vector potential can arise from an artificial gauge field in optical flux lattices with bosons [@Goldman14], or from time-modulated Floquet Hamiltonians [@Eckardt17]. Moreover, the full scope of the Hofstadter model is a topic of current research, with papers on its manifestation in hyperbolic lattices [@Stegmaier22], higher-dimensional lattices [@DiColandrea22], and synthesized materials [@Bodesheim23], all released within the last couple of years.    

2) **Interdisciplinary appeal.** Aside from the Hofstadter model's well-known connection to condensed matter physics and the quantum Hall effect [@Avron03], it is also studied in other fields of research. In mathematics, the Harper equation is a special case of an almost Mathieu operator, which is one of the best-understood types of ergodic Schrödinger operators [@Simon00]. Moreover, the butterfly spectrum has a clear connection to the fields of chaos and fractals [@Cedzich23; @Avila09]. In physics, the Hofstadter model has growing relevance in the field of atomic, molecular, and optical physics [@Cooper19], and can also be manipulated using periodic pumping [@Eckardt17].

3) **Recent experimental realization.** The Hofstadter model was first studied 70 years ago [@Harper55], however it has only been experimentally realized within the last decade. Signatures of the Hofstadter butterfly were first observed in moiré materials [@Dean13], and then later found in optical flux lattice [@Aidelsburger13] and diffraction experiments [@DiColandrea22]. Not only does this spur recent theoretical interest, but it also increases the likelihood of newcomers entering the field, with the need for a user-friendly, well-documented code repository that they can quickly apply to benchmark experimental data or related computations.

A typical use-case of HofstadterTools is for computational condensed matter theorists to utilize as a reference point alongside many-body computations. The Hofstadter model is an infinitely-configurable topological flat-band model and hence, is a popular choice among theorists studying strongly-correlated phenomena, such as the fractional quantum Hall effect [@Andrews20; @Andrews21] and superconductivity [@Shaffer21; @Sahay23]. Since there is a relationship between properties of the single-particle band structures and the stability of exotic strongly-correlated states [@Jackson15; @Andrews23], HofstadterTools may be used to guide theorists who are analyzing quantum many-body systems. This scenario notwithstanding, we hope that HofstadterTools will find many interdisciplinary applications, and we look forward to expanding the package in these directions, with help from the community.   

# Acknowledgements

We thank Gunnar Möller, Titus Neupert, Rahul Roy, Alexey Soluyanov, Mike Zaletel, Johannes Mitscherling, and Mathi Raja, for useful discussions. This project was funded by the Swiss National Science Foundation under Grant No. P500PT_203168.

# References