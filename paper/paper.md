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
date: 23 November 2023
bibliography: paper.bib
---

# Summary

The Hofstadter model successfully describes the behavior of non-interacting quantum particles hopping on a lattice coupled to a perpendicular gauge field, and hence is ubiquitous in many fields of research, including condensed matter, optical, and atomic physics. Motivated by this, we introduce HofstadterTools ([https://hofstadter.tools](https://hofstadter.tools)), a Python package that can be used to analyze the energy spectrum of a generalized Hofstadter model on any regular Euclidean lattice. The package can be applied to compute key properties of the band structure, such as quantum geometry and topology, as well as plot colored Hofstadter butterflies and Wannier diagrams.

# Statement of need

The purpose of HofstadterTools is to consolidate the fragmented theory and code relevant to the Hofstadter model into one well-documented Python package, which can be used easily by non-specialists as a benchmark and springboard for their own research projects. The Hofstadter model [@Harper55; @Azbel64; @Hofstadter76] is one of many tight-binding models that is of interest to physicists and hence, it is often treated as an add-on to larger numerical packages, such as pyqula [@pyqula], DiagHam [@DiagHam], and TeNPy [@TeNPy], or simply included as supplementary code together with research articles [@Bodesheim23]. However, the Hofstadter model's generalizability, interdisciplinary appeal, and recent experimental realization, motivates us to create a dedicated package that can provide a detailed analysis of its band structure, in the general case.

1) **Generalizability.** The Hofstadter model was originally studied in the context of electrons hopping in a periodic potential coupled to a perpendicular magnetic field. However, the model transcends this framework and can be extended in numerous directions. For example, the Peierls phases that arise in the Hamiltonian due to the magnetic field [@Peierls33] can also be generated using artificial gauge fields [@Goldman14] or Floquet modulation [@Eckardt17]. Moreover, the full scope of the Hofstadter model is yet to be unearthed, with papers on its application to hyperbolic lattices [@Stegmaier22], higher-dimensional crystals [@DiColandrea22], and synthesized materials [@Bodesheim23], for example, all published within the last couple of years.    

2) **Interdisciplinary appeal.** Owing to its generalizability, interest in the Hofstadter model goes beyond its well-known connection to condensed matter physics and the quantum Hall effect [@Avron03]. In mathematics, for example, the difference relation arising in the solution of the Hofstadter model, known as the Harper equation [@Harper55], is a special case of an almost Mathieu operator, which is one of the best-understood types of ergodic Schrödinger operator [@Simon00; @Avila09]. Moreover, in other branches of physics, the Hofstadter model has growing relevance in a variety of subfields, including: cold atomic gases [@Cooper19], acoustic metamaterials [@Ni19], and photonics [@Zilberberg18].  

3) **Recent experimental realization.** Although the Hofstadter model was first studied last century [@Peierls33; @Harper55], it has only been experimentally realized within the last decade. Signatures of the Hofstadter spectrum were first observed in moiré materials [@Dean13] and optical flux lattices [@Aidelsburger13], and they have since been reproduced in several other experimental platforms [@Cooper19; @Ni19; @Zilberberg18]. Not only does this spur recent theoretical interest, but it also increases the likelihood of newcomers entering the field, with the need for a user-friendly code repository that can be quickly applied to benchmark experimental data and related computations.

A typical use-case of HofstadterTools is for condensed matter theorists to utilize as a reference point alongside many-body computations. The Hofstadter model is an infinitely-configurable topological flat-band model and hence, is a popular choice among theorists studying strongly-correlated phenomena, such as the fractional quantum Hall effect [@Andrews20; @Andrews21] and superconductivity [@Shaffer21; @Sahay23]. Since there is a relationship between properties of the single-particle band structures and the stability of exotic strongly-correlated states [@Jackson15; @Andrews23], HofstadterTools may be used to guide theorists who are analyzing quantum many-body systems. This scenario notwithstanding, we hope that HofstadterTools will find many interdisciplinary applications, and we look forward to expanding the package in these directions, with help from the community.   

# Acknowledgements

We thank Gunnar Möller, Titus Neupert, Rahul Roy, Alexey Soluyanov, Mike Zaletel, Johannes Mitscherling, and Mathi Raja, for useful discussions. This project was funded by the Swiss National Science Foundation under Grant No. P500PT_203168.

# References