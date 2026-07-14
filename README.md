# TrophicAnalysis.jl

[![Build Status](https://github.com/pjkissack/TrophicAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/pjkissack/TrophicAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## The Trophic Analysis Toolbox in Julia

Trophic levels and incoherence in Julia.

This code [1] comprises two parts:

* ```trophic_tools```: A Julia version of Bazil Sansom's Trophic Analysis Toolbox (see: https://github.com/pjkissack/Trophic-Analysis-Toolbox.git), based on the python version [2], with additional functionality for calculating the original definition of trophic levels as in [3];
* ```trophic_tuning```: Code for the sampling of networks tuned to specified trophic coherence using the preferential preying [4] and generalised preferential preying [5] algorithms.

Please cite [1] upon use of this module, and additionally [2] upon use of the trophic_tools section.

#### REFERENCES

[1] Kissack, P., Drysdale, C. and Johnson, S. (202x).

[2] MacKay, R. S., Johnson, S., Sansom, B. (2020). How directed is a directed network?. R Soc Open Sci. 7 (9): 201138. https://doi.org/10.1098/rsos.201138

[3] Levine, S. (1980). Several measures of trophic structure applicable to complex food webs. 
J Theor Biol. 83 (2): 195-207. 
https://doi.org/10.1016/0022-5193(80)90288-X.

[4] Johnson, S., Domínguez-García, V., Donetti, L., & Muñoz, M. A. (2014). Trophic coherence determines food-web stability. Proc Natl Acad Sci U S A. 111 (50): 17923–17928. https://doi.org/10.1073/pnas.1409077111

[5] Klaise, J., Johnson, S. (2016). From neurons to epidemics: How trophic coherence affects spreading processes. Chaos 26 (6): 065310. https://doi.org/10.1063/1.4953160

