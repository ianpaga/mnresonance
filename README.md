Neutrino flavor evolution - Matter neutrino resonances in neutron star mergers
====

Description: Neutrino flavor evolution in the remnant of two neutron start. This simulation assumes spatial homogeneity, a single neutrino energy, and axial symmetry.

## Journal reference: 
[Ian Padilla-Gay et al JCAP05(2024)037](https://iopscience.iop.org/article/10.1088/1475-7516/2024/05/037), [ePrint:2403.15532](https://arxiv.org/abs/2403.15532)
Results from [ERDA AstroNu Data Archive](https://sid.erda.dk/share_redirect/e2zTyjhG3B/index.html), see also Figures in the paper.

## Figures & animations that appear in my publication:
![multi_angle_iso_per_panels-1](https://github.com/ianpaga/mnresonance/assets/57350668/86e4fa52-4f75-49ef-a4ad-382f74e1b3b9)
![polar_iso_per](https://github.com/ianpaga/mnresonance/assets/57350668/c9d9f41f-522a-4f45-b5e8-c6a45c6b21dd)

## More animations can be found here:
[ERDA AstroNu Data Archive](https://sid.erda.dk/share_redirect/e2zTyjhG3B/index.html)

## Requirements:
- C++ compiler
- [Boost Library](https://www.boost.org/)
- OpenMP
- Python, Matplotlib, NumPy

## Compiling and running:
- Run ./chpc to compile
- Run executable *.out
- Outputs *.raw and *.sum are large files. Use the bash script datfiles.sh to slice the data into smaller *.dat files
- Plot results with python potentials.py (see examples in /plots)
