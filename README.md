SET nextflow pipeline
=====================
Run the "Surface-Enhanced Tractography" (SET) nextflow pipeline

A) Convert images & surfaces files and Register them to the original T1 space
B) Generate surface ROI and masks from labels, and concatenates them in a single files
C) Register surfaces to the diffusion space (fodf), with tractflow ANTs transform
D) Surface smoothing and "Surface Flow"
E) Generate seeding maps and vertices to initialize tractography
F) "Surface-Enhanced Tractography"
G) Combine intersections
H) Compute endpoints maps (connectivity and surface density)


If you use this pipeline, please cite:
```
St-Onge, Etienne, et al. "Surface-enhanced tractography (SET)." NeuroImage 169 (2018): 524-539.
```

Requirements
------------

- [Nextflow](https://www.nextflow.io)
- Dipy
- Scilpy
- TriMeshPy
- VTK (python)

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`
