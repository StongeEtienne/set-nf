SET-FLOW nextflow
====================================================
Apply "Surface Flow", to the given Tractoflow result

A) Convert images to Nifti and surfaces to vtk files
B) Register surfaces to the original T1 space
C) Generate surface ROI and masks from labels, and concatenates them in a single files
D) Register surfaces to the diffusion space (fodf), with tractflow ANTs transform
E) Surface smoothing and "Surface Flow"
F) Compute intersection and combine tractogram with flow

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
