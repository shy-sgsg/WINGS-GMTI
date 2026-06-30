# Stage 3 SAR scene simulation report

This scaffold records the adjusted third-stage plan: high-resolution SAR image scatterer inversion plus new-system 61-beam dual-channel DDC forward generation.

## Requested run

- stage3_config: stage3_config.json
- sar_image: data/sar_highres.tif
- georef_path: 
- scatterer_csv: 
- output_dir: outputs/stage3
- scene_mode: sar
- extract_only: false
- forward_only: false
- target_enabled: true

## Required limitations to retain in final validation

1. Intensity-only SAR images cannot recover true coherent phase; random phase is a statistical approximation.
2. SAR imaging geometry differs from the new 61-beam scan geometry, so extracted scatterers are a reflectivity-scene estimate, not raw echo inversion.
3. Georeferencing accuracy directly controls spatial correctness.
4. More scatterers improve background realism but increase compute cost quickly.
5. The first platform model is ideal_straight; navigation disturbance and attitude error are future extensions.
6. The first beam model is Gaussian; measured antenna patterns can replace it later.
7. The first dual-channel phase model is baseline_approx; a baseline-vector model is the next upgrade.
