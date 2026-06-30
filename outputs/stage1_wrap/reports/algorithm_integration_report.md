# Stage 1 algorithm integration report

- angle grid: -59 to 61
- pulse_resample_mode: normalized
- range_resize_mode: fft_zero_pad
- input period_count generated: 1
- output file size: 1501751680
- suggested smoke command: `./build/GMTI_core outputs/stage1_wrap/config/temp_config_stage1_newsystem.xml`
- POS/velocity header: generated packets write UTC, interpolated lat/lon/alt, and vn/ve/vd derived from adjacent POS position differences. This avoids zero platform speed in squint estimation.
- current status: data generation and static validation complete; runtime algorithm logs should be appended after running GMTI_core or GMTI_pipe_core.
- stage conclusion: no cooperative target is injected in stage 1; out-of-range background beams are reused by forward cyclic angle wrapping for engineering link and statistics validation, not strict wide-angle geometry proof.
