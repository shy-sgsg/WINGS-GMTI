# Stage 1 algorithm integration report

- angle grid: -60 to 60
- pulse_resample_mode: physical_time
- range_resize_mode: fft_zero_pad
- input period_count generated: 1
- output file size: 1501751680
- suggested smoke command: `./build/GMTI_core outputs/stage1/config/temp_config_stage1_newsystem.xml`
- smoke command executed: `./build/GMTI_core --trig-mode math outputs/stage1/config/temp_config_stage1_newsystem.xml > outputs/stage1/logs/run_smoke_test.log 2>&1`
- POS/velocity header: regenerated packets write UTC, interpolated lat/lon/alt, and vn/ve/vd derived from adjacent POS position differences. Header spot checks show speed around 50.14 m/s instead of 0.
- current status: data generation and static validation complete; full runtime algorithm logs should be appended after running GMTI_core or GMTI_pipe_core on a CUDA-capable machine.
- stage conclusion: no cooperative target is injected in stage 1; mirror background is for engineering link and statistics validation, not strict wide-angle geometry proof.
