# Stage 1 parameter audit

## Old system

- data_ch1: /home/shy/AIR/小长/GMTI_Data/Mission068/RawData/ID13_8H13M_GMTI0
- data_ch2: /home/shy/AIR/小长/GMTI_Data/Mission068/RawData/ID13_8H13M_GMTI1
- format: float32 IQ with 256-byte PRT header per channel
- beams per period: 26
- pulses per beam: 256
- complex samples per pulse: 4096
- PRF/fs/fc/Tr: 2000 Hz / 75 MHz / 17 GHz / 38 us

## New system

- scan: -59 to 61 deg, step 2, beam_count 61
- 61 x 130 / 1300 = 6.1 s
- 80 percent realtime gate: 4.88 s
- strict engineering gate: 4.8 s
- sample_window_us x fs_new = 197 x 60 = 11820 points

## File audit

- ch1 bytes: 219807744
- ch2 bytes: 219807744
- bytes per period per channel: 219807744
- period_count: 1
- status: ok
