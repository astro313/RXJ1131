# determine flux density @ 4.88 GHz
noise: line-free regions = 1.3e-5 Jy/B

- Core

1. CASAregion_core1
- CASA viewer 2D fit: CASAregion_core1.log
- CASA imstat:
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(RXJ1131C_R0.FITS)
     Frequency Spectral_Value         Stokes BrightnessUnit       BeamArea 
  4.8851e+09Hz      4.8851GHz              I        Jy/beam        19.2738 
          Npts            Sum    FluxDensity           Mean            Rms 
           247   1.901721e-02   9.866866e-04   7.699276e-05   1.628717e-04 
       Std dev        Minimum        Maximum   region count 
  1.438160e-04  -2.983379e-05   7.394222e-04              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
2. CASAregion_core2 (-- Use this value --)
- CASA viewer 2D fit: CASAregion_core2.log
- CASA imstat: 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(RXJ1131C_R0.FITS)
     Frequency Spectral_Value         Stokes BrightnessUnit       BeamArea 
  4.8851e+09Hz      4.8851GHz              I        Jy/beam        19.2738 
          Npts            Sum    FluxDensity           Mean            Rms 
            86   1.669505e-02   8.662040e-04   1.941285e-04   2.729402e-04 
       Std dev        Minimum        Maximum   region count 
  1.929858e-04   1.435694e-05   7.394222e-04              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
- rms = 1.3e-5 (line-free) * sqrt(86/19.3) = 27.44 µJy

- arc
1. CASAregion_arc
- CASA viewer 2D fit: structure too complex, probably not a good idea to fit
- CASA imstat
- 1.545142e-03 +/- 5.9148825596161243e-05 [Jy]
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(RXJ1131C_R0.FITS)
     Frequency Spectral_Value         Stokes BrightnessUnit       BeamArea 
  4.8851e+09Hz      4.8851GHz              I        Jy/beam        19.2738 
          Npts            Sum    FluxDensity           Mean            Rms 
           399   2.978078e-02   1.545142e-03   7.463855e-05   1.120984e-04 
       Std dev        Minimum        Maximum   region count 
  8.374194e-05  -2.213363e-05   4.782451e-04              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
2. CASAregion_arc2.crtf (traces 3sigma & avoid jet, use this)
- CASA imstat
- 1.273277e-03 +/- 4.2085778291312834e-05 [Jy] 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(RXJ1131C_R0.FITS)
     Frequency Spectral_Value         Stokes BrightnessUnit       BeamArea 
  4.8851e+09Hz      4.8851GHz              I        Jy/beam        19.2738 
          Npts            Sum    FluxDensity           Mean            Rms 
           202   2.454090e-02   1.273277e-03   1.214896e-04   1.518474e-04 
       Std dev        Minimum        Maximum   region count 
  9.131924e-05   1.231700e-05   4.782451e-04              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
3. CASAregion_arc3.crtf
- CASA imstat
- 1.561674e-03 +/- 6.3232754345965799e-05 [Jy]
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(RXJ1131C_R0.FITS)
     Frequency Spectral_Value         Stokes BrightnessUnit       BeamArea 
  4.8851e+09Hz      4.8851GHz              I        Jy/beam        19.2738 
          Npts            Sum    FluxDensity           Mean            Rms 
           456   3.009941e-02   1.561674e-03   6.600747e-05   1.049135e-04 
       Std dev        Minimum        Maximum   region count 
  8.163621e-05  -2.127903e-05   4.782451e-04              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----