Not worth making maps: 
- because the SNR is poor (--> averaging might just be using those high phase noise channels) s.t. a map might not be very representative

Note: 
- checked flux scale (consistent 3C279 @ 1.5mm with online database) for both datasets
- phase calibrator is faint and poor weather (for set 2)
  + can try different solint..., but calibrator is weak, can be too short, also poor weather --> can't be too long
- phase calibrator is far from source (for set3)
  + i.e. gain solutions aren't perfect for target, as they are derived for a different direction




### Directories
imagingD23/
+ contains data products combining D2 and D3 sets

reduced_pretty_D2/
+ contains imaging files for D2 dataset

reduced_pretty/
+ contains imaging files for D3 dataset

reduced_testflag/
+ contains calibrated miriad files ready for imaging
+ *_wide.mir

tmp_D?config_testflag/
+ calibration by-products

raw/
+ raw data from DR
