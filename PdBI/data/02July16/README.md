want to make a de-lensed spectrum.

For the central 800-1000km/s
1: plot the spectrum binned to the velocity reoslution of lens model
2: overplot with full-res spectrum; plot noise as faint color
3: de-lensed
4: mark magnitication factors
--> see 
    code: de-lensedComplicated.map
    figure: de-lensedSpecCO21_complicated.eps

pre-requisite Files for plotting:

copy RXJ1131 delensed from 22May16 (22May16/delensed.spec --> ./de-lensed_RXJ1131.spec)
- keep only chan 126-160
- pad one each at begin and end
--> bin by 5 (de-lensed_RXJ1131_binned.spec)

copy original full-res from 22May16 (22May16/sup127_155_2ndcln_noCont.spec -->./full-res_centralChans.spec)
- keep only chan 126-160
- pad one each at begin and end

copy full_res_centralChans.spec --> binned_centralChans.spec
- bin by 5
- pad one each at begin and end

copy ../22May16/de-lensed_companion.spec to this dir
- keep only chan 126-160
- pad one each at begin and end

- make noise.spec for line-free channels, original full-res
--> not actually needed