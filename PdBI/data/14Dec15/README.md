To figure out the total observing time after flagging.

- Don't know how to do this in GILDAS, but MIRIAD `uvindex` should do it


1. copied co21_trial1-noCont.uvt from /Users/admin/Research/RXJ1131/PdBI/data/04Sep15, where this uvt is after removing continuum.
    + in principle, shouldn't matter  
2. export the uvt into uvfits in GILDAS
3. convert uvfits into .mir inside MIRIAD
4. uvindex vis=co21_trial1-noCont.mir log=log=co21-noCont_totalObsTime.log
5. output: Total observing time is  3.74 hours
    + consistent with PdBI quality log file, where cumulative observing time: 3.75 hours  
