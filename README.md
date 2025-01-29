# HPIES
HPIES processing code


John Dunlap, the instrument engineer and subject matter expert (SME) at UW, has provided code that can convert the L0 data from the uFrame system to create the L1 and L2 data products. This is the code that was originally accessible through an APL server (https://henry.apl.uw.edu/~dunlap/hpies-code/).    

Cronjob examples running to process data automatically:

45 * * * * /data/rsn/hpies.run >> /data/rsn/hpies.log 2>&1
17 22 * * * /data/magobs/getmagobs.py new tuc hon frn >> /data/magobs/getmagobs.log 2>&1
