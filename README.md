# HPIES
HPIES processing code


This repo includes code from John Dunlap to process HPIES data.  Cronjob examples as follows:

45 * * * * /data/rsn/hpies.run >> /data/rsn/hpies.log 2>&1
17 22 * * * /data/magobs/getmagobs.py new tuc hon frn >> /data/magobs/getmagobs.log 2>&1
