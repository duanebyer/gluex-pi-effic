#!/bin/tcsh
# Helium target runs only
python launch.py config-hdroot/RunPeriod-2021-11.config 90001 90206 -v True
python launch.py config-hdroot/RunPeriod-2021-11.config 90603 90662 -v True

