#!/bin/bash

source $1

Rscript ${ScriptDir}/SaveResults.R $OutPrefix $CisPvalue $TransPvalue $Dist $TransCrossChr $SavecsvCis $SavecsvTrans

