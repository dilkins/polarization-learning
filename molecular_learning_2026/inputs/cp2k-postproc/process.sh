#!/bin/bash

RUNDIR=$(dirname $(readlink -f "$0"))

${RUNDIR}/associate-wannier-centres.sh
${RUNDIR}/get-results.sh
