#!/bin/bash
set -e

export BLASTDB="/usr/local/ncbi-blast-2.6.0+/db/"
export PATH="/usr/local/ncbi-blast-2.6.0+/bin:$PATH"

cd /data/

if [ "$1" = 'update' ]; then
    shift;
    exec /refineRDPdatabase.py $@
fi

exec /identification16S.py $1 16S_named_RDP
