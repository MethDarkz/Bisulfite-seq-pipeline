#!/bin/bash -e
gunzip -c "$1" | awk '{if (match($0,"^chr")) {chr=$0} else print chr","$0}' | awk 'BEGIN {FS=","} {printf "%s,%i,%i,%i,%i,%i,%i\n", $1,$2,$3,$4,$5,$6,$7}'