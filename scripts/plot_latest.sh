#!/usr/bin/env bash
set -euo pipefail

if [ ! -d "results" ]; then
  echo "No results/ directory found. Run the app first."
  exit 1
fi

gnuplot -persist scripts/plot_latest.gp