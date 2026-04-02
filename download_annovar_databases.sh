#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="/export/home/public_test1/annovar/annovar/annovar-databse"

_d="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANNOVAR_PL="${ANNOVAR_PL:-$_d/annotate_variation.pl}"


mkdir -p "$OUT_DIR"

run() {
  echo "+ $*"
  "$@"
}

echo "Downloading ANNOVAR databases into: $OUT_DIR"

run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar 1000g2014oct "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar avsnp142 "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar clinvar_20150330 "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar cosmic70 "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar caddindel "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar cosmic74 "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar gnomad_exome "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb gwasCatalog "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg19 -downdb -webfrom annovar ljb26_all "$OUT_DIR"

run "$ANNOVAR_PL" -buildver hg38 -downdb -webfrom annovar avsnp150 "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg38 -downdb -webfrom annovar clinvar_20190305 "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg38 -downdb -webfrom annovar gnomad_exome "$OUT_DIR"
run "$ANNOVAR_PL" -buildver hg38 -downdb -webfrom annovar ljb26_all "$OUT_DIR"
