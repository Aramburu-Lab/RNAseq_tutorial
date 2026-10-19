#!/usr/bin/env bash
# F2: what the figures image holds; run inside it:
#   docker run --rm --platform linux/amd64 -v "$PWD":/w --entrypoint bash ctx-figures:f2 /w/docs/plans/logs/2026-10-19_f2_image_contents.sh
env=/app/.pixi/envs/laptop
echo "conda packages: $(ls "$env"/conda-meta/*.json | wc -l)"
cd "$env/lib/python3.12/site-packages" || exit 1
echo "python distributions: $(ls -d ./*.dist-info | wc -l)"
echo "largest site-packages entries, MB:"
du -sm ./* 2>/dev/null | sort -rn | head -15
echo "largest shared libraries in lib/, MB:"
du -sm "$env"/lib/*.so* 2>/dev/null | sort -rn | head -8
