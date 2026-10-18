cd "$CTX_REPO" && pixi run -e laptop ctx --root ~/ctx_tcga_hnsc setup --dataset-file conf/datasets/tcga_hnsc.toml
pixi run -e laptop ctx --root ~/ctx_tcga_hnsc run --dataset tcga_hnsc --dataset-file conf/datasets/tcga_hnsc.toml --gdc-dir ~/ctx_tcga_hnsc/data/raw/tcga_hnsc --mirna-dir ~/ctx_tcga_hnsc/data/raw/tcga_hnsc_mirna --isoform-dir ~/ctx_tcga_hnsc/data/raw/tcga_hnsc_isoform --check
