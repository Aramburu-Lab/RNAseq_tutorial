cd "$CTX_REPO"
pixi install -e laptop --frozen
cp -n conf/dataset_template_gdc.toml conf/datasets/tcga_hnsc.toml
