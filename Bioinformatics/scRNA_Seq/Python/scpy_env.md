

```{sh}

mamba create -n scpy scanpy python=3 r-base r-seurat

mamba activate scpy
pip install 'scib[bbknn]'
mamba install scanorama
pip install -U "scvi-tools[cuda]"
pip install --quiet scvi-colab
mamba install ipython
pip install memento-de

```


