
# 配置单细胞分析python环境

## 安装软件

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

## 参考


[scvi-tools](https://docs.scvi-tools.org/en/latest/installation.html)
[scib](https://scib.readthedocs.io/en/latest/)
