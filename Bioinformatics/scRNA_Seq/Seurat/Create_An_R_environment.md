
# How to create A useful Conda R virtual Environment

当你已经安装过conda或者mamba之后，你就可以用它来创建属于你的虚拟环境，R环境自然不在话下。

对于单细胞分析，R 语言包Seurat 是最为常用的，且R 擅长绘图和统计分析。在创建 R虚拟环境的时候可以一并将我们想要的软件也装上，避免软件版本之间发生冲突。

## Create R Environment

```{sh}
mamba create -n Seurat4  r-base r-ggsci r-ggseqlogo r-ggsignif r-ggthemes r-ggupset r-ggrepel r-ggraph r-igraph r-patchwork r-pak r-reshape2 r-rmarkdown r-seurat=4.4.0 r-seuratobject=5.0.1  r-vegan r-tidyverse python=3.7.12 umap-learn=0.4.6 prompt-toolkit=3.0.39 numba=0.51.2 numpy=1.21.6 radian
```

## install R packages

pak 是一个很棒的用于安装R包的工具，可以先安装上，后续的需要安装R包就可以用pak来安装了。

```{R}
install.packages("pak")
```

## Code Server R

如果一个服务器不支持使用Rstudio-Server，那么Code Server 将是一个不错的替代品。Code Server 的安装另有教程，不再赘述。

安装两个R语言相关的插件：R 和 R Debugger

```{R}
pak::pkg_install("httpgd")
pak::pkg_install("languageserver")
```

接着我们在Code Server 的设置中输入 r.plot.useHttpgd，启用httpgd。

![](https://files.mdnice.com/pic/1df43cab-f266-47de-b48c-0bfee4531869.png)





