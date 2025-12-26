# SCOP install


```{r}

pak::pak("mengxu98/scop")


scop::PrepareEnv(
    envname = "scop_py",
  conda = "/soft/sim/mambaforge/bin/conda")

pak::pkg_install("SingleCellExperiment")
```

```{sh}

pak::pkg_install("plger/scDblFinder")

pak::pkg_install("proxyC")

```



## reference

+ [scop](https://github.com/mengxu98/scop)

