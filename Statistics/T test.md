# <p align="center"> T test </p>

T-test，即t检验，是一种统计方法，用于判断两组数据的均值是否存在显著性差异。t检验的几种主要类型包括单样本t检验、独立样本t检验（又分为两种情况），以及配对样本t检验。

1. **单样本t检验（One-Sample T-Test）**
   - **原理**：比较一组样本的均值与一个已知或假设的总体均值之间是否存在显著差异。
   - **适用场景**：当你有一个样本数据集和一个理论上的或已知的总体均值时，例如，检验一个班级的平均考试成绩是否显著高于全国平均水平。

2. **独立样本t检验（Independent Two-Sample T-Test）**
   - **原理**：比较两个独立样本的均值之间是否存在显著差异。
   - **适用场景**：适用于两组相互独立的样本数据比较，例如，比较两个不同学校学生的考试成绩。
   - **注意**：这种检验又分为两种情况：
     - **方差齐性**（假设两组数据的方差相等）。
     - **方差不齐**（假设两组数据的方差不等）。

3. **配对样本t检验（Paired Sample T-Test）**
   - **原理**：比较来自同一组受试者在两个不同条件或时间点的均值之间是否存在显著差异。
   - **适用场景**：当你有成对的数据或者重复测量的数据时，例如，检验一个人在接受治疗前后的血压变化。

在实际应用中，选择哪种类型的t检验取决于数据的特点和研究问题。在使用时，还需关注数据是否符合正态分布和方差齐性等假设条件。在R语言中，可以使用内置的`t.test()`函数来执行这些检验，根据不同的需要设置不同的参数。例如，对于独立样本t检验，可以设置参数`var.equal = TRUE`或`FALSE`以适应方差齐性或不齐性的情况。

## Paired T-Test
配对样本t检验（Paired Sample T-Test）的数学基础建立在以下几个关键点：

1. **差异分数（Difference Scores）**：计算每一对样本间的差异。如果有n对样本，则会有n个差异分数。

2. **均值和标准偏差**：计算这些差异分数的均值（\($\bar{d}\$)）和标准偏差（\($s_d\$)）。

3. **t统计量**：计算t值，其公式为：
   $\[ t = \frac{\bar{d}}{s_d / \sqrt{n}} \]$
   其中，\(\bar{d}\)是差异分数的均值，\(s_d\)是标准偏差，n是样本对数。

4. **自由度（Degrees of Freedom, df）**：配对t检验的自由度是\(n-1\)，其中n是样本对数。

5. **p值**：使用t值和对应的自由度，通过t分布表或计算软件来找到p值，以判断结果的显著性。

### R语言实例

假设我们有一个简单的数据集，其中包含10位患者在治疗前后的某项指标测量值：

```R
# 创建数据
pre_treatment <- c(8, 7, 6, 9, 10, 5, 7, 8, 9, 6)
post_treatment <- c(5, 6, 5, 6, 5, 4, 5, 6, 6, 5)

# 执行配对样本t检验
t.test(pre_treatment, post_treatment, paired = TRUE)
```

在这个例子中，`pre_treatment`和`post_treatment`分别是治疗前后的测量值，我们使用`paired = TRUE`参数来指定这是一个配对样本t检验。

### 图形解释

我们可以通过绘制差异分数的条形图来直观展示这些差异。这样的图形有助于理解配对样本之间的变化情况：

```R
# 计算差异
differences <- pre_treatment - post_treatment

# 绘制条形图
barplot(differences, main = "Differences in Measurements", 
        xlab = "Patient", ylab = "Difference in Measurement")
```

这个条形图将展示每位患者治疗前后测量值的差异。正值表示治疗前的测量值高于治疗后，负值则相反。这样的视觉展示可以帮助理解数据的变化趋势和差异的一般情况。


## Paired T-Test p值计算

### 推导过程

1. **计算差异分数**：对于每一对样本（如前后治疗的数据），计算它们的差异（\($d_i = X_{i1} - X_{i2}\$)）。

2. **计算差异分数的均值和标准差**：
   - 均值：\($\bar{d} = \frac{\sum{d_i}}{n}\$)
   - 标准差：\($s_d = \sqrt{\frac{\sum{(d_i - \bar{d})^2}}{n - 1}}\$)

3. **计算t统计量**：
   $\[ t = \frac{\bar{d}}{s_d / \sqrt{n}} \]$
   其中，n是样本对的数量。

4. **计算p值**：使用t值和自由度（df = n - 1），在t分布表中找到相应的p值，或者使用统计软件进行计算。

### R代码解释

假设我们有一组前后治疗的数据：

```R
# 创建数据
pre_treatment <- c(8, 7, 6, 9, 10, 5, 7, 8, 9, 6)
post_treatment <- c(5, 6, 5, 6, 5, 4, 5, 6, 6, 5)

# 计算差异
differences <- pre_treatment - post_treatment

# 计算差异的均值和标准差
mean_diff <- mean(differences)
sd_diff <- sd(differences)
n <- length(differences)

# 计算t值
t_value <- mean_diff / (sd_diff / sqrt(n))

# 计算p值
p_value <- 2 * pt(-abs(t_value), df = n - 1)
```

在这个R代码中：

- `differences` 是前后治疗数据的差异。
- `mean_diff` 和 `sd_diff` 分别是这些差异的均值和标准差。
- `t_value` 是基于这些差异计算的t统计量。
- `p_value` 是通过`t_value`和自由度计算出的p值，使用`pt`函数来获取t分布的累积分布函数值。因为我们进行的是双尾检验，所以乘以2。

最后，`p_value`就是我们的结果，表示治疗前后数据均值差异的显著性。如果p值小于显著性水平（如0.05），则认为治疗前后存在显著差异。












