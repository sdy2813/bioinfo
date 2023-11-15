# <p align=center> Fisher's Exact Test </p>

## Fisher's Exact Test 主要类型和适用场景

_Fisher's Exact Test_ 是一种用于分析两个分类变量之间是否相关的统计检验方法，特别适用于样本量较小的情况。它的主要用途是确定两个比例或发生率是否存在显著差异。此检验通常用于 2x2 的列联表。Fisher's Exact Test 分为单侧（one-sided）和双侧（two-sided）两种形式，它们的主要区别和适用场景如下：

### Two-sided Fisher's Exact Test（双侧检验）

**目的**：检验两个分类变量之间是否存在任何形式的关联（不论方向性）。

**应用场景**：当你不确定两个变量之间的关系是增加还是减少，或者当你对两种可能的关系都感兴趣时使用。例如，如果你想知道某种治疗是否影响疾病恢复的几率，但事先没有具体的假设是提高还是降低。

**计算方式**：计算观察到的数据以及比观察到的数据更极端的数据出现的概率之和。

### One-sided Fisher's Exact Test（单侧检验）

**目的**：检验两个分类变量之间是否存在特定方向的关联。

**应用场景**：当你有具体的假设预测两个变量之间的关系时使用。例如，如果你有理由相信某种治疗会提高疾病恢复的几率，那么你会使用单侧检验来仅检验这个方向的关系。

**计算方式**：只计算观察到的数据以及比观察到的数据在特定方向上更极端的数据出现的概率。

总结来说，选择单侧还是双侧检验取决于研究的假设和目的。如果你有具体的方向性假设，那么单侧检验是合适的；如果你只是想知道两个变量是否相关，而不关心具体的方向性，那么双侧检验是更好的选择。

## 示例

让我们来看两个例子，一个适合使用单侧 Fisher's Exact Test，另一个适合使用双侧 Fisher's Exact Test，并且展示如何使用 R 语言来执行这些检验。

例子 1: 单侧 Fisher's Exact Test

场景：假设你正在研究一种新药对于治疗某种疾病的有效性。你的假设是这种新药会提高治愈率。

数据：在一个小规模的临床试验中，10位接受新药治疗的患者中有7人治愈，而在10位接受安慰剂的患者中只有2人治愈。

列联表：
|     | 治愈 | 无效 | 总计 |
| --- | -- | -- | -- |
| 新药  | 7  | 3  | 10 |
| 安慰剂 | 2  | 8  | 10 |
| 总计  | 9  | 11 | 20 |

```Rs
# 创建列联表
data <- matrix(c(7, 3, 2, 8), nrow = 2, byrow = TRUE)
rownames(data) <- c("新药", "安慰剂")
colnames(data) <- c("治愈", "未治愈")

# 执行单侧 Fisher's Exact Test
fisher.test(data, alternative = "greater")

	Fisher's Exact Test for Count Data

data:  data
p-value = 0.03489
alternative hypothesis: true odds ratio is greater than 1
95 percent confidence interval:
 1.155327      Inf
sample estimates:
odds ratio 
  8.153063
```

例子 2: 双侧 Fisher's Exact Test

场景：假设你正在研究男性和女性在某种特定决策上的差异，而你没有特定的预期哪个性别更倾向于某个决策。

数据：在一个调查中，15名男性中有10人选择了方案A，5人选择了方案B。在20名女性中，有5人选择了方案A，15人选择了方案B。

列联表：
| |    方案A | 方案B | 总计 |
| --- | --- | -- | -- |
| 男性  | 10  | 5  | 15 |
| 女性  | 5   | 15 | 20 |
| 总计  | 15  | 20 | 35 |

```Rs
# 创建列联表
data <- matrix(c(10, 5, 5, 15), nrow = 2, byrow = TRUE)
rownames(data) <- c("男性", "女性")
colnames(data) <- c("方案A", "方案B")

# 执行双侧 Fisher's Exact Test
fisher.test(data)

	Fisher's Exact Test for Count Data

data:  data
p-value = 0.01923
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  1.12163 33.98348
sample estimates:
odds ratio 
  5.657108
```

在这两个例子中，R 代码创建了一个 2x2 的列联表并使用 fisher.test 函数来进行检验。对于单侧检验，我们使用了 alternative = "greater" 参数来指定我们的假设方向。在双侧检验中，我们没有指定 alternative 参数，因为默认情况下 fisher.test 执行的是双侧检验。

## 结论

当 Fisher's Exact Test 显示统计显著性时，我们可以从结果中得出一些关键的结论。但在具体说明之前，重要的是要理解“统计显著性”意味着什么：在统计显著的情况下，我们观察到的数据模式出现的概率，在无效假设（即两个变量无关联）为真的情况下非常低。

结论 1: 单侧 Fisher's Exact Test

在第一个例子中，单侧 Fisher's Exact Test 用于测试新药相对于安慰剂是否显著提高了治愈率。

如果结果统计显著，这意味着在假设新药和安慰剂效果相同的前提下，观察到的或更极端的结果（即新药治愈率显著高于安慰剂）出现的概率非常低。

结论：可以合理推断新药治疗比安慰剂效果更好，至少在统计上是如此。

结论 2: 双侧 Fisher's Exact Test

在第二个例子中，双侧 Fisher's Exact Test 用于测试男性和女性在选择方案A和方案B上是否存在显著差异。

如果结果统计显著，这意味着在假设性别和方案选择无关的前提下，观察到的或更极端的结果出现的概率非常低。

结论：可以合理推断男性和女性在这项决策上存在显著差异。

## 注意事项
因果关系：虽然 Fisher's Exact Test 可以告诉我们两个变量之间的关联是显著的，但它并不能证明因果关系。例如，在第一个例子中，虽然新药与更高的治愈率相关，但这并不一定意味着新药导致治愈率提高。
样本大小和效应大小：统计显著性并不等同于实际（或临床）显著性。一个结果可以在统计上显著，但实际上的效应大小可能很小，因此在解释结果时需要考虑这一点。
其他因素：应当考虑可能影响结果的其他因素，如样本偏差、实验设计问题等。
总之，统计显著性是一个重要的指标，但应当结合实验设计、效应大小和其他相关因素来全面理解和解释结果。


