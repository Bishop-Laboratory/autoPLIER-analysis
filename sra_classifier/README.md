# SRA Classifier 

This should be a software package which is used to classify samples from SRA.

Ideal usage:

As a python module:

```
import sra_classifier
import pandas as pd

meta = pd.read_csv("example_metadata.csv")
expr = pd.read_csv("example_gene_counts.csv")

classifier = sra_classifier.SRAClassifier()
predictions = classifier.fit_transform(meta = meta, expr = expr)

```

As an R package:

```
library(sra_classifier)

meta <- read.csv("example_metadata.csv")
expr <- read.csv("example_gene_counts.csv")

classifierfit <- sra_classifier(meta = meta, expr = expr)
predictions <- predict(classifierfit, list(meta = meta, expr = expr))

```


As a CLI tool:

```commandline
sra_classifier -m example_metadata.csv -e example_gene_counts.csv -o predictions.csv
```


