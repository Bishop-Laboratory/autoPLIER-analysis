# Cancer - Normal Transformer

Repo structure:

```

```


## TLDR

TLDR: This repo is for a research paper that develops 3 approaches for addressing one question: "How do cancers relate to normal tissues?" In another sense, this question is asking:

1. What are the normal cell states that a tissue can inhabit?
2. What are the normal states which cancer most resembles?
3. How do interventions, such as chemotherapy, change this relationship?
4. What new drug targets does this reveal? Does this tell us anything about the cell of origin?

Here are the three approaches:

1. PHATE atlas: We use manifold learning approaches to embed normal tissues and cancers in a low-dimensional space and then we can see how similar or different the cancer is from various normal tissues (e.g., PHATE, H-PHATE, clustering, etc). **Pro**: Very easy to visualize and interpret. **Con**: It's not that accurate as dimensionality reduction causes lots of information loss.
2. CellO Plier: We use NMF to find pathway-based latent dimensions that we can use to map cancer samples to the cell ontology and derive normal tissue scores for each cancer sample. E.g., Cancer sample X is 50% T-Cell and 50% Stem Cell. 
3. Cancer scores for normal tissues: This Dr. Zheng's idea (this is a tertiary approach we might not have time to get to) -- rather than finding "normal tissue scores in cancer" he wants to find "cancer scores in normal tissues" in order to find the normal tissues which most resemble a particular cancer. This approach is also optimized to deal with an issue relating to the number of genes expressed influencing the mapping of cancer to normal. 

## Background

This repository is based on a simple idea:

There is are a discrete number of **Normal** cellular states which are defined by the expression of specific **genes** in specific **pathways** -- and that we can learn these states using RNA-Sequencing. One popular depiction of this idea is the *Tabula Muris* - [article](https://www.nature.com/articles/s41586-018-0590-4):

![alt](https://ds.czbiohub.org/images/Introduction/facs_tsne_by_tissue.png)

From our recent Cancers [paper](https://www.mdpi.com/2072-6694/12/4/948), we showed that at least one cancer (Ewing sarcoma) has a cellular state which is similar to the cellular states found during normal embryonic development. More improtant, we learned how certain interventions change the state of Ewing sarcoma tumor cells. This led us to realize new drug targets which were subsequently [demonstrated](https://www.nature.com/articles/s41389-020-00294-8) (albeit not by us...). However, the method we used (PHATE) is more qualitative with respect to answer questions of the normal-cancer relationship. For a larger-scale implementation of this approach, we need more quantitative tools in addition to qualitative ones. 

## The approach

### A Cancer PHATE Atlas

We plan to extend the method we showed previously to all other cancers. Then, we hope to make a web application where users can explore their cancer of interest and see how it relates to normal tissues. 

### Cancers scored with normal tissue levels

Using a method similar to PLIER, we can incorporate pathway-level information to score cancer samples based on their expression of pathways which are related to particular normal tissues. For example, we might have found that Ewing sarcoma has a score of 50% Mesodermal and 50% pluripotency. How does an intervention change these scores? This would give us an innovative and quantitative way to assess the same questions posed above. 

### Normal tissues scores based on cancer signatures

This Dr. Zheng's idea which may provide an interesting supplement to our other analyses -- but he hasn't yet disclosed the details of his approach to me. 


