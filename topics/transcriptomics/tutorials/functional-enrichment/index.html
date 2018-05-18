---
layout: tutorial_hands_on
topic_name: training-test
tutorial_name: new-tutorial-test
---

# Introduction
{:.no_toc}
---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial
<!-- This is a comment. -->

In this tutorial you will learn about on how work and analyse functional enrichment data.
Analyse gene expression data the enrichment analysis method is based on the functional annotation of the differentially expressed genes. The enrichment analysis consists of the application of statistical testes in order to verify if the sample set of entities is enriched in comparison with the whole population. When says enrichment, means that the frequency of the sample will be higher than would be expected by chance given the population frequency.
> 
> 
> 
What is the Gene Ontology?
> 
The [Gene Ontology](http://www.geneontology.org "Gene Ontology Homepage") (GO) gives an structured, controlled vocabularies and classifications of many domains of molecular and cellular biology. Each GO type apresents an structure of directed acyclic graph (a hierarchy with multi-pareting). Contain three types of independent ontologies (GO terms):  biological process (e.g., signal transduction), molecular function (e.g., catalytic activity) and cellular component (e.g., ribossome). 

![QuickGO - http://www.ebi.ac.uk/QuickGO](/home/jpaulas/Pictures/Gene_ontology_example.png)
> 
What is GO annotation?
> 
The genes are relacioned with GO terms via annotations, and is possible to have multiple annotations associated. According to the true path rule, a gene annotated to a term is implicity annotated to each ancestor of that term. The GO annotations have evidence codes that encode the type of evidence supporting them.    

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. Functional Enrichment Analysis
> 2. 
> {:toc}
>
{: .agenda}

> 
> 
# Functional Enrichment Analysis

How work the functional enrichment analysis? When we have an RNA-seq experiment, a functional classification shema and the corresponding functional annotations, the analysis have a certain order: is necessary to determine what should be the sample and population sets of genes; and then, the program will compute all inferred annotations, compute n and N (total count of genes in the sample/population); for each functional annotation get the counts k and K (count of genes in the sample/population) and for last step, compute the statistical test.

The appropriate statistical test is the one-tailed variant of Fisher’s exact test, also known as hypergeometric test for over-representation. When the one-tailed version is applies, this test will compute the probability of observing at least the sample frequency, giben the population frequency. 
The hypergeometric distribution consistes in the probability of k successes in n draws, without  replacement, from a finite population of size N that contains exactly K successful objects:  (FALTA A FORMULA)

> ### {% icon question %} Question
> 
> Question?
>Shouldn’t we correct for multiple testing?
> <details>
> <summary>Click to view answers</summary>
> In general should, in cases we’re testing multiple functional aspects: the statistical testing is based on the probability of erroneous rejection of the null hypothesis being low; but if you make an multiple related tests, the probability of at least one of them being a false positive increases. 
However, the stochastic event (sampling of genes), normally it was already been the subject of statistical testing and multiple test correction. 
> </details>
{: .question}

<!--
{% icon hands_on %} will render the hands_on icon as specified in
_config.yml in the root of this repository.
-->

For the first exercice we go to use the results of differential expression Drosophila melanogaster from the article of Trapnell et al. (link for article). 

> ### <i class⁼"fa fa-pencil" aria-hidden="true"></i> Hands-on: 
>
> 1. **Create a new history** 
> 2. **Upload to the Galaxy** the following files:
> ..* go.obo  (FALTA O LINK)
> ..* drosophila_gene_association.fb (FALTA O LINK)
> ..* Trapnell_diff_genes.txt (FALTA O LINK)
> ..* Trapnell_Diff_Table.tab (FALTA O LINK)
> 
>    > ### <i class="fa fa-lightbulb-o" aria-hidden="true"></i> Tip: Upload data to Galaxy [1](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-peaks2genes/tutorial.html)
>    >![](/home/jpaulas/Pictures/galaxy_upload.png)
>    > * **Click** on the upload button in the upper left of the interface.
>    > 
>    > * Press **Choose local file** and search for your file.
>    > * Press **Start** and wait for the upload to finish. Galaxy will automatically unpack the file.
>    {: .tip}
> **Rename** the *go.obo* file to <span style="color:red">some **GO** text</span>, *drosophila_gene_association.fb* file to <span style="color:red">some **GO annotations Drosophila melanogaster** text</span>, *Trapnell_diff_genes.txt* file to <span style="color:red">some **Trapnell_study_sample** text</span> and *Trapnell_Diff_Table.tab* file to <span style="color:red">some **Trapnell_population** text</span>.
> ...After you upload the files, and if you press the eye icon you should look someting like this:
> ![](/home/jpaulas/Pictures/trapnell_file.png)
> **Figure 1** Trapnell file
> 
> ...The both files have the same information, the little difference between the two files is the number of genes.
> 
>    > ### <i class="fa fa-commenting-o" aria-hidden="true"></i> Comments
>    > To create the study sample file was necessary to select the best genes, and for that, you need to regulated by the column of p-value, and choose the level of significance do you want. In this case, i select the genes with a p-value < 0,05. 
>    {: .comment} 
> 
> 5. **GOEnrichment** <i class="fa fa-wrench" aria-hidden="true"></i>: Run the GOEnrichment tool with the four files.
> ..* Select “No” in the option Summarize Output. 
> ![](/home/jpaulas/Pictures/galaxy_trapnell.png)
> 6. This will generate 6 files with the respective default names:  MF_Result.txt, BP_Result.txt, CC_Result.txt, MF_Graph, BP_Graph and CC_Graph. Rename MF_Trapnell, BP_Trapnell, CC_Trapnell, MF_Graph_Trapnell, BP_Graph_Trapnell and CC_Graph_Trapnell. 
> ... If you press the eye icon of the three graphs you should look someting like this:
> ..* <span style="color:blue">some **Biological Process** text</span>
> ![](/home/jpaulas/Pictures/BP_trapnell.png)
> ..* <span style="color:blue">some **Molecular Function** text</span>
> ![](/home/jpaulas/Pictures/MF_trapnell.png)
> ..* <span style="color:blue">some **Cellular Component** text</span>
> ![](/home/jpaulas/Pictures/CC_trapnell.png)
> 
>    > ### {% icon question %} Question
>    >
>    > Question?
>    > How this method apply to RNA-seq gene sets? In the sample and population?
>    > <details>
>    > <summary>Click to view answers</summary>
>    > When we apply in the sample, the set of differentially or over- and under-expressed genes, it always go to depent on the biological question being addressed. And the case of the population, is the transcriptome, i.e., all genes existing in the RNA-seq experiment with meaningful counts. 
>    > </details>
>    {: .question}
{: .hands_on}
> 
# Interpretation of the results
The interpretation we make of the results will depend on the information we want to retrieve and the motivation to do the analysis. During the analysis, we must not forget that the results that are statistically significant are different from the biologically meaningful. However, there may be situations where the statistically enriched terms provide some biological or technical insight into the underlying experiment, even if not immediately apparent.
Terms that are very generic tend to be difficult to interpret, since the very specific terms are generally not integrative. It is intended that terms that are sufficiently specific to impart substantial biological, however, that are generic enough to integrate multiple genes.

>    > ### {% icon question %} Question
>    >
>    > Question?
>    >When can be applie the enrichment analysis?
>    > <details>
>    > <summary>Click to view answers</summary>
>    > Enrichment analysis can be used in validation (e.g., of a protocol for extracting membrane proteins), characterization (e.g., of the effects of a stress in a organism) and elucidation (e.g., of the functions impacted by the knock-out of a transcription factor). 
>    > </details>
>    {: .question}

This produces an output for each GO type and the interpretation varies. The output results a table with the p-values and frequencies. Also gives, based of the semantics the GO terms, an graph. From the graphs is possible to visualize the enrichment results and highlight enriched ontology branches. But sometimes the graph can become overwhelmed because the size of the results.

Now we go to use differential expression results from an study of a mouse. 

> ### <i class⁼"fa fa-pencil" aria-hidden="true"></i> Hands-on:
>
> 1. **Upload to the Galaxy** the mouse_brain_vs._heart.txt (LINK) and the Mus_musculus_annotations_biomart_e92.tab (FILE) files.
> 2. **Rename** the *mouse_brain_vs._heart.txt* file to <span style="color:red">some **Mouse population** text</span> and  *Mus_musculus_annotations_biomart_e92.tab* file to <span style="color:red">some **GO annotations _Mus musculus_** text</span>.
> ... After you upload and if you press the eye icon you should look someting like this:
> ![](/home/jpaulas/Pictures/mouse_file.png)
> 3. **Upload to the Galaxy** the two study samples files, one with overexpressed genes (*mouse_up.txt* LINK) and the other with underexpressed genes (*mouse_down.txt* LINK).
>    > ### {% icon question %} Question
>    >
>    > Question?
>    >How do you know which genes are up- and downregulated?
>    > <details>
>    > <summary>Click to view answers</summary>
>    > It is through the logFC values that we derive the information whether the gene is up- or downregulated. If the logFC value is positive it means that the gene is upregulated, and if it is negative the gene is downregulated.
>    > </details>
>    {: .question}
> 4. **Rename** the *mouse_up.txt* file to <span style="color:red">some **Mouse overexpressed** text</span> and the mouse_down.txt file to <span style="color:red">some **Mouse underexpressed** text</span>.
> 5. **GOEnrichment** <i class="fa fa-wrench" aria-hidden="true"></i>: Run GOEnrichment for the both files **Mouse overexpressed** and **Mouse underexpressed**, with the **GO annotations _Mus musculus_** file you downloaded earlier.
> ..* Select “No” in the option Summarize Output. 
>    > ### {% icon comment %} Comments
>    > A comment
>    {: .comment}
>
>    > ### {% icon tip %}Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}




# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
