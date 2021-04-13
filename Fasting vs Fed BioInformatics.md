

**Background**

Fasting is a practice as old as humanity, yet little is known about the biological implications of its

implementation. Ancient accounts of fasting being used as a medical intervention begins with

Hippocrates’ observation that fasting was the only effective remedy for epilepsy [1]. Since the

5th century, the application of fasting has expanded vastly. Fasting has been shown to increase

the mean lifespan of multiple model species including c.elegans and mice by 50% and 15%,

respectively [2,3]. Not only has fasting been shown to increase lifespan, but also have positive

effects on certain pathologies that plague us to this day. Neurodegenerative diseases such as

Alzheimer’s Disease and Parkinson’s Disease are also beneficially affected by fasting [4]. When

applied properly, fasting can have therapeutic effects on metabolic diseases such as Type 2

Diabetes [4]. It’s clear that fasting can have medicinal effects on multiple different conditions,

however it is imperative that more research is done to illuminate the metabolic processes that

underlie this phenomenon. Similar to all medical interventions, fasting has the potential to be

dangerous. Therefore, it is necessary that we understand the pathways and metabolic

consequences of fasting to navigate through obscure treatment intervention. With a better

understanding of fasting, not only will we have more methods to combat these diseases but we

will also have more potential drug targets for future drug development.

This study is motivated by a previous experiment done by Defour and colleagues [5]. The data

was downloaded from the NCBI GEO website. There were 11 participants in this study who had

perennial subcutaneous fat samples extracted from them at different points in time. The first

sample from each patient was obtained 2 hours after their last meal, while the second sample

was obtained 26 hours after their last meal. Subjects were all healthy individuals between the

ages of 40-70. The samples were then sequenced through Microarray Analysis using an

Affymetrix Chip.

**Hypothesis and Predictions**

In this study we are examining which genes are differentially expressed between a fed state and

a fasted state. Furthermore, we are looking to uncover the biological processes that are either

enhanced and/or perturbed due to being fasted for 26 hours compared to 2 hours. We anticipate

that there will be a plethora of genes that have their expressions modulated due to the fasting

regimen. We expect to see an enhancement in genes that regulate fatty acid oxidation,

autophagy, and protective pathways. When it comes to genes that have their expressions

repressed, we expect to see these genes regulate glucose metabolism, growth factors, and

insulin release.

**Methods**

To analyze the microarray data obtained from the study, the data had to be downloaded and

then loaded into R for further investigation. Due to inherent differences in base expression levels

and potential systematic biases [6], we must first normalize and filter the data. To normalize the

data, Bioconductor’s package *oligo* was used. Then the data was filtered using a soft median





filtering technique to remove probes that had low expression levels (shown in Figure 1),

recommended by the authors of the package limma [7]. Then genes were annotated using the

*hugene21transcriptcluster*, provided by Bioconductor [8]. An MA plot was used to initially

visualize the difference in log expressions of the samples that were kept. To identify genes that

were differentially expressed, a linear model was used with the help of Bioconductor’s package

*limma*. To determine whether or not a genes’ differential expression was significant or not,

values with a magnitude greater than 0.5 for log-fold change (푙표푔 퐹퐶) with a p-value less than

2

0.05 was used as a threshold. The threshold value for 푙표푔 퐹퐶is arbitrary, it is common practice

2

to use a threshold value of 1 (absolute value). However, due to the relatively short fasting

protocol that participants underwent, we don’t expect to see very large differences in gene

expression. Therefore, we will use a lower threshold value in determining whether or not a gene

is differentially expressed. Once genes were determined to be significantly differentially

expressed the remaining data was visualized with a volcano plot, with the help of Bioconductor’s

package *limma*. Finally, to extract the biological processes related to the genes that had their

expressions enhanced or repressed GO analysis was used. This was done using the Ma’ayan

Labs online software, *ENRICHR* [13,14,15].

**Results**

Upon analysis of the microarray data there was a significant amount of genes that were

differentially expressed. This was visualized with an MA plot (Figure 2), which shows the

relationship between two samples’ log expression ratio and their log mean ratio. Probes with a

푙표푔 퐹퐶 of greater than 0.5 were designated as differentially expressed. It is important to note

2

that an MA plot shows genes that are differentially expressed, without consideration of whether

it is significant or not. Once it was determined that there were in fact differentially expressed

genes, further statistical procedures were implemented to determine the statistical significance

behind the change in expression. By fitting a linear model to the obtained data we are able to

determine the variability in expression levels, thus allowing us to estimate the variability

between each sample and any background noise. Using an empirical Bayes Method, we are

then able to identify the standard errors between each sample’s log expression. Figure 3 shows

the relationship between a probe’s 푙표푔 퐹퐶and its computed p-value. Volcano plots are useful

2

because they are easy to understand and provide information on the experimental data as a

whole.

Using gene ontology analysis software, we were able to take genes that were significantly and

differentially expressed and predict the biological functions that are affected by the genes.

Taking the genes that were found to be upregulated, we find that many different pathways are

affected. There was an increase in processes that correspond to the regulation of glucose

metabolism, cell catabolism, and autophagosome maturation — to name a few. This is shown in

the Gene Ontology table pictured in Figure 4A. The top 5 genes that had their expressions

significantly increased the most include: LINC01474, DEPP1, PDK4, CFAP69, SLC19A3. There

was a decrease in processes that correspond to the regulation of lipid biosynthesis, protein





maturation, and glucose breakdown. This is shown in Figure 4B. The top 5 genes that had their

expressions decreased include: PNPLA3, ANGPTL8, SREBF1, SLED1, FFAR4.

These results help confirm some of our initial predictions that genes which would be

upregulated would have effects on glucose metabolism regulation, fatty acid oxidation, and

protective pathways (in response to damage signals). For our predictions about downregulated

genes, we saw a confirmation in the decrease of fatty acid biosynthetic pathways. However, we

did not see a decrease in genes corresponding to the release of insulin.

**Figures and Visualizations**

Figure 1) Applying a soft median filtering technique a histogram was created to identify a cutoff point where

probes that had a median intensity value of less than 2 were removed. This was done to remove

unnecessary points that do not have high expression levels which can lead to an increase in extreme

values of variance. Figure 2) An MA plot was created to visualize the differences in average log

expressions vs log ratio of fasted samples to fasted samples. Probes that had an absolute value of 0.5 or

greater for the log ratio of expressions was deemed as differentially expressed. Here we see that we have

a plethora of probes (genes) that fit the criteria of being differentially expressed.





Figure 3) A volcano plot illustrating the log fold change (fasted vs. fed) of genes that were filtered. A linear

model was used to characterize genes that were significantly differentially expressed. Genes that had their

expressions changed the most can be seen through a larger magnitude in log fold change values. The top

5 genes that had their expressions decreased include: PNPLA3, ANGPTL8, SREBF1, SLED1, FFAR4. The

top 5 genes that had their expressions increased include: LINC01474, DEPP1, PDK4, CFAP69, SLC19A3.

Figure 4A [left]) Gene Ontology analysis was done to obtain the biological processes that are impacted by

genes noted. On the left is a GO table that was created from genes that had their expressions increased.

Genes that had their expressions increased had an effect on processes that include cell catabolism,

glucose metabolism, and autophagosome synthesis (not shown).

Figure 4B [right]) GO table created from genes that had their expression levels increased significantly.

Many processes that were impacted negatively include lipid biosynthesis, protein maturation, and glucose

metabolism (not shown).





**Discussion**

The purpose of this study was to analyze which genes are differentially expressed in human

adipose tissue when a subject is fasted (26 hours since last meal) vs. fed (2 hours since last

meal). Furthermore, we wanted to infer the biological processes that were affected when

participants were subjected to 26 hours of fasting. It was found that processes that regulate

glucose metabolism, cell catabolism, and autophagosome maturation were upregulated in

fasted samples. On the other hand, it was found that processes that regulate lipid biosynthesis,

protein maturation, and glucose breakdown were downregulated in fasted samples. This helps

build our understanding behind the processes that are affected by fasting. This study does not

exhaustively identify all of the genes and processes that are affected by fasting.

The next step for researchers includes further experiments investigating the potential effects of

the differential expressions of these genes. Determining the physiological effects of these genes

can provide insight into whether or not fasting and/or its metabolic effects can be used in a

clinical setting. Investigators can also identify the metabolic pathways that these genes affect

and use these mechanisms to elucidate future dug leads.

One gene that had its gene expression downregulated is SREBF1, with a log-fold change value

of -2.1 and a p-value of 9 × 10−5. SREBF1 encodes for the SREBP-1(a/c) protein, which are

transcription factors that play a key role in glucose metabolism and lipogenesis [9]. SREBP-1c is

the dominant isoform in human liver and adipose tissue [9]. Because we are analyzing human

adipose tissue, it only makes sense that we try and understand the potential fasting plays on it.

This is thought to allow for constant cell proliferation by providing enough lipids for membrane

biogenesis. Therefore SREBF1 is essential for cell proliferation and survival. SREBF1’s role in

cell proliferation can even be analyzed from a cancer cell perspective. A major player in cancer

cell survival is an increase of de novo lipogenesis, which is needed for a growing cell

membrane. KO SREBF1 models have been shown to decrease cell proliferation in colon cancer

cells by decreasing the expression of its downstream targets such as: FAS (fatty acid synthase),

and HMGCR (HMG-CoA reductase) [10]. Because SREBP-1c is a pivotal factor in lipogenesis, it

would be interesting to investigate the efficacy of using forms of SREBF1 inhibition as another

weapon in our arsenal against cancer, imagine a combination of SREBF1 inhibition (fasting)

with standard cancer treatments.

To contrast, one gene that had its gene expression increased is ACER2, with a log-fold change

value of 1.5 and a p-value of 3 × 10−4. ACER2 is an alkaline ceramidase, which is responsible

for the conversion of ceramides into sphingosines. Sphingosines are lipids found on cell

membranes across all eukaryotic organisms and their generation produces byproducts such as

reactive oxygen species (ROS). ACER2 has been shown to induce autophagy and apoptosis

through the generation of ROS. p53 helps mediate the induction of apoptosis by acting as a

transcription factor for ACER2 [11]. Autophagy is a cellular process that plays an integral part of

cell recycling and pathology prevention. Dysfunction of autophagy has been linked to many

different diseases including cancer, neurodegenerative diseases, and metabolic diseases.

Interestingly enough, while ACER2 has been shown to induce autophagy and apoptosis, it also





has been proposed to play a positive role in Hepatocellular carcinoma (HCC) survival [12]. This

points out the complexity of cancer biology and leads to future questions that scientists may

investigate. Does the increase of ACER2 in HCC cells increase ROS concentration? If so, do

HCC cells have increased biomarkers for autophagy (increase in autophagosome formation or

increase in LCIII concentration)? Are HCC cells immune to autophagy and/or apoptosis?

This is just a stepping stone in understanding the complex effects of fasting on the human body.

This study was meant to provide a better understanding of the genes that are differentially

expressed in humans that are fasting vs. not fasting. But this does not remotely come close to

fully elucidating the mechanisms behind fasting. Many questions are left on the table, due to

many different factors such as experiment design. For example, this study was only done with

11 patients. This begs the question, was the sample representative of the population? Would we

get different results if we were to expand the patient cohort to a much larger number? Secondly,

this experiment only obtained samples of patients when they were fed and fasted. Investigating

the effects of refeeding after a fast would also provide us more information on the benefits of

fasting. Third, we only obtained samples of patients who fasted for 26 hours. We know that

glycogen levels aren’t usually depleted until around ~24 hours after a meal, and this can vary

greatly. It would be interesting to analyze samples of patients who have fasted for longer

amounts of time such as 72+ hours, giving us a greater understanding of the effects of

prolonged fasting. Fourth, the samples were obtained from subcutaneous adipose tissue. How

does each different tissue type respond to nutrient deprivation? It would be beneficial to analyze

how different tissues respond to fasting and the different pathways that are involved.





References

\1. Bailey EE, Pfeifer H, Thiele EA. The Use of Diet in the Treatment of Epilepsy. Epilepsy and Behavior.

2005; 6: 4-8.

\2. Kaeberlein, T.L., Smith, E.D., Tsuchiya, M., Welton, K.L., Thomas, J.H., Fields, S., Kennedy, B.K. and

Kaeberlein, M. (2006), Lifespan extension in Caenorhabditis elegans by complete removal of food.

Aging Cell, 5: 487-494. <https://doi.org/10.1111/j.1474-9726.2006.00238.x>

\3. Mitchell SJ, Bernier M, Mattison JA, et al. Daily Fasting Improves Health and Survival in Male Mice

Independent of Diet Composition and Calories. Cell Metab. 2019;29(1):221-228.e3.

doi:10.1016/j.cmet.2018.08.011

\4. Mattson MP, Longo VD, Harvie M. Impact of intermittent fasting on health and disease processes.

Ageing Res Rev. 2017;39:46-58. doi:10.1016/j.arr.2016.10.005

\5. Defour M, Michielsen CCJR, O'Donovan SD, Afman LA, Kersten S. Transcriptomic signature of fasting

in human adipose tissue. Physiol Genomics. 2020 Oct 1;52(10):451-467. doi:

10.1152/physiolgenomics.00083.2020. Epub 2020 Aug 31. PMID: 32866087.

\6. Chang, Kai-Ming & Harbron, Chris & South, Marie. (2006). An Exploration of Extensions to the RMA

Algorithm.

\7. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential

expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47.

doi: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007).

\8. MacDonald JW (2017). hugene21sttranscriptcluster.db: Affymetrix hugene21 annotation data (chip

hugene21sttranscriptcluster). R package version 8.7.0.

\9. Oberkofler H, Fukushima N, Esterbauer H, Krempler F, Patsch W. Sterol regulatory element binding

proteins: relationship of adipose tissue gene expression with obesity in humans. Biochim Biophys Acta.

2002 May 3;1575(1-3):75-81. doi: 10.1016/s0167-4781(02)00279-8. PMID: 12020821.

\10. Ferré P, Foufelle F: SREBP-1c Transcription Factor and Lipid Homeostasis: Clinical Perspective. Horm

Res 2007;68:72-82. doi: 10.1159/000100426

\11. Wang, Y., Zhang, C., Jin, Y. et al. Alkaline ceramidase 2 is a novel direct target of p53 and induces

autophagy and apoptosis through ROS generation. Sci Rep 7, 44573 (2017).

<https://doi.org/10.1038/srep44573>

\12. Liu B, Xiao J, Dong M, Qiu Z, Jin J. Human alkaline ceramidase 2 promotes the growth, invasion, and

migration of hepatocellular carcinoma cells via sphingomyelin phosphodiesterase acid-like 3B. Cancer

Sci. 2020;111(7):2259-2274. doi:10.1111/cas.14453

\13. Chen EY, Ta n CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A. Enrichr: interactive

and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14)





\14. Kuleshov M V, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL,

Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A. Enrichr: a

comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016;

gkw377 .

\15. Xie Z, Bailey A, Kuleshov M V, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz

ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A. Gene set knowledge discovery with Enrichr.

Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90

