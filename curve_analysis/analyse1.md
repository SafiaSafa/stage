---
title: "analyse1"
author: "Safia Safa-tahar-henni"
date: "January 30, 2018"
output: pdf_document
---

## Analyse sur les composantes principales.


### Classification (Cas vs. Controle)

#### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.

test avec niter=10,nb_pc=3


![Resultat courbe PC](curve_analysis/Pvalue/result_rank_pvalue.png)

* _Calcul pente / plateau_

#### __Prediction VS. #PC__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de PCs selon le nombre de SNPs étudié (500/100/200/500/100/200).


test avec niter=5,nb_pc=3


![Resultat courbe SNPs](script/fct_SNP.pdf)

### Regression (Age at onset)

#### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.
* _Calcul pente / plateau_

#### __Prediction VS. #PC__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de PCs selon le nombre de SNPs étudié (500/100/200/500/100/200).

