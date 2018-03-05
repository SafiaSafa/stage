---
title: "analyse1"
author: "Safia Safa-tahar-henni"
date: "January 30, 2018"
output: pdf_document
---

## Analyse sur les composantes principales.


### Classification (Cas vs. Controle)
#### __Ranking: Pvalue__
##### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.


![Resultat courbe PC](Pvalue/result/curve_all_snpsRplot.png)

* _Calcul pente / plateau_

![Resultat courbe PC](Pvalue/result/curve_zoom_all_snpRplot.png)

##### __Prediction VS. #PC__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de PCs selon le nombre de SNPs étudié (500/100/200/500/100/200).



![Resultat courbe SNPs](Pvalue/result/curve_pcRplot.png)

#### __Ranking: I-score__
##### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.


![Resultat courbe PC](iscore/prunning/index.png)

* _Calcul pente / plateau_

![Resultat courbe PC](iscore/prunning/zoomRplot.png)

##### __Prediction VS. #PC__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de PCs selon le nombre de SNPs étudié (500/100/200/500/100/200).



![Resultat courbe SNPs](iscore/prunning/fct_pc.png)

#### __Ranking: Pvalue + Iscore__
##### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.




![Resultat courbe PC](combined/index.png)

* _Calcul pente / plateau_

![Resultat courbe PC](combined/zoomRplot.png)

##### __Prediction VS. #PC__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de PCs selon le nombre de SNPs étudié (500/100/200/500/100/200).



![Resultat courbe SNPs](combined/fct_pc.png)
#### __Ranking: Beta__
##### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.




![Resultat courbe PC](beta/index.png)

### Regression (Age at onset)

#### __Prediction VS. #SNPs__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de SNPs selon le nombre de Composantes Principale (PC) ajouté en covariables.
* _Calcul pente / plateau_

#### __Prediction VS. #PC__
Mesure de la prédiction (ou coefficents de prédictivité) en fonction du nombre de PCs selon le nombre de SNPs étudié (500/100/200/500/100/200).

