![alt_text][logo]

# Perception in Functional Motor Disorders: *public repository*

###### Julian Matthews*, Kanae Nagao*, Catherine Ding, Rachel Newby, Peter Kempster, Jakob Hohwy

***

Functional dystonia is the most diagnostically challenging of the functional motor disorders (FMDs), stubbornly resistent to criteria that emphasise psychiatric symptoms. However, novel frameworks from cognitive neuroscience show promise. Drawing on writing by [Edwards et al. (2012)](https://www.ncbi.nlm.nih.gov/pubmed/22641838); [Stenner & Haggard (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27719833); and [Newby, Alty, & Kempster (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27753149) that speaks to FMDs by drawing on the [*predictive processing* framework of brain function](https://global.oup.com/academic/product/the-predictive-mind-9780199682737?cc=au&lang=en&), this project employs a [contemporary extension](http://rstb.royalsocietypublishing.org/content/373/1755/20170352) of a [psychophysical dual-task paradigm](https://www.ncbi.nlm.nih.gov/pubmed/25973773) to examine the relationship between perceptual and active inference in motor disordered patient groups. 

In particular, we were interested in four domains implicated in the *predictive processing* account:
* **Attention** —Precision optimisation
* **Expectations** —Predictive beliefs
* **Sensation** —Here, perceptual contrast sensitivity
* **Metacognition** —Introspective access to behavioural accuracy

Our study found evidence that perceptual sensitivity (i.e., *sensation*) is a domain of interest in motor disorders despite their ostensibly movement-oriented symptoms. Beyond this, our study is important because it applies cutting-edge psychophysics to examine these domains in a single paradigm and contrasts matched-controls with functional and organic patient groups within the same experiment.

## Task design
We employed an extended version of the dual-task:

![alt_text][methods]

> Method for letter and gabor task including 8AFC response to Gabor presence (named *grill* in the experiment to aid description for subjects). To manipulate attention, participants performed the Gabor Task alone (Full Attention condition) or in conjunction with a Letter Task (Diverted Attention condition). To manipulate expectations the presence of the peripheral gabor was altered between blocks (25%, 50%, and 75% likelihood). Subjects were instructed before each block of trials and after each trial about the probability of gabor presence. **δt** is stimulus-onset asynchrony (SOA), equivalent to the time between presentation of the letter stimulus and mask. In order to standardise the difficulty of the Letter Task, SOA timing was adjusted psychometrically for each participant during training to maintain 79.4% performance when conducted alone (see QUEST staircasing).

## What is this?
In the interests of open science, we provide the code and de-identified data from our perceptual study of functional and organic motor disorders. The supplied code includes critical scripts (written for MATLAB, R, and JASP) that can be used to run and analyse our experiment. I have endeavoured to provide extensive comments in this code so it is accessible to others. De-identified data includes the entire collection of behavioural responses (43,200 trials in total) from our participants along with preprocessed summary statistics used for behavioural and perceptual analysis. I have tried to use common language to describe this data but [please contact me](mailto:julian.r.matthews@gmail.com?subject=FMD%20study%20enquiry) if you have questions.

## You will need:
* You will need:
 
 [MATLAB](https://www.mathworks.com/products/matlab.html)
 
 [Psychtoolbox](http://psychtoolbox.org/)
* We recommend:
 
 [R](https://www.r-project.org/)
 
 [JASP](https://jasp-stats.org/)

## Getting started:
[`runExp.m`](./fmd-perceptual-study/scripts/experiment/) is the critical file for running the experiment.

***

![alt_text][avatar]

###### J. Matthews and K. Nagao contributed equally to this work

[logo]: https://cogphillab.files.wordpress.com/2018/08/header1.jpg "Cognition and Philosophy Lab"
[methods]: https://github.com/julian-matthews/fmd-public-repository/blob/master/fmd-perceptual-study/figures/figure1.png
[avatar]: https://avatars0.githubusercontent.com/u/18410581?v=3&s=96 "Julian Matthews"
