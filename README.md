rdd_econ_stats
==============

## Description

Controlling for higher degrees of polynomials for the forcing variable is a common 
practice in regression discontinuity analysis. Such models can overfit, resulting 
in substantively implausible causal inferences that randomly contribute variation 
to the high-degree polynomial and can lead to false estimates. The aim of this 
paper is to investigate the three arguments presented in Gelman and Imbens (2019). 

The arguments recommends against the use of high-order polynomial approximations 
are: 

* the implicit weights are inappropriate
* the estimates are sensitive to the degree of polynomial, and
* poor coverage of confidence intervals. 

Instead of using global high-order polynomial, the paper suggests using local 
linear or quadratic polynomial approximations.

## Introduction

Regression discontinuity designs have significantly gained importance in econometric 
and causal inference over the past few decades. Introduced in 1960 by Donald L.
Thistlethwaite and Donald T. Campbell, the regression discontinuity design aims to 
estimate the impact of a treatment policy on the variable of interest. The simplicity
of the design enables its application across research fields, like health, labor, 
education, political economy, etc.
    
The basic idea behind the regression discontinuity (RD) design is that the allotment 
of the treatment is based on the value of the predictor (also known as assignment/
forcing variable) being on either side of a fixed threshold. Further, the assignment
variable may be related to the potential outcome, however, this relation is assumed 
to be smooth which implies that any discontinuity at the cutoff point in the 
conditional expectation of the outcome as a function of the assignment variable is 
because of the treatment and this discontinuity gap represents the causal effect of 
the treatment. For instance, Thistlethwaite and Campbell studied the role of merit 
awards on future academic outcomes by allocating awards only to the students who had 
scored above a threshold score value in a test, while the students who scored below 
the cutoff point did not receive any award.
    
In regression discontintuity, the average treatment effect can be estimated as the 
discontinuity in the conditional expectation of the outcome given the assignment 
variable at the threshold. The equaiton then can be represented as:

$$ \tau_{SRD} = \lim_{x\downarrow c} \mathbb{E}[Y_{i} | X_{i}=x] - \lim_{x\uparrow c} \mathbb{E}\left[Y_{i} | X_{i}=x\right]$$

which can be simply written as:

$$ \tau_{SRD} = \mathbb{E}\left[Y_{i}(1) - Y_{i}(0)| X = c\right]$$

Here, $\tau_{SRD}$ is the 'Average Treatment Effect' at threshold $c$.

This can also be seen in figure below, which illustrates a sharp RD design with the threshold 
at 6. Individuals right to the threshold are given the treatment, while those to the left 
are the control group. The dashed line represents the value of the potential outcome for the 
treated group and control group if they were observable beyond their treatment status. The 
smooth lines indicate the observed outcome conditional on the assignment variable. Under the 
continuity assumption, therefore, the average treatment effect, τSRD is given by the jump in 
outcome value at threshold 6.

    

## Credits

This project was created with [cookiecutter](https://github.com/audreyr/cookiecutter)
and the
econ-project-templates](https://github.com/OpenSourceEconomics/econ-project-templates).
