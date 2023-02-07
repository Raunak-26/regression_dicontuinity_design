rdd_econ_stats
==============

# Description

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

# Introduction

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

<p align="center">
  <img src="https://github.com/vanshajbindlish/regression_dicontuinity_design/blob/main/paper/Intro%20Sharp%20RDD.png" alt="Sharp Regression Discontinuity">
  <br>
  <em>Sourced from Imbens and Lemieux (2007). X-Axis is the forcing variable, while Y-Axis is the outcome value</em>
</p> 

In the above equation the problem is the estimate the values of $\mu_+$ and $\mu_-$

$$\qquad \mu_{+} = \lim_{x\downarrow c} \mathbb{E}\left[Y_{i} | X_{i}=x\right] \\qquad \text{and}, \mu_{-} = \lim_{x\uparrow c} \mathbb{E}\left[Y_{i} | X_{i}=x\right]$$

Since researchers are often unsure about the efficiency of the estimates based on global linear function, 
they often use polynomial functions of order fifth or sixth of the assignment variable, $X$ in the regression 
model. The appropriate degree of polynomial is mainly selected based on some statistical information 
criteria (like goodness-of-fit) or cross-validation.
    
The other approach to estimate $\mu_{+}$ and $\mu_{-}$ and consequently the treatment effect is to use local 
low-order (linear or quadratic) polynomial approximation. In this method, the values of $x_{i}$ around the 
threshold $c$ that lie within the bandwidth $h$ are used to estimate a linear or quadratic function, while 
the rest of the values of $x_{i}$ are discarded. Several methods have been developed to estimate the value 
of $g$ (see Imbens and Kalyanaraman(2012), Hahn, Todd, and Van Der Klaauw (2001), Lee and Lemieux(2010)).

Gelman and Imbens(2019)\cite{gelman2019high} believe that the former approach, that is, to use the global 
high-order polynomials to estimate the causal treatment effect is flawed and rather the inference should 
be based on local low-order polynomials. The authors bring forth threefold arguments to show global 
high-order polynomials as a poor model choice for regression discontinuity analysis.

\begin{itemize}
\item The first argument is that the estimates of the treatment effect can be driven by the values of 
the forcing variable that are away from the threshold. Since the RD estimate is the difference between 
the weighted average of the outcome for the treated and untreated, wherein the weights depend only on 
the forcing variable; the weights in a global high-order polynomial approximation are extreme and 
unattractive compared to a local linear/quadratic approximation.
\item The second argument put forward is that the estimates in the case of global high-order polynomials 
are highly sensitive to the chosen degree of the polynomial of the forcing variable $X$. Further, a lack 
of just criteria to select the order of polynomials that is suitable for the estimator to satisfy its 
objectives of causal inference exacerbates the issue.
\item Finally, the third argument is that the confidence interval based on the estimates of global 
high-order polynomials is often misleading. The confidence intervals can be too narrow that they reject 
the null hypothesis even when it shouldn't have, that is, the probability of Type 1 error rate is much 
higher than the nominal rate of 5\%. This implies that there is a higher bias in the case of global 
high-order polynomials to detect a discontinuity even if there isn't one, and hence, generate poor inference.
\end{itemize}

The target of our paper is to test the veracity of these arguments through a simulation study and 
empirical application (wherever possible). The paper hereby proceeds as follows. In the next section, 
we introduce the model specification of the RD design that is used to critique the aforementioned arguments. 
Thereafter, we evaluate our results based on the simulation and empirical evidence so as to comment on 
the feasibility of global high-order polynomials in RD design.

This introduction is the brief overview of the regression discontibuity design and why high-order polynomials
should not be used in regression discontuinity designs. For forther exploration of the paper you can see
the rdd_paper.pdf in paper repository.
