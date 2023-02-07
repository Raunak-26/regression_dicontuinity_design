rdd_econ_stats
==============

# Abstract

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

* The first argument is that the estimates of the treatment effect can be driven by the values of 
the forcing variable that are away from the threshold. Since the RD estimate is the difference between 
the weighted average of the outcome for the treated and untreated, wherein the weights depend only on 
the forcing variable; the weights in a global high-order polynomial approximation are extreme and 
unattractive compared to a local linear/quadratic approximation.

* The second argument put forward is that the estimates in the case of global high-order polynomials 
are highly sensitive to the chosen degree of the polynomial of the forcing variable $X$. Further, a lack 
of just criteria to select the order of polynomials that is suitable for the estimator to satisfy its 
objectives of causal inference exacerbates the issue.

* Finally, the third argument is that the confidence interval based on the estimates of global 
high-order polynomials is often misleading. The confidence intervals can be too narrow that they reject 
the null hypothesis even when it shouldn't have, that is, the probability of Type 1 error rate is much 
higher than the nominal rate of 5\%. This implies that there is a higher bias in the case of global 
high-order polynomials to detect a discontinuity even if there isn't one, and hence, generate poor inference.

The target of our paper is to test the veracity of these arguments through a simulation study and 
empirical application (wherever possible). The paper hereby proceeds as follows. In the next section, 
we introduce the model specification of the RD design that is used to critique the aforementioned arguments. 
Thereafter, we evaluate our results based on the simulation and empirical evidence so as to comment on 
the feasibility of global high-order polynomials in RD design.

# Methodology

In practice, it is convenient to run a pooled regression to obtain the average treatment effect of the SRD 
rather than to estimate the values of $\mu_{+}$ and $\mu_{-}$ and then take their difference. This general 
 approach of pooled regression based on Lee and Lemieux (2010) can be represented as:This general approach 
 of pooled regression based on Lee and Lemieux (2010) can be represented as:

$$\qquad Y = \beta_{1} + \tau_{SRD} D + f(X-c) + \epsilon$$

Where, $\beta_{1}$ is the intercept term and equal to $\beta_{1}^{+} - \beta_{1}^{-}$\
$\tau_{SRD}$ represents the estimate of the treatment effect.\
$D$ is the dummy variable that takes value $1$ if $\{X \geq c \}$ and $0$ if $\{X < c \}$.\
$f(X-c)$ represents the functional form of the assignment variable of the regression function and is equivalent to \
$$\qquad f(X-c) = f^{-}(X-c) + D \left[f^{+}(X-c) - f^{-}(X-c) \right].$$

Hence, we can write the regression function in the form of: \
$$\qquad Y = \beta_{1} + \tau_{SRD}D + \beta_{2}(X-c) + \beta_{3}(X-c)^{2} + \beta_{4}D(X-c) + \beta_{5}D(X-c)^2 + \epsilon$$

The equation stated above can therefore be also considered the regression function for the global 
high-order polynomial approach. The local linear/quadratic approach is distinct from the former 
approach in the sense that for a certain bandwidth value $h$ and the consequent range around the 
threshold point $(c-h \leq X \leq c+h)$, the set of values of the forcing variable that lie outside 
of this range are dropped for the estimation of ATE. This implies that the regression function changes to:

$$\qquad Y = \beta_{1} + \tau_{SRD}D + \beta_{2}(X-c) + \beta_{3}(X-c)^{2} + \beta_{4}D(X-c) + \beta_{5}D(X-c)^2 + \epsilon$$\
where, $(c-h)\leq X \leq (c+h)$.

The advantage of such a regression form is that it introduces varying slopes on both sides of the threshold
through the interaction term between $X$ and $D$, which is often the case in empirical applications.

## Data Generating Process for Simulaiton

The objective of the simulation study is to test the arguments presented by Gelman and Imbens (2019)
in favor of the local low-order polynomial approximation for an RD analysis. We conduct the study for 
a Sharp Regression Discontinuity framework.

To conduct a Monte Carlo study, we first develop an SRD model. We consider a model with a constant 
average treatment effect, which means that the only change for the treated group compared to the control 
is of an additive shift  by a value equal to the treatment effect. Further, we consider a functional form 
of up to sixth-order polynomial. The model developed hereby is based on the data generating process in 
Imbens and Kalyanaraman 2012

The data generating process for $Y$, then can be defined as:

$$Y = 0.42 + 0.56D + 0.84X - 3.00X^{2} + 7.99X^{3} - 9.01X^{4} + 3.56X^{5} + 6.11X^{6} + e$$

where, the forcing variable is derived from a beta distribution $Z = B(\alpha = 2, \beta = 5)$. 
The forcing variable is then of the form:

$$X \sim 2^{*}Z - 1 $$

We introduce the threshold point at $c = 0$, which implies that all the individuals above this point 
are given treatment, while those below the cutoff point are not. We consider the true value of the 
estimate of treatment effect as \textbf{0.56}. And the treatment status is introduced in the data 
generating process through the dummy variable, $D$ where:

$$D = \begin{cases}1, & x_{i}\geq 0\\
0, & x_{i} < 0\ \end{cases}$$

And, the random error term follows a normal distribution with mean $0$ and standard deviation $0.1295$, 
that is,
 
$$e \sim \mathbb{N}(\mu = 0, \sigma = 0.1295)$$

Based on the above specification, we simulate a dataset of sample size, $n = 100,000$, which is then 
used to conduct an analysis of our objectives.

# Argument 1: Noisy Weigths

In an RD design, researchers are interested in estimating the impact of treatment at the cutoff point. 
Since the treatment effect estimates can be represented as the difference between the weighted average 
of outcomes for the treated and the control, one should expect an allocation of higher weights to the 
observations closer to the cutoff point. However, Gelman and Imbens (2019) argue that for the global 
polynomial approach, high and erratic weights are allotted to observations that are far from the 
threshold. This implies that the estimates calculated from such a model framework might be biased.

The weights are just a function of the assignment variable $X$. Thus, we can calculate the weights allotted 
to each observation of the forcing variable for different bandwidth values $h$ and the order of 
polynomial $K$ and examine them without even considering the functional form of the data or the outcome 
value. The weights can then be derived in a similar manner as Gelman and Imbens (2019) which is as follows:

$$\qquad w_{i} = 1_{c \leq x_{i} \leq h}.e^\prime_{K+1,1}
\begin{pmatrix}
\sum_{j:c \leq x_j \leq h}
\begin{pmatrix}
1 & x_j & \cdots & x_j^K \\
x_j & x_j^2 & \cdots & x_j^{K+1} \\
\vdots & \vdots & \ddots & \vdots \\
x_j^K & x_j^{K+1} & \cdots & x_j^{2K}
\end{pmatrix}
\end{pmatrix}^{-1} 
\begin{pmatrix}
1 \\  x_i \\ \vdots \\ x_i^K
\end{pmatrix}$$

Where, 
* The values of the forcing variable, $X$, are above the threshold $c$, and within the bandwidth $h$
* $e$ is the K+1 column vector, with the first element equal to 1 and the rest elements equal to $0$
* 1$ denotes the sum vector with its length equal to the number of values of forcing variable between 
  $c$ and $h$.
* $w_{i}$ is the weight allotted to each individual $i$ of the assignment variable $X$ for bandwidth 
  value $h$ and order of polynomial $K$. We can then obtain the weight vector, $w$, with its elements 
  equal to the weight allotted to each individual $i$ of the assignment variable $X$ for bandwidth 
  value $h$ and order of polynomial $K$. The average weight of the vector has been normalized to 1.
  
Given the data generating process and the weights functions defined in the above formula, we calculated
the weigths results for global high-order polynomial and local polynomial. In case of the local linear 
and quadratic regressoin the bandwidth is calculated based on Imbens and Kalyanaraman (2012) optimal 
bandwidth selector. These weight vectors are thereby plotted with the observations of forcing variables 
above the threshold on the X-Axis and the corresponding weights on the Y-Axis which can be shown in the 
figure below.

<p align="center">
  <img src="https://github.com/vanshajbindlish/regression_dicontuinity_design/blob/main/paper/arg1_sim.png" alt="Argument 1 Simulation Global Weights"> 
  <br>
  <em>Weights for Global High-Order Polynomial - Simulation Study</em>
</p> 

As can be seen from figure, the weights associated with global high-order polynomials are indeed high 
and erratic irrespective of the polynomial order. Further, it is also noticeable that higher weights 
are allotted to the observations away from the threshold. We also observe that for large values of the 
forcing variable, the associated weights are quite sensitive to the order of polynomials. Finally, similar
to the global polynomial approach, we plot the derived weight vectors in order to inspect the weights 
assigned to each observation of the forcing variable that lies above the threshold.

<p align="center">
  <img src="https://github.com/vanshajbindlish/regression_dicontuinity_design/blob/main/paper/arg1_sim_local.png" alt="Argument 1 Simulation Local Weights"> 
  <br>
  <em>Weights for Local Linear and Quadratic Polynomial - Simulation Study</em>
</p> 

Further, we perform the above stated steps for empirical application as well. We use the New York: 
MDRC, 2012 dataset for this study. we arrive at similar results as of the simulation study. The detailed
analysis of this empirical application can be found in the <a href= "https://github.com/vanshajbindlish/regression_dicontuinity_design/tree/main/paper">paper</a> itself.

# Team

Rauank Mehrotra -- s6rkmehr@uni-bonn.de \
Vanshaj Bindlish -- s6vabind@uni-bonn.de
