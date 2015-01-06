vcfassoc
========

perform phenotype::genotype-association tests on a VCF with logistic regression.

Quickstart
==========

(also see Installation)

You will need:
 1. A multi-sample VCF
 2. A covariates file with columns you wish to model. (currently the first column
    of this file must match the sample headings in the VCF.)
 3. A formula in R (patsy) syntax that will pull from 2. Must contain 'genotype'
    which will be added to the covariates on-the-fly.

Example run with **additive** model:

```Bash
vcfassoc.py samples.vcf covariates.txt 'status ~ genotype + age + gender + PC1 + PC2'
```

Example run with **dominant** model:

```Bash
vcfassoc.py samples.vcf covariates.txt 'status ~ I(genotype > 0) + age + gender + PC1 + PC2'
```

Example run with **recessive** model:

```Bash
vcfassoc.py samples.vcf covariates.txt 'status ~ I(genotype == 2) + age + gender + PC1 + PC2'
```

The output will look like:

site                | pvalue   | REF/ALT  | OR              | OR\_CI            | z       | p\_chi\_additive  | p\_chi\_dominant  | df\_resid
------------------  | -------- | -------  | --------------- | ---------------- | --------| --------------  | --------------- | --------
1:256769(rs166498)  | 0.1418   | C/T      | 5.000           | 0.5842..42.7971  | 1.469   | 0.174           | 0.174           | 31
1:257063(rs194882)  | 0.9994   | G/A      | 0.000           | 0.0000..inf      |-0.001  | 1               | 1               | 31
1:257334(rs3959)    | 0.8998   | A/G      | 1.032           | 0.6298..1.6919   | 0.126   | 0.875           | 1               | 31
1:257454(rs70148)   | 0.5714   | C/G      | 0.500           | 0.0453..5.5141   |-0.566  | 1               | 1               | 31
1:257735(rs134134)  | 0.9994   | T/G      | 2322868131.974  | 0.0000..inf      | 0.001   | 1               | 1               | 31
1:257925(rs39944)   | 0.05584  | C/T      | 4.076           | 0.9656..17.2084  | 1.912   | 0.0465          | 0.0431          | 27
1:258014(rs3283)    | 0.4235   | A/G      | 2.000           | 0.3663..10.9192  | 0.800   | 0.651           | 0.651           | 31
1:258064(rs3091)    | 0.1863   | G/C      | 0.427           | 0.1207..1.5090   |-1.321  | 0.206           | 0.202           | 30
1:258091(rs3092)    | 0.438    | A/G      | 0.730           | 0.3294..1.6175   |-0.776  | 0.473           | 0.6             | 30

Where the `p_chi_*` columns are the result of running a chisq test under the additive or dominant (and now recessive) 
models. `pvalue` is for the specified model under logistic regression. the `OR_CI` is the 95% confidence interval
for the odds-ratio of genotype from the model.

The INFO column from the VCF and contingency table for the chisq test are also printed, but not shown here.

If *--as-vcf*, is specified, the above info will be crammed into the INFO field and the output will be
in VCF format.

Installation
============

`vcfassoc` requires a number of python modules: numpy, scipy, statsmodels, toolshed,
	patsy, click, scikit-learn

All dependencies can be installed with pip via:

    pip install -r requirements.txt

For users who do not have the scientific python stack set up, it is recommended
to download and install [anaconda](https://store.continuum.io/cshop/anaconda/),
a free distribution that contains most of these modules.

You will also need to install [futures](https://pypi.python.org/pypi/futures).

Features
========

In addition to the obvious features, vcfassoc will also perform a feature-selection
using **'l1' penalized regression** and report the top features from that selection
to STDERR at the end of the run.

Further, it will **allow correlated data** via the use of generalized estimating
equations (**GEE**)'s. Just specify, e.g. `--group` 'family\_id' to indicate the grouping
and GEE will automatically be used instead of traditional logistic regression.

See: https://github.com/brentp/vcfassoc/blob/master/test/test.sh

`vcfassoc` has been tested against R.

TODO
====

 + option to output genotype matrix.
 + allow continuous dependent variable
 + parallelization (naive version didn't improve speed much)
 + weighted regression by genotype qualities
