"""
association test on a vcf.
please run:  python vcfassoc.py --help
for usage
"""
from __future__ import print_function
import os.path as op
import sys
import toolshed as ts
import pandas as pd
import statsmodels
import statsmodels.api as sm
from statsmodels.genmod.dependence_structures import Exchangeable
from scipy.stats import chi2_contingency

import numpy as np
import patsy
import click
from itertools import chain
import re


def _get_gt(gt, splitter = re.compile("/|\|")):
    """
    return 0, 1, or 2. max of 2 because additional alts are just used as '1'

    >>> _get_gt("0/1")
    1
    >>> _get_gt("2/1")
    2
    >>> _get_gt("1/1")
    2
    >>> _get_gt("0/0")
    0
    >>> _get_gt(".")
    nan
    """
    if gt == ".": return np.nan
    return sum(min(int(n), 1) for n in splitter.split(gt))

def tfloat(n):
    try: return float(n)
    except ValueError: return np.nan

def _get_genotypes(vcf, min_qual, min_genotype_qual, min_samples, as_vcf):

    fh = ts.nopen(vcf)
    if as_vcf:
        for header in fh:
            print(header.rstrip("\r\n"))
            if header.startswith("#CHROM"):
                header = header.split("\t")
                header[0] = "CHROM"
                break
        else:
            1/0
        vcf_iter = ts.reader(chain(["\t".join(header)], fh), header="ordered")
    else:
        vcf_iter = ts.reader(vcf, skip_while=lambda l: l[0] !="#CHROM",
                            header="ordered")

    for i, variant in enumerate(vcf_iter):
        if float(variant['QUAL']) < min_qual: continue
        genotypes = variant.values()[9:]

        # it was a fraction
        if i == 0 and min_samples < 1:
            min_samples *= len(genotypes)

        # too few samples with data.
        if len(genotypes) - sum(gt.startswith('./.') or gt == "." or all ("." == v for v in
            gt.split(":")) for gt in genotypes) < min_samples:
            continue
        gs = [dict(zip(variant['FORMAT'].split(":"), gt.split(":"))) for gt in
                genotypes]

        gqs = [tfloat(d.get('GQ', 20)) for d in gs]
        gqs = [g if g > min_genotype_qual else np.nan for g in gqs]
        if (~np.isnan(gqs)).sum() < min_samples: continue
        gts = [_get_gt(d['GT']) for d in gs]
        gts = [g if q > min_genotype_qual else np.nan for g, q
                                                      in zip(gts, gqs)]
        yield variant.keys()[9:], gts, gqs, variant

def xtab(formula, covariate_df):
    y, X = patsy.dmatrices(str(formula), covariate_df)
    X = patsy.dmatrix('genotype', covariate_df)
    ix = get_genotype_ix(X)

    tbl = pd.crosstab(X[:, ix], y.ravel())
    tbl.columns = ['%s_%i' % (y.design_info.column_names[-1], j) for j in range(2)]
    tbl.index = ['%i_alts' % i for i in tbl.index]
    alts = set(tbl.index)
    if len(alts) < 2 or not '0_alts' in alts:
        tbl_dom = None
    else:
        tbl_dom = pd.DataFrame({'0_alts': tbl.ix['0_alts', :], 'n_alts': tbl.ix[list(alts - set(['0_alts'])), :].sum()}).T

    if not '2_alts' in alts or len(alts) < 2:
        tbl_rec = None
    else:
        tbl_rec = pd.DataFrame({'lt2_alts': tbl.ix[['0_alts', '1_alts'], :].sum(), '2_alts': tbl.ix['2_alts', :]})

    for name, xtbl in (('additive', tbl), ('dominant', tbl_dom), ('recessive', tbl_rec)):
        if xtbl is None:
            d['p.chi.%s' % name] = 'nan'
            continue

        chi, p, ddof, e = chi2_contingency(xtbl)
        if name == 'additive':
            d = xtbl.to_dict()
        d['p.chi.%s' % name] = "%.3g" % p
    return d

def get_genotype_ix(X):
    ix = [i for i, x in enumerate(X.design_info.column_names) if "genotype" in x]
    #print(X[:10])
    cols = X.design_info.column_names
    assert len(ix) == 1 or ('False' in cols[0] and 'True' in cols[1]), cols
    ix = ix[-1]
    return ix

def vcfassoc(formula, covariate_df, groups=None):

    y, X = patsy.dmatrices(str(formula), covariate_df)
    # get the column containing genotype
    ix = get_genotype_ix(X)


    model = sm.GLM(y, X, missing='drop', family=sm.families.Binomial())
    if groups:
        #covariate_df['grps'] = map(str, range(len(covariate_df) / 8)) * 8
        cov = Exchangeable()
        model = sm.GEE(y, X, groups=covariate_df[groups], cov_struct=cov,
                    family=sm.families.Binomial())
    result = model.fit(maxiter=1000)
    res = {'OR': np.exp(result.params[ix]),
           'pvalue': result.pvalues[ix],
           'z': result.tvalues[ix],
           'OR_CI': tuple(np.exp(result.conf_int()[ix, :])),
           }
    try:
        res['df_resid'] = result.df_resid
    except AttributeError:
        pass
    return res

def print_result(res, variant, as_vcf, i):
    """if as_vcf is True, append to the info field
    otherwise, print out tab-delimited info"""
    fmt = "{site}\t{pvalue}\t{REF/ALT}\t{OR}\t{OR_CI}\t{z}\t{p_chi_additive}\t{p_chi_dominant}\t{p_chi_recessive}"
    if 'df_resid' in res:
        fmt += "\t{df_resid}"
    fmt += "\t{xtab}\t{INFO}"
    if i == 0 and not as_vcf:
        print(fmt.replace("}", "").replace("{", ""))
    res['site'] = "{CHROM}:{POS}".format(**variant)
    if not variant["ID"] in (".", ""):
        res['site'] += "(%s)" % variant['ID']
    res['pvalue'] = "%.4g" % res['pvalue']
    res['INFO'] = variant['INFO']
    res['REF/ALT'] = variant['REF'] + "/" + variant['ALT']
    res['OR_CI'] = "%.4f..%.4f" % res['OR_CI']
    for m in 'additive dominant recessive'.split():
        res['p_chi_%s' % m] = res['xtab'].pop('p.chi.%s' % m)


    for k in 'z', 'OR':
        if k in res:
            res[k] = "%.3f" % res[k]
    if not as_vcf:
        print(fmt.format(**res))
    else:
        keys = fmt.replace("}", "").replace("{", "").split("\t")[1:]

        variant['INFO'] += ";" + ";".join("%s=%s" % (k, res[k]) for k in keys)
        print("\t".join(variant.values()))

# weighted is one of False, 'GQ', log, log10 to do weighted regression be
# genotype quality
@click.command(help="""\
        vcfassoc: association testing with logistic regression on a vcf

first column in [covariates] must be sample ids that match the sample columns
in the [vcf].

[formula] is an R (or patsy)-like formula.

Example:

    python vcfassoc.py multisample.vcf study-covariates.txt \
        "status ~ genotype + age" --min-samples 10

A 'genotype' column is added automatically by this program. 'status'
must be a binary variable defined in [covariates] age must also be in
[covariates] the returned p-value and odds-ratio will always be fore
the genotype variable.

To do a dominant model, use:

        "status ~ I(genotype > 0) + age"

The default is to output a tab-delimited table with 1 row per variable
with only the regression information. But, the regression values can be
added to the [vcf] with --as-vcf.
       """)
@click.argument('vcf', type=str)
@click.argument('covariates', type=click.Path(exists=True))
@click.argument('formula', type=str)
@click.option('--min-qual', default=1, help="skip variants with QUAL < this",
        type=click.IntRange(0, 255))
@click.option('--min-genotype-qual', default=1, help="set samples with"
        "genotype qual less than this to 'NA'", type=click.IntRange(0, 255))
@click.option('--min-samples', default=2, help='skip variants where fewer than'
        'this many samples have usable genotypes')
#@click.option('weighted', type=click.Choice(('log10', 'log', 'GQ', 'FALSE')),
#        default='FALSE')
@click.option('--as-vcf', is_flag=True, help="output the entire VCF with model stuff"
        "added to the INFO field instead of default of writing a table of results")
@click.option('--groups', help="a grouping column in `covariates` used in "
        "fitting a GEE with exchangeable correlation. Useful for families or "
        "paired samples. "
        "(see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3736986/)")
def main(vcf, covariates, formula, min_qual, min_genotype_qual, min_samples,
        weighted=False, as_vcf=False, groups=None):
    #if weighted == "FALSE": weighted = False
    #else:
    #    weight_fn = {'log10': np.log10, 'log': np.log, 'GQ': np.array}[weighted]
    if covariates.endswith('.csv'):
        covariate_df = pd.read_csv(covariates, index_col=0)
    else:
        covariate_df = pd.read_table(covariates, index_col=0)
    covariate_df.index = [str(x) for x in covariate_df.index]
    gmatrix = {}
    for i, (samples, genos, quals, variant) in enumerate(
            _get_genotypes(vcf, min_qual, min_genotype_qual, min_samples,
                           as_vcf)):
        if i == 0:
            # make sure we have covariates for all samples in the vcf
            assert not set(samples).difference(covariate_df.index),\
                        set(samples).difference(covariate_df.index)
            covariate_df = covariate_df.ix[samples,:]
        covariate_df['genotype'] = genos
        #if weighted:
        #    covariate_df['weights'] = weight_fn(quals)
        try:
            res = vcfassoc(formula, covariate_df, groups)
            gmatrix['{CHROM}:{POS}'.format(**variant)] = genos
        except np.linalg.linalg.LinAlgError:
            continue
        except statsmodels.tools.sm_exceptions.PerfectSeparationError:
            print("WARNING: perfect separation, too few samples(?)",
                  ": setting to -9: {CHROM}:{POS}".format(**variant),
                  file=sys.stderr)
            res['z'] = res['OR'] = np.nan
            res['pvalue'] = -9.0 # blech.
            res['OR_CI'] = np.nan, np.nan
            gmatrix['{CHROM}:{POS}'.format(**variant)] = genos
        except IndexError:
            continue
        res['xtab'] = xtab(formula, covariate_df)
        print_result(res, variant, as_vcf, i)

    l1_regr(pd.DataFrame(gmatrix), covariate_df, formula)

def l1_regr(genotypes, covariate_df, formula, C=0.1):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import Imputer, StandardScaler
    assert genotypes.shape[0] == covariate_df.shape[0]

    genotypes.index = covariate_df.index
    genotypes.ix[:, :] = Imputer(missing_values="NaN",
                                 strategy="most_frequent",
                                 axis=0).fit_transform(genotypes)
    y, covariates = patsy.dmatrices(str(formula), covariate_df,
            return_type='dataframe')

    covariates.drop([x for x in covariates.columns if 'genotype' in x], inplace=True, axis=1)
    covariates.ix[:, :] = StandardScaler().fit_transform(covariates)
    X = pd.concat((covariates, genotypes), axis=1)
    #X.to_csv('/tmp/X.csv')
    #y.to_csv('/tmp/y.csv')

    assert X.shape[0] == covariate_df.shape[0]
    # or drop variants with missing data.
    #X = X.dropna(axis=1)
    tol, do_break = 1e-4, True
    for p in ('l1', 'l2'):
        for c in (1e-3, 0.01, 0.1, 0.25, 0.5, 1, 2, 10, 1000, 1e8):
            clf = LogisticRegression(C=c, penalty=p, tol=tol,
                        fit_intercept=False)
            try:
                clf.fit_transform(X, y.squeeze())
            except ValueError:
                continue
            coefs = clf.coef_.ravel()
            variants = list(X.columns[coefs != 0])

            scores = list(coefs[coefs != 0])
            #n = X_new.shape[1]
            n = len(variants)
            print("penalty: {p} C: {c} n_covariates: {n}".format(**locals()), file=sys.stderr)
            print([(v, "%.3f" % s) for s, v in sorted(zip(np.abs(scores), variants),
                reverse=True)], file=sys.stderr)
            break
        else:
            do_break = False
        if do_break: break

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    main()
