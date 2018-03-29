# nscores
Executes the Van der Waerden version of non parametric tests (Normal scores tests) <br/>
Named for the Dutch mathematician Bartel Leendert van der Waerden, the Van der
Waerden test is a statistical test that k population distribution functions
are equal. The Van Der Waerden test converts the ranks to quantiles of the
standard normal distribution. These are called normal scores and the test is
computed from these normal scores. The standard ANOVA assumes that the errors
(i.e., residuals) are normally distributed. If this normality assumption is
not valid, an alternative is to use a non-parametric test. The advantage of
the Van Der Waerden test is that it provides the high efficiency of the
standard ANOVA analysis when the normality assumptions are in fact satisfied,
but it also provides the robustness of the non-parametric test when the
normality assumptions are not satisfied.
This function compute the Normal Scores of 5 tests:
- Levene, Mann-Whitney-Wilcoxon and Wilcoxon tests when there are 2 groups;
- Kruskal-Wallis and Friedman test whene there are more than 2 groups.

The function will use a GUI to select the proper test.
Moreover, the GUI will ask which version of Normal score do you want to
use: Blom, Tukey, Rankit, Van der Waerden

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2010). NSCORES: Normal scores version of several non-parametric tests.
http://www.mathworks.com/matlabcentral/fileexchange/26855
