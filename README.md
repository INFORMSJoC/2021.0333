[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Asymptotically Optimal Sampling Policy for Selecting Top-m Alternatives

This repository contains supporting material for the paper [Asymptotically Optimal Sampling Policy for Selecting Top-m Alternatives](https://doi.org/????) by Gongbo Zhang, Yijie Peng<sup>*</sup>, Jianghua Zhang, and Enlu Zhou.

## Disclaimer

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

WARNING: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability, clarity, generality, and efficiency has not been well considered.

## Cite

To cite this code, please cite the paper using its DOI and the code itself, using the following DOI: DOI:10.1287/ijoc.2021.0333.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{ZhangPengZhangZhou2021.0333,
  author = {Zhang, Gongbo and Peng, Yijie and Zhang, Jianghua and Zhou, Enlu},
  publisher = {INFORMS Journal on Computing},
  title = {Asymptotically Optimal Sampling Policy for Selecting Top-m Alternatives},
  year = {2023},
  doi = {10.1287/ijoc.2021.0333.cd},
  url = {https://github.com/INFORMSJoC/2021.0333},
  note={available for download at https://github.com/INFORMSJoC/2021.0333}
}
```

## Content

This repository includes

1. The [main](main) folder contains the MATLAB implementations of the numerical experiments in the mian text.
    * The [Exp](main/Exp) folder contains codes for compared sampling policies;
    * The [Exp5](main/Exp5) folder contains codes for numerical experiments in Section 5.2 of the main text.
2. The [App](App) folder contains the MATLAB implementations of the numerical experiments in the online Appendix.
    * The [OCBA](App/OCBA) folder contains codes for numerical experiments in Section A.1 of the online Appendix;
    * The [Sampling ratio](App/Sampling%20ratio) folder contains codes for numerical experiments in Section A.3 of the online Appendix;
    * The [Non-normal](App/Non-normal) folder contains codes for numerical experiments in Section A.4.2 of the online Appendix;
    * The [Top-m arms identification](App/Top-m%20arms%20identification) folder contains codes for numerical experiments in Section A.4.3 of the online Appendix.
3. The [Results](Results) folder contains raw results of our numerical experiments.


## Installation

* The codes were written and run in MATLAB R2021a, on Windows 10 Home 64-bit OS, with Intel i5-11300H CPU @ 3.10 GHz, 16 GB RAM.

* To install the MATLAB codes, just copy the entire folder [main](main) and [App](App), respectively, into your MATLAB directory, or change the path of MATLAB to the folder [main](main) and [App](App), respectively.

## Details

#### Numerical Experiments in Section 5.1 of the main text

* The function [AOAPm.m](/main/Exp/AOAPm.m) is the proposed policy of our work;
* The function [EAm.m](/main/Exp/EAm.m) is a compared policy called equal allocation (EA);
* The function [OCBAmsa.m](/main/Exp/OCBAmsa.m) is a modified compared policy from Chen et al. (2008), where a better separating parameter that improves performance of the original policy is sequentially implemented;
* The function [OCBAmjia.m](/main/Exp/OCBAmjia.m) is a compared policy from Zhang et al. (2012,2015);
* The function [OCBAss.m](/main/Exp/OCBAss.m) is a compared policy from Gao and Chen (2015);
* The function [OCBASSS.m](/main/Exp/OCBASSS.m) is a compared policy from Gao and Chen (2016).

Remark:
* The policies can adjust the ascending or descending performance of alternatives;
* The policies can incoporate prior information or not;
* The policies can adjust sampling distributions;
* Update parameters in the same way when comparing;
* The independent macro experiments can run in parallel.

Get into the folder [Exp](/main/Exp). Run [EAm.m](/main/Exp/EAm.m) for *EA*, [OCBAmsa.m](/main/Exp/OCBAmsa.m) for *OCBAm*, [OCBAmjia.m](/main/Exp/OCBAmjia.m) for *OCBAm+*, [OCBAss.m](/main/Exp/OCBAss.m) for *OCBAss*, [OCBASSS.m](/main/Exp/OCBASSS.m) for *OCBASS*, and [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*.

* Set corresponding input parameters for all policies in the folder [Exp](/main/Exp).

#### Numerical Experiments in Section 5.2 of the main text

* The function [evacuation.m](/main/Exp5/evacuation.m) is a simulator of the evacuation network;
* The function [Untitled1.m](/main/Exp5/Untitled1.m) estimates the mean performance of each evacuation plan.

Get into the folder [Exp5](/main/Exp5). Run [EAm.m](/main/Exp/EAm.m) for *EA*, [OCBAmsa.m](/main/Exp/OCBAmsa.m) for *OCBAm*, [OCBAmjia.m](/main/Exp/OCBAmjia.m) for *OCBAm+*, [OCBAss.m](/main/Exp/OCBAss.m) for *OCBAss*, [OCBASSS.m](/main/Exp/OCBASSS.m) for *OCBASS*, and [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*.

* Set corresponding input parameters for all policies in the folder [Exp5](/main/Exp5);
* The input parameter *truemu* of each policy is calculated by running [Untitled1.m](/main/Exp5/Untitled1.m), which calls [evacuation.m](/main/Exp5/evacuation.m) during execution.

#### Numerical Experiments in Section A.1 of the online Appendix

* The function [OCBAm.m](/App/OCBA/OCBAm.m) is a compared policy from Chen et al. (2008);

* The function [OCBAms.m](/App/OCBA/OCBAms.m) is a modified compared policy from Chen et al. (2008), where the original policy is implemeted in a sequential manner.

Get into the folder [OCBA](/App/OCBA). Run [EAm.m](/main/Exp/EAm.m) for *EA*, [OCBAm.m](/App/OCBA/OCBAm.m) for *OCBAm(two-stage)*, [OCBAms.m](/App/OCBA/OCBAms.m) for *OCBAm(sequential)*, [OCBAmsa.m](/main/Exp/OCBAmsa.m) for *OCBAm*, and [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*.

* Set corresponding input parameters for all policies in the folder [OCBA](/App/OCBA) and [Exp](/main/Exp).

#### Numerical Experiments in Section A.3 of the online Appendix

* The function [AOAPm.m](/App/Sampling%20ratio/AOAPm.m) is our proposed policy, where the output is the sampling ratio for each simulation budget.

Get into folder [Sampling ratio](/App/Sampling%20ratio). Run [AOAPm.m](/App/Sampling%20ratio/AOAPm.m) for calculating its sampling ratios.

#### Numerical Experiments in Section A.4.1 of the online Appendix

Run [EAm.m](/main/Exp/EAm.m) for *EA*, [OCBAmsa.m](/main/Exp/OCBAmsa.m) for *OCBAm*, [OCBAmjia.m](/main/Exp/OCBAmjia.m) for *OCBAm+*, [OCBAss.m](/main/Exp/OCBAss.m) for *OCBAsst*, [OCBASSS.m](/main/Exp/OCBASSS.m) for *OCBASSt*, and [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*.

* Set corresponding input parameters for all policies.

#### Numerical Experiments in Section A.4.2.1 and A.4.2.2 of the online Appendix

Run [EAm.m](/main/Exp/EAm.m) for *EA*, [OCBAmsa.m](/main/Exp/OCBAmsa.m) for *OCBAm*, [OCBAmjia.m](/main/Exp/OCBAmjia.m) for *OCBAm+*, [OCBAss.m](/main/Exp/OCBAss.m) for *OCBAss*, [OCBASSS.m](/main/Exp/OCBASSS.m) for *OCBASS*, and [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*.

* Set corresponding input parameters for all policies.

#### Numerical Experiments in Section A.4.2.3 of the online Appendix

* The function [queueing.m](/App/Non-normal/Queueing/queueing.m) is a simulator of a two-node tandem queueing system;

* The function [truevalue.m](/App/Non-normal/Queueing/truevalue.m) estimates the mean performance of each worker allocation plan.

Get into folder [Queueing](/App/Non-normal/Queueing). Run [EAm.m](/main/Exp/EAm.m) for *EA*, [OCBAmsa.m](/main/Exp/OCBAmsa.m) for *OCBAm*, [OCBAmjia.m](/main/Exp/OCBAmjia.m) for *OCBAm+*, [OCBAss.m](/main/Exp/OCBAss.m) for *OCBAss*, [OCBASSS.m](/main/Exp/OCBASSS.m) for *OCBASS*, and [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*.

* Set corresponding input parameters for all policies;
* The input parameter *truemu* of each policy is calculated by running [truevalue.m](/App/Non-normal/Queueing/truevalue.m), which calls [queueing.m](/App/Non-normal/Queueing/queueing.m) during execution.

#### Numerical Experiments in Section A.4.3 of the online Appendix

* The function [pSAR.m](/App/Top-m%20arms%20identification/pSAR.m) is a compared policy from Bubeck et al. (2013);

* The function [pSR.m](/App/Top-m%20arms%20identification/pSR.m) is a compared policy from Audibert et al. (2010) and is then modified by Bubeck et al. (2013);

* The function [pGapE.m](/App/Top-m%20arms%20identification/pGapE.m) is a compared policy from Gabillon et al. (2011) and is then modified by Bubeck et al. (2013).

Remark:
* The policies can adjust sampling distributions;
* Update parameters in the same way when comparing;
* The independent macro experiments can run in parallel.

Get into folder [Top-m arms identification](/App/Top-m%20arms%20identification). Run [EAm.m](/main/Exp/EAm.m) for *EA*, [AOAPm.m](/main/Exp/AOAPm.m) for *AOAm*, [pGapE.m](/App/Top-m%20arms%20identification/pGapE.m) for *Gap-E*, [pSAR.m](/App/Top-m%20arms%20identification/pSAR.m) for *SAR*, [pSR.m](/App/Top-m%20arms%20identification/pSR.m) for *SR*.

* Set corresponding input parameters for all policies.

## Results

* The folder [Results](Results) contains the results of each figure in the numerical experiments in the main body of the paper.
* Due to the random nature of the simulation experiments and the absence of fixed random seeds in the code, the results may vary slightly when reproducing the experiment.
* Each row of the matrix corresponds to the numerical results of *EA*, *OCBAm*, *OCBAm+*, *OCBAss*, *OCBASS*, and *AOAm*, respectively.

## References

* Chen CH, He D, Fu M, Lee LH (2008) Efficient Simulation Budget al.location for Selecting an Optimal Subset. INFORMS Journal on Computing 20(4): 579–595.
* Zhang S, Lee LH, Chew EP, Chen CH, Jen HY (2012) An Improved Simulation Budget al.location Procedure to Efficiently Select the Optimal Subset of Many Alternatives. 2012 IEEE International Conference on Automation Science and Engineering (CASE), 230–236 (IEEE).
* Zhang S, Lee LH, Chew EP, Xu J, Chen CH (2015) A Simulation Budget al.location Procedure for Enhancing the Efficiency of Optimal Subset Selection. IEEE Transactions on Automatic Control 61(1): 62–75.
* Gao S, Chen W (2015) A Note on the Subset Selection for Simulation Optimization. Proceedings of the 2015 Winter Simulation Conference, 3768–3776 (IEEE).
* Gao S, Chen W (2016) A New Budget al.location Framework for Selecting Top Simulated Designs. IIE Transactions 48(9): 855–863.
* Bubeck S, Wang T, Viswanathan N. Multiple identifications in multi-armed bandits[C] // In International Conference on Machine Learning (ICML). PMLR, 2013: 258-265.
* Gabillon V, Ghavamzadeh M, Lazaric A, et al. Multi-bandit best arm identification[J]. In Advances in Neural Information Processing Systems (NIPS), 2011, 24.
* Audibert J Y, Bubeck S, Munos R. Best arm identification in multi-armed bandits[C]// In Proceedings of the 23rd Annual Conference on Learning Theory (COLT). Citeseer, 2010: 41-53.
