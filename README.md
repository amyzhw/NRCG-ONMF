## Introduction

Orthogonal Non-negative Matrix Factorization (ONMF) approximates a data matrix X by the product of two lower-dimensional factor matrices: X=UVT, with one of them orthogonal. ONMF has been widely applied for clustering, but it often suffers from high computational cost due to the orthogonality constraint. In this paper, we propose a method, called Nonlinear Riemannian Conjugate Gradient ONMF (NRCG-ONMF), which updates U and V alternatively and preserves the orthogonality of U while achieving fast convergence speed. Specifically, in order to update U, we develop a Nonlinear Riemannian Conjugate Gradient (NRCG) method on the Stiefel manifold using Barzilai-Borwein (BB) step size. For updating V, we use a closed-form solution under non-negativity constraint. Extensive experiments on both synthetic and real-world data sets show consistent superiority of our method over other approaches in terms of orthogonality preservation, convergence speed and clustering performance.

This work has been accepted as a full paper by CIKM 2016. 



## References
- [*Efficient Orthogonal Non-negative Matrix Factorization over Stiefel Manifold*](https://dl.acm.org/citation.cfm?doid=2983323.2983761). Wei Emma Zhang, Mingkui Tan, Quan Z. Sheng, Lina Yao and Qinfeng Shi. The 25th ACM International Conference on Information and Knowledge Management (CIKM 2017), Indianapolis, USA, October 2016. [bibTex](https://dblp.uni-trier.de/rec/bibtex/conf/cikm/ZhangTSYS16)

Please cite this paper when you use the codes. Thanks!