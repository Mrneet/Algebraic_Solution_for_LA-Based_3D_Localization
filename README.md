# Algebraic_Solution_for_LA-Based_3D_Localization

**Abstract**

Localizing a three-dimensional (3D) source using linear arrays (LAs) is a promising new localization technology. Existing solutions are either designed for specific LA deployments, are computationally intensive, or rely on iterative methods that do not guarantee convergence. This paper presents a novel algebraic solution algorithm for 3D source localization using space angle (SA) measurements from LAs. We propose a new formulation of the SA measurement equation, which leads to a constrained weighted least squares (CWLS) problem. Solving it by Lagrangian multipliers, the
optimal estimation is obtained with an error correction. The solution does not require specific arrangement and placement of LAs and effectively balances accuracy with computational efficiency. We analyze the performance and complexity of the proposed solution, demonstrating its ability to achieve the Cramér-Rao Lower Bound (CRLB) in the small error region under Gaussian noise with a low computational load. Simulations validate the analysis and confirm the superiority of the proposed solution compared to existing ones.

---

**If you use any of the following codes in your research, please cite the paper as a reference in your publication. Thank you!**

## Algebraic Solution for Linear Array-Based 3D Localization Without Deployment Limitations (IEEE Signal Processing Letters)

### <u>Reference</u>
>C. Li, B. Tang, Y. Yang, L. Chen and Y. Sun, "Algebraic Solution for Linear Array-Based 3D Localization Without Deployment Limitations," in IEEE Signal Processing Letters, vol. 32, pp. 1326-1330, 2025.

### Code List:
- Algebraic Solution:
  - LA3DLoc_CF.m
- Cramer Rao Bound: SA3DLocLACRLB.m
- Example: main_AOALocLA_CF_noise.m
