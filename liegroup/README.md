# Evaluate Lie group functions safely in Fortran
This project contains several modules that allow easy and safe evaluation of Lie group functions related to the following Lie groups of unit quaternions 𝕊³ and its semidirect product 𝕊³⋉ℝ³. In order to do so, there is a module that evaluates singular functions, minimizing round-off errors. Additionally, some functionality of quaternions are included.

## Contents

### Module `cross_functions`
Contains just three functions related to cross products:

 * `cross(v,w)`: Calculates the cross product 𝑣×𝑤.
 * `skw(v)`: Calculates the skew-symmetric matrix skw(𝑣) such that skw(𝑣)⋅𝑤=𝑣×𝑤.
 * `dynadic_product(x,y)`: Calculates the dyadic product 𝑥⋅𝑦ᵀ of column vectors 𝑥, 𝑦 of arbitrary (but equal) length.

### Module `quaternion_functions`
Contains functions related to quaternions in general. Quaternions are represented as 4D vectors.

 * `qp(p,q)`: Calculates the quaternion product 𝑝*𝑞.
 * `apply_quat(p,v)`: Applies the rotation represented by a (unit) quaternion 𝑝 to a 3D vector 𝑣.
 * `apply_conj_quat(p,v)`: Applies the inverse rotation represented by a (unit) quaternion 𝑝 to a 3D vector 𝑣. The name refers to the fact that the inverse of a unit quaternion is its conjugate.
 * `conj_quat(p)`: Calculate the conjugate of a quaternion.

### Module `s3_functions`
Contains functions related to the Lie group of unit quaternions 𝕊³:

 * `lp_s3(p1,p2)`: Calculates the Lie product 𝑝₁*𝑝₂ of two unit quaternions. This is just the quaternion product.
 * `inv_s3(p)`: Calculates the inverse of a unit quaternion. This is just the conjugate.
 * `lie_bracket_s3(v1,v2)`: Calculates the result of ̂𝑣₁⋅𝑣₂ with the hat operator. The name comes from the fact that the hat operator is defined by the Lie bracket of the Lie algebra 𝔰³ and the tilde operator chosen for 𝕊³. Also, the mapping (𝑣₁,𝑣₂)↦̂𝑣₁⋅𝑣₂ fulfills the definition of a Lie bracket on ℝ³.
 * `expt_s3(v)`: Calculates the exponential of the tilde operator expt(𝑣)=exp(̃𝑣).
 * `tan_op_s3(v)`: Calculates the tangent operator 𝐓(𝑣).
 * `tan_tr_mult_s3(v,w)`: Calculates the result of 𝐓ᵀ(𝑣)⋅𝑤.
 * `tan_op_inv_s3(v)`: Calculates the inverse matrix 𝐓⁻¹(𝑣)=(𝐓(𝑣))⁻¹ of the tangent operator.
 * `tan_tr_inv_s3(v)`: Calculates 𝐓⁻ᵀ(𝑣)=((𝐓(𝑣))ᵀ)⁻¹, the inverse matrix of the transposed tangent operator.
 * `tan_op_inv_mult_s3(v,w)`: Calculates 𝐓⁻¹(𝑣)⋅𝑤.
 * `tan_tr_inv_mult_s3(v,w)`: Calculates 𝐓⁻ᵀ(𝑣)⋅𝑤.
 * `logt_s3(p)`: Calculates the logarithm logt(𝑝) with results in ℝ³ of a quaternion 𝑝 sufficiently close to the identity quaternion. Note that it holds tilde(logt(𝑝))=log(𝑝).
 * `d_tan_tr_inv_s3(Om, A)`: Calculates the Jacobi matrix of 𝐓⁻ᵀ(Ω)⋅𝐴 with respect to the vector Ω∈ℝ³, where 𝐴∈ℝ³.

It is recommended to use the functions with `_mult` in the function name whereever possible for increased speed.

### Module `s3sdr3_functions`
Contains functions related to the semidirect product Lie group 𝕊³⋉ℝ³:

 * `lp_s3sdr3(q1,q2)`: Calculates the Lie product 𝑞₁*𝑞₂ for 𝑞₁,𝑞₂∈𝕊³⋉ℝ³.
 * `inv_s3sdr3(q)`: Calculates the inverse of an element 𝑞∈𝕊³⋉ℝ³.
 * `lie_bracket_s3(v1,v2)`: Calculates the result of ̂𝑣₁⋅𝑣₂ with the hat operator and 𝑣₁,𝑣₂∈ℝ⁶. The name is explained in the module `s3_functions`.
 * `hat_tr_s3sdr3(v)`: Calculates the 6×6-matrix ̂𝑣ᵀ.
 * `hat_tr_mult_s3sdr3(v,w)`: Calculates the result of ̂𝑣ᵀ⋅𝑤 for 𝑣,𝑤∈ℝ⁶.
 * `expt_s3sdr3(v)`: Calculates the exponential of the tilde operator expt(𝑣)=exp(̃𝑣).
 * `tan_op_s3sdr3(v)`: Calculates the tangent operator 𝐓(𝑣).
 * `tan_op_inv_s3sdr3(v)`: Calculates the inverse matrix 𝐓⁻¹(𝑣)=(𝐓(𝑣))⁻¹ of the tangent operator.
 * `tan_tr_inv_s3sdr3(v)`: Calculates 𝐓⁻ᵀ(𝑣)=((𝐓(𝑣))ᵀ)⁻¹, the inverse matrix of the transposed tangent operator.
 * `tan_tr_mult_s3sdr3(v,w)`: Calculates the result of 𝐓ᵀ(𝑣)⋅𝑤.
 * `tan_op_inv_mult_s3sdr3(v,w)`: Calculates 𝐓⁻¹(𝑣)⋅𝑤.
 * `tan_tr_inv_mult_s3sdr3(v,w)`: Calculates 𝐓⁻ᵀ(𝑣)⋅𝑤.
 * `logt_s3sdr3(q)`: Calculates the logarithm logt(𝑞) with results in ℝ⁶ of a 𝑞∈𝕊³⋉ℝ³ sufficiently close to the identity. Note that it holds tilde(logt(𝑞))=log(𝑞).
 * `d_tan_tr_inv_s3sdr3(v, w)`: Calculates the Jacobi matrix of 𝐓⁻ᵀ(Ω)⋅𝐴 with respect to the vector Ω∈ℝ³, where 𝐴∈ℝ³.

### Module `singular_functions`
This is a module that implements functions that are used in the other modules of this project, which have a singularity. In all of these functions, we assume that the input variable 𝑥 is nonnegative, although there will be no error if 𝑥 is actually negative. Usually, 𝑥 will be the norm of some vector. All function names try to describe the function it evaluates by stringing together all letters and numbers by underscores.
If 𝑥 is near a singular value, the function will be evaluated by a Taylor approximation. All thresholds and Taylor polynomial degrees are determined to minimize the errors as well as the loss-of-significance errors, which could occasionally become very large for 𝑥 very close to a singular value and no Taylor approximation would be used.
Note that there is a file `singular_functions_1e-6.F90`, which implements the same module but with higher precision. This will, however, result in more computing time.

### Testing
The file `test.F90` contains a bunch of tests. At this point in time, all tests should pass, except for the test of `d_tan_tr_inv_s3sdr3` which is slightly above the error threshold.

## Usage
Call `make` in order to build the project. Object and module files will be created in the directory `obj/`. Call `make test` in order to run the test.

## Validation and finding error bounds
In order to validate the implementations there is the Mathematica file `misc/Lie_group_functions.nb`. It implements most of the Lie group functions.
Finding the error bounds that are used in the module `singular_functions` is not easy. Here, the Mathematica file `misc/singular_coefficients.nb` was used. It defines a function `singAnal[f,OptionsPattern[]]`, where a function `f` can be passed. In the options pattern, the critical point, the threshold, the order to the Taylor polynomial (among other things) can be specified. Calling `singAnal` on a function will give an overview of the function near the critical value. For this, the function is evaluated as given with machine precision and with extremely high precision. Futhermore, the Taylor approximation given and is evaluated and the three are compared. The goal is to adjust the threshold and the order of the Taylor polynomial in such a way that neither the Taylor approximation is too imprecise inside the threshold nor the as-is-evaluation outside the threshold has loss of significance that is too big.

## Related projects
Integrators:

 * [The Lie group generalized-α method `gena`](https://github.com/StHante/gena)
 * [The Lie group BDF method `BLieDF`](https://github.com/StHante/BLieDF)
 * [The Lie group RATTLE method `RATTLie`](https://github.com/StHante/RATTLie)
 * [The Lie group SHAKE method `SHAKELie`](https://github.com/StHante/SHAKELie)
 * [The nonholonomic RATTLie method `RATTLie_nonhol`](https://github.com/StHante/RATTLie_nonhol)

Test problems:

 * [The heavy top example `heavy_top`](https://github.com/StHante/heavy_top)
 * [The constrained Cosserat beam model `crmS3R3`](https://github.com/StHante/crmS3R3)
 * [The rolling disk example `rolling_disk`](https://github.com/StHante/rolling_disk)

Miscellaneous:

 * [Implementation of Lie group functions `liegroup`](https://github.com/StHante/liegroup)
 * [Expand a config file with different configurations to several files `expandconfig`](https://github.com/StHante/expandconfig)
 * [Read lua files in Matlab and Octave `readLua`](https://github.com/StHante/readLua-for-Matlab-and-Octave)

Third party projects:

 * [Reading lua files in Fortran `aotus`](https://geb.sts.nt.uni-siegen.de/doxy/aotus/)
 * [GFortran](https://gcc.gnu.org/fortran/)
 * [GNU Parallel](https://www.gnu.org/software/parallel/)
