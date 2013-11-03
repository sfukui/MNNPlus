MNNPlus
=================

MNNPlus is an extension of opensource numerical library "Math.NET Numerics"(The explapation of Math.NET Numerics is in [README_Original.md](README_Original.md) in detail).

MNNPlus provides some features frequently requiered in econometrics,
finance, and other numerical analysis in economics.

This fork append following features to original library "Math.Net Numerics".
* Numerical differentiation(Not include Hessian calculation yet).
* Numerical optimization(BFGS and Nelder-Mead method)
* MCMC sampling(Adaptive Rejection Metropolis Sampling)

The core part of MNNPlus is written in F#.
Of course we are able to use MNNPlus from C#through wrapper classes.

Currently, MNNPlus supports .NET 4.0/4.5/4.5.1 on Windows only.
We never guarantee the valid operation in the other systems.

MNNPlus is covererd under terms of the MIT/X11 license. 
You may therefore link to it and use it in both opensource and proprietary
software projects. See also the [COPYRIGHT](COPYRIGHT.markdown) file in the root folder.

Installation Instructions
-------------------------

First, download the *MathNet.Numerics.Appendix.FSharp.dll*(core library) and *MathNet.Numerics.Appendix.dll*(library of the wrapper classes) assemblies.
Then, add references to them in your project.

Author
-------------------
Fukui Shogo

Contributors of original MathNet.Numericws library are listed in [CONTRIBUTORS.md](CONTRIBUTORS.md).

