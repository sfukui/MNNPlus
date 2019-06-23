// <copyright file="ConstrainedHMC.cs" company="Fukui Shogo">
// Math.NET Numerics Appendix Library License (MIT/X11)
// ===================================================
// 
// Copyright (c) 2019 Fukui Shogo
// 
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// The copyright of original Math.Net Numerics is as follows.
//
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
//
// Copyright (c) 2009-2010 Math.NET
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using System;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;

namespace MathNet.Numerics.Appendix.Statistics.Mcmc
{
    /// <summary>
    /// A hybrid Monte Carlo sampler for multivariate distributions.
    /// </summary>
    public class ConstrainedHMC : ConstrainedHMCGeneric<double[]>
    {
        /// <summary>
        /// Number of parameters in the density function.
        /// </summary>
        private readonly int _length;

        /// <summary>
        /// Distribution to sample momentum from.
        /// </summary>
        private Normal _pDistribution;

        /// <summary>
        /// Standard deviations used in the sampling of different components of the
        /// momentum.
        /// </summary>
        private double[] _mpSdv;

        /// <summary>
        /// Gets or sets the standard deviations used in the sampling of different components of the
        /// momentum.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">When the length of pSdv is not the same as Length.</exception>
        public double[] MomentumStdDev
        {
            get { return (double[])_mpSdv.Clone(); }
            set
            {
                CheckVariance(value);
                _mpSdv = (double[])value.Clone();
            }
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviation
        /// 1 using the default <see cref="System.Random"/> random
        /// number generator. A three point estimation will be used for differentiation.
        /// This constructor will set the burn interval.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval = 0)
            : this(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, new double[x0.Length], SystemRandomSource.Default, Grad, Enumerable.Repeat(Double.NegativeInfinity, x0.Length).ToArray(), Enumerable.Repeat(Double.PositiveInfinity, x0.Length).ToArray())
        {
            for (int i = 0; i < _length; i++)
            {
                _mpSdv[i] = 1;
            }
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviation
        /// 1 using the default <see cref="System.Random"/> random
        /// number generator. A three point estimation will be used for differentiation.
        /// This constructor will set the burn interval.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <param name="xInfimums">The infimums of the sample values.</param>
        /// <param name="xSupremums">The supremums of the sample values.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval, double[] xInfimums, double[] xSupremums)
            : this(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, new double[x0.Length], SystemRandomSource.Default, Grad, xInfimums, xSupremums)
        {
            for (int i = 0; i < _length; i++)
            {
                _mpSdv[i] = 1;
            }
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviation
        /// specified by pSdv using the default <see cref="System.Random"/> random
        /// number generator. A three point estimation will be used for differentiation.
        /// This constructor will set the burn interval.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <param name="pSdv">The standard deviations of the normal distributions that are used to sample
        /// the components of the momentum.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval, double[] pSdv)
            : this(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, pSdv, SystemRandomSource.Default)
        {
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviation
        /// specified by pSdv using the default <see cref="System.Random"/> random
        /// number generator. A three point estimation will be used for differentiation.
        /// This constructor will set the burn interval.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <param name="pSdv">The standard deviations of the normal distributions that are used to sample
        /// <param name="xInfimums">The infimums of the sample values.</param>
        /// <param name="xSupremums">The supremums of the sample values.</param>
        /// the components of the momentum.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval, double[] pSdv, double[] xInfimums, double[] xSupremums)
            : this(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, pSdv, SystemRandomSource.Default, Grad, xInfimums, xSupremums)
        {
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviation
        /// specified by pSdv using the a random number generator provided by the user.
        /// A three point estimation will be used for differentiation.
        /// This constructor will set the burn interval.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <param name="pSdv">The standard deviations of the normal distributions that are used to sample
        /// the components of the momentum.</param>
        /// <param name="randomSource">Random number generator used for sampling the momentum.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval, double[] pSdv, System.Random randomSource)
            : this(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, pSdv, randomSource, Grad, Enumerable.Repeat(Double.NegativeInfinity, x0.Length).ToArray(), Enumerable.Repeat(Double.PositiveInfinity, x0.Length).ToArray())
        {
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviation
        /// specified by pSdv using the a random number generator provided by the user.
        /// A three point estimation will be used for differentiation.
        /// This constructor will set the burn interval.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <param name="pSdv">The standard deviations of the normal distributions that are used to sample
        /// the components of the momentum.</param>
        /// <param name="randomSource">Random number generator used for sampling the momentum.</param>
        /// <param name="xInfimums">The infimums of the sample values.</param>
        /// <param name="xSupremums">The supremums of the sample values.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval, double[] pSdv, System.Random randomSource, double[] xInfimums, double[] xSupremums)
            : this(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, pSdv, randomSource, Grad, xInfimums, xSupremums)
        {
        }

        /// <summary>
        /// Constructs a new Hybrid Monte Carlo sampler for a multivariate probability distribution.
        /// The components of the momentum will be sampled from a normal distribution with standard deviations
        /// given by pSdv. This constructor will set the burn interval, the method used for
        /// numerical differentiation and the random number generator.
        /// </summary>
        /// <param name="x0">The initial sample.</param>
        /// <param name="pdfLnP">The log density of the distribution we want to sample from.</param>
        /// <param name="frogLeapSteps">Number frog leap simulation steps.</param>
        /// <param name="stepSize">Size of the frog leap simulation steps.</param>
        /// <param name="burnInterval">The number of iterations in between returning samples.</param>
        /// <param name="pSdv">The standard deviations of the normal distributions that are used to sample
        /// the components of the momentum.</param>
        /// <param name="randomSource">Random number generator used for sampling the momentum.</param>
        /// <param name="diff">The method used for numerical differentiation.</param>
        /// <param name="xInfimums">The infimums of the sample values.</param>
        /// <param name="xSupremums">The supremums of the sample values.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the number of burnInterval iteration is negative.</exception>
        /// <exception cref="ArgumentOutOfRangeException">When the length of pSdv is not the same as x0.</exception>
        /// <exception cref="ArgumentException">When the all of the lengths of x0, xInfimums, and xSupremums are not same.</exception>
        /// <exception cref="ArgumentException">When one or more infimums are greater than or equal to the paired supremums.</exception>
        /// <exception cref="ArgumentException">When one or more infimums and/or supremums are NaN.</exception>
        public ConstrainedHMC(double[] x0, Numerics.Statistics.Mcmc.DensityLn<double[]> pdfLnP, int frogLeapSteps, double stepSize, int burnInterval, double[] pSdv, System.Random randomSource, DiffMethod diff, double[] xInfimums, double[] xSupremums)
            : base(x0, pdfLnP, frogLeapSteps, stepSize, burnInterval, randomSource, diff, xInfimums, xSupremums)
        {
            if (x0.Length != xInfimums.Length || x0.Length != xSupremums.Length)
                throw new ArgumentException("All of the lengths of initial samples, the infimums, and the supremums are not same.");

            for(int i = 0; i < x0.Length; i++)
            {
                if (xInfimums[i] >= xSupremums[i])
                    throw new ArgumentException("Some infimums are greater than or equal to the paired supremums.");
                if (Double.IsNaN(xInfimums[i]) || Double.IsNaN(xSupremums[i]))
                    throw new ArgumentException("Some infimums and/or supremums are NaN.");
            }

            _length = x0.Length;
            MomentumStdDev = pSdv;

            Initialize(x0);

            Burn(BurnInterval);
        }

        /// <summary>
        /// Initialize parameters.
        /// </summary>
        /// <param name="x0">The current location of the sampler.</param>
        private void Initialize(double[] x0)
        {
            Current = (double[])x0.Clone();
            _pDistribution = new Normal(0.0, 1.0, RandomSource);
        }

        /// <summary>
        /// Checking that the location and the momentum are of the same dimension and that each component is positive.
        /// </summary>
        /// <param name="pSdv">The standard deviations used for sampling the momentum.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the length of pSdv is not the same as Length or if any
        /// component is negative.</exception>
        /// <exception cref="ArgumentNullException">When pSdv is null.</exception>
        private void CheckVariance(double[] pSdv)
        {
            if (pSdv == null)
            {
                throw new ArgumentNullException(nameof(pSdv), "Standard deviation cannot be null.");
            }
            if (pSdv.Length != _length)
            {
                throw new ArgumentOutOfRangeException(nameof(pSdv), "Standard deviation of momentum must have same length as sample.");
            }
            if (pSdv.Any(sdv => sdv < 0))
            {
                throw new ArgumentOutOfRangeException(nameof(pSdv), "Standard deviation must be positive.");
            }
        }

        /// <summary>
        /// Use for copying objects in the Burn method.
        /// </summary>
        /// <param name="source">The source of copying.</param>
        /// <returns>A copy of the source object.</returns>
        protected override double[] Copy(double[] source)
        {
            var destination = new double[_length];
            Array.Copy(source, 0, destination, 0, _length);
            return destination;
        }

        /// <summary>
        /// Use for creating temporary objects in the Burn method.
        /// </summary>
        /// <returns>An object of type T.</returns>
        protected override double[] Create()
        {
            return new double[_length];
        }

        ///<inheritdoc/>
        protected override void DoAdd(ref double[] first, double factor, double[] second)
        {
            for (int i = 0; i < _length; i++)
            {
                first[i] += factor * second[i];
            }
        }

        /// <inheritdoc/>
        protected override void DoSubtract(ref double[] first, double factor, double[] second)
        {
            for (int i = 0; i < _length; i++)
            {
                first[i] -= factor * second[i];
            }
        }

        /// <inheritdoc/>
        protected override double DoProduct(double[] first, double[] second)
        {
            double prod = 0;
            for (int i = 0; i < _length; i++)
            {
                prod += first[i] * second[i];
            }
            return prod;
        }

        /// <inheritdoc/>
        protected override bool IsSampleWithinBounds(double[] xValues)
        {
            for(int i = 0; i < xValues.Length; i++)
            {
                if (xValues[i] <= _xInfimums[i] || xValues[i] >= _xSupremums[i] || Double.IsNaN(xValues[i]))
                    return false;
            }
            return true;
        }

        /// <inheritdoc/>
        protected override bool IsGradientValid(double[] gradients)
        {
            foreach(var g in gradients)
            {
                if (Double.IsInfinity(g) || Double.IsNaN(g))
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Samples the momentum from a normal distribution.
        /// </summary>
        /// <param name="p">The momentum to be randomized.</param>
        protected override void RandomizeMomentum(ref double[] p)
        {
            for (int j = 0; j < _length; j++)
            {
                p[j] = _mpSdv[j] * _pDistribution.Sample();
            }
        }

        /// <summary>
        /// The default method used for computing the gradient. Uses a simple three point estimation.
        /// </summary>
        /// <param name="function">Function which the gradient is to be evaluated.</param>
        /// <param name="x">The location where the gradient is to be evaluated.</param>
        /// <returns>The gradient of the function at the point x.</returns>
        static double[] Grad(Numerics.Statistics.Mcmc.DensityLn<double[]> function, double[] x)
        {
            int length = x.Length;
            var returnValue = new double[length];
            var increment = new double[length];
            var decrement = new double[length];

            Array.Copy(x, 0, increment, 0, length);
            Array.Copy(x, 0, decrement, 0, length);

            for (int i = 0; i < length; i++)
            {
                double y = x[i];
                double h = Math.Max(10e-4, (10e-7) * y);
                increment[i] += h;
                decrement[i] -= h;
                returnValue[i] = (function(increment) - function(decrement)) / (2 * h);
                increment[i] = y;
                decrement[i] = y;
            }

            return returnValue;
        }
    }
}


