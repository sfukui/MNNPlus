// Math.NET Numerics Appendix Library License (MIT/X11)
// ===================================================
// 
// Copyright (c) 2013 Fukui Shogo
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

using System;
using System.Collections.Generic;
using MathNet.Numerics.Properties;
using MathNet.Numerics.Random;

namespace MathNet.Numerics.Distributions
{
    /// <summary>
    /// Generalized Beta Distribution of 2nd Kind.
    /// </summary>
    /// <remarks>
    /// <para>The generalized beta distribution of second kind(GB2) has three shape parameters
    /// (alpha, p, and q) and one scale parameter(beta). </para>
    /// <para>Random numbers for GB2 are generated with slice sampling. 
    /// In addition, the <see cref="System.Random"/> is used in slice sampling by default. 
    /// Users can get/set the random number generator by using the <see cref="RandomSource"/> property.</para>
    /// </remarks> 
    public class GeneralizedBeta2 : IContinuousDistribution
    {
        System.Random _random;

        double _alpha;
        double _beta;
        double _p;
        double _q;

        /// <summary>
        /// Initializes a new instance of generalized beta distribution of 2nd kind.
        /// </summary>
        /// <param name="alpha">One of the shape parameters "alpha". Range: alpha > 0.</param>
        /// <param name="beta">The scale parameter "beta". Range: beta > 0.</param>
        /// <param name="p">One of the shape parameters "p". Range: p > 0.</param>
        /// <param name="q">One of the shape parameters "q". Range: q > 0.</param>
        public GeneralizedBeta2(double alpha, double beta, double p, double q)
        {
            SetParameters(alpha, beta, p, q);
            _random = MersenneTwister.Default;
        }

        /// <summary>
        /// Initializes a new instance of generalized beta distribution of 2nd kind.
        /// </summary>
        /// <param name="theta">Parameters "alpha", "beta", "p", and "q". Range: alpha > 0, beta > 0, p > 0, and q > 0.</param>
        public GeneralizedBeta2(double[] theta)
        {
            if (theta.Length != 4)
            {
                throw new ArgumentException(Resources.ArgumentArrayWrongLength);
            }

            SetParameters(theta[0], theta[1], theta[2], theta[3]);
            _random = MersenneTwister.Default;
        }

        /// <summary>
        /// Initializes a new instance of generalized beta distribution of 2nd kind.
        /// </summary>
        /// <param name="theta">Parameters "alpha", "beta", "p", and "q". Range: alpha > 0, beta > 0, p > 0, and q > 0.</param>
        public GeneralizedBeta2(MathNet.Numerics.LinearAlgebra.Vector<double> theta)
        {
            if (theta.Count != 4)
            {
                throw new ArgumentException(Resources.ArgumentArrayWrongLength);
            }

            SetParameters(theta[0], theta[1], theta[2], theta[3]);
            _random = MersenneTwister.Default;
        }

        /// <summary>
        /// Sets the parameters of the distribution after checking their validity.
        /// </summary>
        /// <param name="alpha">One of the shape parameters "alpha". Range: alpha > 0.</param>
        /// <param name="beta">The scale parameter "beta". Range: beta > 0.</param>
        /// <param name="p">One of the shape parameters "p". Range: p > 0.</param>
        /// <param name="q">One of the shape parameters "q". Range: q > 0.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the parameters are out of range.</exception>
        void SetParameters(double alpha, double beta, double p, double q)
        {
            if (alpha <= 0.0 || beta <= 0.0 || p <= 0.0 || q <= 0.0 ||
                Double.IsNaN(alpha) || Double.IsNaN(beta) || Double.IsNaN(p) || Double.IsNaN(q) )
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            _alpha = alpha;
            _beta = beta;
            _p = p;
            _q = q;
        }

        /// <summary>
        /// A string representation of the distribution.
        /// </summary>
        /// <returns>a string representation of the distribution.</returns>
        public override string ToString()
        {
            return "GB2(alpha = " + _alpha + ", beta = " + _beta + ", p = " + _p + ", q = " + _q + ")";
        }

        /// <summary>
        /// Gets or sets the parameter "alpha" of the generalized beta distribution of second kind. Range: alpha > 0.
        /// </summary>
        public double Alpha
        {
            get
            {
                return _alpha;
            }
            set
            {
                SetParameters(value, _beta, _p, _q);
            }
        }

        /// <summary>
        /// Gets or sets the parameter "beta" of the generalized beta distribution of second kind. Range: beta > 0.
        /// </summary>
        public double Beta
        {
            get
            {
                return _beta;
            }
            set
            {
                SetParameters(_alpha, value, _p, _q);
            }
        }

        /// <summary>
        /// Gets or sets the parameter "p" of the generalized beta distribution of second kind. Range: p > 0.
        /// </summary>
        public double P
        { 
            get
            {
                return _p;
            }
            set
            {
                SetParameters(_alpha, _beta, value, _q);
            }
        }

        /// <summary>
        /// Gets or sets the parameter "q" of the generalized beta distribution of second kind. Range: q > 0.
        /// </summary>
        public double Q
        {
            get
            {
                return _q;
            }
            set
            {
                SetParameters(_alpha, _beta, _p, value);
            }
        }

        public double[] Parameters
        {
            get
            {
                return new double[4] { _alpha, _beta, _p, _q };
            }
            set
            {
                double[] parameters = value;
                SetParameters(parameters[0], parameters[1], parameters[2], parameters[3]);
            }
        }

        /// <summary>
        /// Gets the mode of the generalized beta distribution of second kind.
        /// </summary>
        double IContinuousDistribution.Mode
        {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        /// Gets the median of the generalized beta distribution of second kind.
        /// </summary>
        double IContinuousDistribution.Median
        {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        /// Gets the minimum of the generalized beta distribution of second kind.
        /// </summary>
        public double Minimum
        {
            get
            {
                return 0.0;
            }
        }

        /// <summary>
        /// Gets the maximum of the generalized beta distribution of second kind.
        /// </summary>
        public double Maximum
        {
            get
            {
                return Double.PositiveInfinity;
            }
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x with given parameters.
        /// The value is calculated as exp(DensityLn(x)).
        /// </summary>
        /// <param name="x">The location at which to compute the density.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="PDF"/>
        public double Density(double x)
        {
            return System.Math.Exp(DensityLn(x));
        }

        /// <summary>
        /// Computes the log probability density of the distribution (lnPDF) at x with given parameters.
        /// The functional form is as follows:
        /// ln(alpha) + (alpha * p - 1.0) * ln(x) - (alpha * p) * ln(beta) - BetaLn(p,q) - (p + q) * ln(1.0 + (x / beta)^alpha)
        /// </summary>
        /// <param name="x">The location at which to compute the log density.</param>
        /// <returns>the log density at <paramref name="x"/>.</returns>
        /// <seealso cref="PDFLn"/>
        public double DensityLn(double x)
        {
            return PDFLn(_alpha, _beta, _p, _q, x);
        }

        /// <summary>
        /// Generates a sample from the generalized beta distribution of second kind.
        /// The method is not implemented at this time.
        /// If you want to sample from this distribution, please use some MCMC method.
        /// </summary>
        /// <returns>a sample from the distribution.</returns>
        double IContinuousDistribution.Sample()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Generates a sequence of samples from the Gamma distribution.
        /// The method is not implemented at this time.
        /// If you want to sample from this distribution, please use some MCMC method.
        /// </summary>
        /// <returns>a sequence of samples from the distribution.</returns>
        IEnumerable<double> IContinuousDistribution.Samples()
        {
            throw new NotImplementedException();
        }

        double IUnivariateDistribution.Mean
        {
            get { throw new NotImplementedException(); }
        }

        double IUnivariateDistribution.Variance
        {
            get { throw new NotImplementedException(); }
        }

        double IUnivariateDistribution.StdDev
        {
            get { throw new NotImplementedException(); }
        }

        double IUnivariateDistribution.Entropy
        {
            get { throw new NotImplementedException(); }
        }

        double IUnivariateDistribution.Skewness
        {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        /// Gets the mean of the generalized beta distribution of second kind.
        /// The value is derived with numerical integration.
        /// </summary>
        /// <param name="xmax">Upper limit of numerical integration.</param>
        /// <returns>Mean of the distribution.</returns>
        public double Mean(double xmax)
        {
            return Integration.DoubleExponentialTransformation.Integrate(y => y * this.Density(y), 0.0, xmax, 0.001);
        }

        /// <summary>
        /// Gets the variance of the generalized beta distribution of second kind.
        /// The value is derived with numerical integration.
        /// </summary>
        /// <param name="xmax">Upper limit of numerical integration.</param>
        /// <returns>Variance of the distribution.</returns>
        public double Variance(double xmax)
        {
            double mean = this.Mean(xmax);
            return Integration.DoubleExponentialTransformation.Integrate(y => (y - mean) * (y - mean) * this.Density(y), 0.0, xmax, 0.001);
        }

        /// <summary>
        /// Gets the standard deviation of the generalized beta distribution of second kind.
        /// The value is calculated as square root of Variance(xmax).
        /// </summary>
        /// <param name="xmax">Upper limit of numerical integration.</param>
        /// <returns>Standard deviation of the distribution.</returns>
        public double StdDev(double xmax)
        {
            return System.Math.Sqrt(Variance(xmax));
        }

        /// <summary>
        /// Computes the cumulative distribution (CDF) of the distribution at x with given parameters.
        /// The value is derived from numerical integration.
        /// </summary>
        /// <param name="x">The location at which to compute the cumulative distribution function.</param>
        /// <returns>The cumulative distribution at location <paramref name="x"/>.</returns>
        /// <seealso cref="CDF"/>
        public double CumulativeDistribution(double x)
        {
            return CDF(_alpha, _beta, _p, _q, x);
        }

        /// <summary>
        /// Gets the gini coefficient of the generalized beta distribution of second kind.
        /// The value is derived with numerical integration.
        /// </summary>
        /// <param name="xmax">Upper limit of numerical integration.</param>
        /// <returns>Gini coefficient of the distribution.</returns>
        public double Gini(double xmax)
        {
            double mean = Mean(xmax);
            return 1.0 - (1.0 / mean) *
                Integration.DoubleExponentialTransformation.Integrate(y => System.Math.Pow((1.0 - this.CumulativeDistribution(y)), 2.0), 0.0, xmax, 0.001);
        }

        /// <summary>
        /// Gets the Theil measure of the generalized beta distribution of second kind.
        /// The value is derived with numerical integration.
        /// </summary>
        /// <param name="xmax">Upper limit of numerical integration.</param>
        /// <returns>Theil index of the distribution.</returns>
        public double Theil(double xmax)
        {
            double mean = this.Mean(xmax);
            return Integration.DoubleExponentialTransformation.Integrate(y => (y / mean) * System.Math.Log(y / mean) * this.Density(y), 0.0, xmax, 0.001);
        }

        System.Random IDistribution.RandomSource
        {
            get
            {
                throw new NotImplementedException();
            }
            set
            {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x with parameters given by arguments;
        /// alpha, beta, p, and q.
        /// The functional form is as follows:
        /// ln(alpha) + (alpha * p - 1.0) * ln(x) - (alpha * p) * ln(beta) - BetaLn(p,q) - (p + q) * ln(1.0 + (x / beta)^alpha)
        /// </summary>
        /// <param name="alpha">One of the shape parameters "alpha". Range: alpha > 0.</param>
        /// <param name="beta">The scale parameter "beta". Range: beta > 0.</param>
        /// <param name="p">One of the shape parameters "p". Range: p > 0.</param>
        /// <param name="q">One of the shape parameters "q". Range: q > 0.</param>
        /// <param name="x">The location at which to compute the density.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="DensityLn"/>
        public static double PDFLn(double alpha, double beta, double p, double q, double x)
        {
            double[] theta = new double[4] { alpha, beta, p, q };

            return PDFLn(theta, x);
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x with parameters given by arguments;
        /// alpha, beta, p, and q.
        /// The value is calculated as exp(PDFLn(alpha, beta, p, q, x)).
        /// </summary>
        /// <param name="alpha">One of the shape parameters "alpha". Range: alpha > 0.</param>
        /// <param name="beta">The scale parameter "beta". Range: beta > 0.</param>
        /// <param name="p">One of the shape parameters "p". Range: p > 0.</param>
        /// <param name="q">One of the shape parameters "q". Range: q > 0.</param>
        /// <param name="x">The location at which to compute the density.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="Density"/>
        public static double PDF(double alpha, double beta, double p, double q, double x)
        {
            return System.Math.Exp(PDFLn(alpha, beta, p, q, x));
        }

        /// <summary>
        /// Computes the cumulative distribution (CDF) of the distribution at x with parameters given by arguments;
        /// alpha, beta, p, and q.
        /// The value is derived from numerical integration.
        /// </summary>
        /// <param name="alpha">One of the shape parameters "alpha". Range: alpha > 0.</param>
        /// <param name="beta">The scale parameter "beta". Range: beta > 0.</param>
        /// <param name="p">One of the shape parameters "p". Range: p > 0.</param>
        /// <param name="q">One of the shape parameters "q". Range: q > 0.</param>
        /// <param name="x">The location at which to compute the cumulative distribution function.</param>
        /// <returns>The cumulative distribution at location <paramref name="x"/>.</returns>
        /// <seealso cref="CumulativeDistribution"/>
        public static double CDF(double alpha, double beta, double p, double q, double x)
        {
            return Integration.DoubleExponentialTransformation.Integrate(y => PDF(alpha, beta, p, q, y), 0.0, x, 0.001);
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x with parameters given by the argument "theta".
        /// The functional form is as follows:
        /// ln(theta[0]) + (theta[0] * theta[2] - 1.0) * ln(x) - (theta[0] * theta[2]) * ln(theta[1]) 
        /// - BetaLn(theta[2],theta[3]) - (theta[2] + theta[3]) * ln(1.0 + (x / theta[1])^theta[0])
        /// </summary>
        /// <param name="theta">Parameters "alpha", "beta", "p", and "q". Range: alpha > 0, beta > 0, p > 0, and q > 0.</param>
        /// <param name="x">The location at which to compute the density.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="DensityLn"/>
        public static double PDFLn(double[] theta, double x)
        {
            if (theta.Length != 4)
            {
                throw new ArgumentException(Resources.ArgumentArrayWrongLength);
            }

            if (Array.Exists(theta, ot => ot <= 0.0) || Array.Exists(theta, ot => Double.IsNaN(ot)))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            double lnumer = System.Math.Log(theta[0]) + (theta[0] * theta[2] - 1.0) * System.Math.Log(x);
            double ldenom = (theta[0] * theta[2]) * System.Math.Log(theta[1])
                + SpecialFunctions.BetaLn(theta[2], theta[3])
                + (theta[2] + theta[3]) * System.Math.Log(1.0 + System.Math.Pow((x / theta[1]), theta[0]));

            return (lnumer - ldenom);
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x with parameters given by the argument "theta".
        /// The value is calculated as exp(PDFLn(theta, x)).
        /// </summary>
        /// <param name="theta">Parameters "alpha", "beta", "p", and "q". Range: alpha > 0, beta > 0, p > 0, and q > 0.</param>
        /// <param name="x">The location at which to compute the density.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="Density"/>
        public static double PDF(double[] theta, double x)
        {
            return System.Math.Exp(PDFLn(theta, x));
        }

        /// <summary>
        /// Computes the cumulative distribution (CDF) of the distribution at x with parameters given by the argument "theta".
        /// The value is derived from numerical integration.
        /// </summary>
        /// <param name="theta">Parameters "alpha", "beta", "p", and "q". Range: alpha > 0, beta > 0, p > 0, and q > 0.</param>
        /// <param name="x">The location at which to compute the cumulative distribution function.</param>
        /// <returns>The cumulative distribution at location <paramref name="x"/>.</returns>
        /// <seealso cref="CumulativeDistribution"/>
        public static double CDF(double[] theta, double x)
        {
            return Integration.DoubleExponentialTransformation.Integrate(y => PDF(theta, y), 0.0, x, 0.001);
        }
    }
}
