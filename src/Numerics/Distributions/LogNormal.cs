﻿// <copyright file="LogNormal.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
//
// Copyright (c) 2009-2013 Math.NET
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
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Properties;
using MathNet.Numerics.Random;
using MathNet.Numerics.Statistics;

namespace MathNet.Numerics.Distributions
{
    /// <summary>
    /// Continuous Univariate Log-Normal distribution.
    /// For details about this distribution, see 
    /// <a href="http://en.wikipedia.org/wiki/Log-normal_distribution">Wikipedia - Log-Normal distribution</a>.
    /// </summary>
    /// <remarks><para>The distribution will use the <see cref="System.Random"/> by default. 
    /// Users can get/set the random number generator by using the <see cref="RandomSource"/> property.</para>
    /// <para>The statistics classes will check all the incoming parameters whether they are in the allowed
    /// range. This might involve heavy computation. Optionally, by setting Control.CheckDistributionParameters
    /// to <c>false</c>, all parameter checks can be turned off.</para></remarks>
    public class LogNormal : IContinuousDistribution
    {
        System.Random _random;

        double _mu;
        double _sigma;

        /// <summary>
        /// Initializes a new instance of the <see cref="LogNormal"/> class. 
        /// The distribution will be initialized with the default <seealso cref="System.Random"/>
        /// random number generator.
        /// </summary>
        /// <param name="mu">The log-scale (μ) of the logarithm of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the logarithm of the distribution. Range: σ ≥ 0.</param>
        public LogNormal(double mu, double sigma)
        {
            _random = MersenneTwister.Default;
            SetParameters(mu, sigma);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogNormal"/> class. 
        /// The distribution will be initialized with the default <seealso cref="System.Random"/>
        /// random number generator.
        /// </summary>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <param name="randomSource">The random number generator which is used to draw random samples.</param>
        public LogNormal(double mu, double sigma, System.Random randomSource)
        {
            _random = randomSource ?? MersenneTwister.Default;
            SetParameters(mu, sigma);
        }

        /// <summary>
        /// Constructs a log-normal distribution with the desired mu and sigma parameters.
        /// </summary>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <param name="randomSource">The random number generator which is used to draw random samples. Optional, can be null.</param>
        /// <returns>A log-normal distribution.</returns>
        public static LogNormal WithMuSigma(double mu, double sigma, System.Random randomSource = null)
        {
            return new LogNormal(mu, sigma, randomSource);
        }

        /// <summary>
        /// Constructs a log-normal distribution with the desired mean and variance.
        /// </summary>
        /// <param name="mean">The mean of the log-normal distribution.</param>
        /// <param name="var">The variance of the log-normal distribution.</param>
        /// <param name="randomSource">The random number generator which is used to draw random samples. Optional, can be null.</param>
        /// <returns>A log-normal distribution.</returns>
        public static LogNormal WithMeanVariance(double mean, double var, System.Random randomSource = null)
        {
            var sigma2 = Math.Log(var/(mean*mean) + 1.0);
            return new LogNormal(Math.Log(mean) - sigma2 / 2.0, Math.Sqrt(sigma2), randomSource);
        }

        /// <summary>
        /// Estimates the log-normal distribution parameters from sample data with maximum-likelihood.
        /// </summary>
        /// <param name="samples">The samples to estimate the distribution parameters from.</param>
        /// <param name="randomSource">The random number generator which is used to draw random samples. Optional, can be null.</param>
        /// <returns>A log-normal distribution.</returns>
        /// <remarks>MATLAB: lognfit</remarks>
        public static LogNormal Estimate(IEnumerable<double> samples, System.Random randomSource = null)
        {
            var muSigma2 = samples.Select(s => Math.Log(s)).MeanVariance();
            return new LogNormal(muSigma2.Item1, Math.Sqrt(muSigma2.Item2), randomSource);
        }

        /// <summary>
        /// A string representation of the distribution.
        /// </summary>
        /// <returns>a string representation of the distribution.</returns>
        public override string ToString()
        {
            return "LogNormal(μ = " + _mu + ", σ = " + _sigma + ")";
        }

        /// <summary>
        /// Sets the parameters of the distribution after checking their validity.
        /// </summary>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the parameters are out of range.</exception>
        void SetParameters(double mu, double sigma)
        {
            if (sigma < 0.0 || Double.IsNaN(mu) || Double.IsNaN(sigma))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            _mu = mu;
            _sigma = sigma;
        }

        /// <summary>
        /// Gets or sets the log-scale (μ) (mean of the logarithm) of the distribution.
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set { SetParameters(value, _sigma); }
        }

        /// <summary>
        /// Gets or sets the shape (σ) (standard deviation of the logarithm) of the distribution. Range: σ ≥ 0.
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set { SetParameters(_mu, value); }
        }

        /// <summary>
        /// Gets or sets the random number generator which is used to draw random samples.
        /// </summary>
        public System.Random RandomSource
        {
            get { return _random; }
            set { _random = value ?? MersenneTwister.Default; }
        }

        /// <summary>
        /// Gets the mu of the log-normal distribution.
        /// </summary>
        public double Mean
        {
            get { return Math.Exp(_mu + (_sigma*_sigma/2.0)); }
        }

        /// <summary>
        /// Gets the variance of the log-normal distribution.
        /// </summary>
        public double Variance
        {
            get
            {
                var sigma2 = _sigma*_sigma;
                return (Math.Exp(sigma2) - 1.0)*Math.Exp(_mu + _mu + sigma2);
            }
        }

        /// <summary>
        /// Gets the standard deviation of the log-normal distribution.
        /// </summary>
        public double StdDev
        {
            get
            {
                var sigma2 = _sigma*_sigma;
                return Math.Sqrt((Math.Exp(sigma2) - 1.0)*Math.Exp(_mu + _mu + sigma2));
            }
        }

        /// <summary>
        /// Gets the entropy of the log-normal distribution.
        /// </summary>
        public double Entropy
        {
            get { return 0.5 + Math.Log(_sigma) + _mu + Constants.LogSqrt2Pi; }
        }

        /// <summary>
        /// Gets the skewness of the log-normal distribution.
        /// </summary>
        public double Skewness
        {
            get
            {
                var expsigma2 = Math.Exp(_sigma*_sigma);
                return (expsigma2 + 2.0)*Math.Sqrt(expsigma2 - 1);
            }
        }

        /// <summary>
        /// Gets the mode of the log-normal distribution.
        /// </summary>
        public double Mode
        {
            get { return Math.Exp(_mu - (_sigma*_sigma)); }
        }

        /// <summary>
        /// Gets the median of the log-normal distribution.
        /// </summary>
        public double Median
        {
            get { return Math.Exp(_mu); }
        }

        /// <summary>
        /// Gets the minimum of the log-normal distribution.
        /// </summary>
        public double Minimum
        {
            get { return 0.0; }
        }

        /// <summary>
        /// Gets the maximum of the log-normal distribution.
        /// </summary>
        public double Maximum
        {
            get { return Double.PositiveInfinity; }
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x, i.e. ∂P(X ≤ x)/∂x.
        /// </summary>
        /// <param name="x">The location at which to compute the density.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="PDF"/>
        public double Density(double x)
        {
            if (x < 0.0)
            {
                return 0.0;
            }

            var a = (Math.Log(x) - _mu) / _sigma;
            return Math.Exp(-0.5*a*a)/(x*_sigma*Constants.Sqrt2Pi);
        }

        /// <summary>
        /// Computes the log probability density of the distribution (lnPDF) at x, i.e. ln(∂P(X ≤ x)/∂x).
        /// </summary>
        /// <param name="x">The location at which to compute the log density.</param>
        /// <returns>the log density at <paramref name="x"/>.</returns>
        /// <seealso cref="PDFLn"/>
        public double DensityLn(double x)
        {
            if (x < 0.0)
            {
                return Double.NegativeInfinity;
            }

            var a = (Math.Log(x) - _mu) / _sigma;
            return (-0.5*a*a) - Math.Log(x*_sigma) - Constants.LogSqrt2Pi;
        }

        /// <summary>
        /// Computes the cumulative distribution (CDF) of the distribution at x, i.e. P(X ≤ x).
        /// </summary>
        /// <param name="x">The location at which to compute the cumulative distribution function.</param>
        /// <returns>the cumulative distribution at location <paramref name="x"/>.</returns>
        /// <seealso cref="CDF"/>
        public double CumulativeDistribution(double x)
        {
            return x < 0.0 ? 0.0
                : 0.5*SpecialFunctions.Erfc((_mu - Math.Log(x))/(_sigma*Constants.Sqrt2));
        }

        /// <summary>
        /// Computes the inverse of the cumulative distribution function (InvCDF) for the distribution
        /// at the given probability. This is also known as the quantile or percent point function.
        /// </summary>
        /// <param name="p">The location at which to compute the inverse cumulative density.</param>
        /// <returns>the inverse cumulative density at <paramref name="p"/>.</returns>
        /// <seealso cref="InvCDF"/>
        public double InverseCumulativeDistribution(double p)
        {
            return p <= 0.0 ? 0.0 : p >= 1.0 ? double.PositiveInfinity
                : Math.Exp(_mu - _sigma*Constants.Sqrt2*SpecialFunctions.ErfcInv(2.0*p));
        }

        /// <summary>
        /// Generates a sample from the log-normal distribution using the <i>Box-Muller</i> algorithm.
        /// </summary>
        /// <returns>a sample from the distribution.</returns>
        public double Sample()
        {
            return Math.Exp(Normal.Sample(_random, _mu, _sigma));
        }

        /// <summary>
        /// Generates a sequence of samples from the log-normal distribution using the <i>Box-Muller</i> algorithm.
        /// </summary>
        /// <returns>a sequence of samples from the distribution.</returns>
        public IEnumerable<double> Samples()
        {
            return Normal.Samples(_random, _mu, _sigma).Select(Math.Exp);
        }

        /// <summary>
        /// Computes the probability density of the distribution (PDF) at x, i.e. ∂P(X ≤ x)/∂x.
        /// </summary>
        /// <param name="x">The location at which to compute the density.</param>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <returns>the density at <paramref name="x"/>.</returns>
        /// <seealso cref="Density"/>
        /// <remarks>MATLAB: lognpdf</remarks>
        public static double PDF(double mu, double sigma, double x)
        {
            if (sigma < 0.0) throw new ArgumentOutOfRangeException("sigma", Resources.InvalidDistributionParameters);

            if (x < 0.0)
            {
                return 0.0;
            }

            var a = (Math.Log(x) - mu) / sigma;
            return Math.Exp(-0.5*a*a)/(x*sigma*Constants.Sqrt2Pi);
        }

        /// <summary>
        /// Computes the log probability density of the distribution (lnPDF) at x, i.e. ln(∂P(X ≤ x)/∂x).
        /// </summary>
        /// <param name="x">The location at which to compute the density.</param>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <returns>the log density at <paramref name="x"/>.</returns>
        /// <seealso cref="DensityLn"/>
        public static double PDFLn(double mu, double sigma, double x)
        {
            if (sigma < 0.0) throw new ArgumentOutOfRangeException("sigma", Resources.InvalidDistributionParameters);

            if (x < 0.0)
            {
                return Double.NegativeInfinity;
            }

            var a = (Math.Log(x) - mu) / sigma;
            return (-0.5*a*a) - Math.Log(x*sigma) - Constants.LogSqrt2Pi;
        }

        /// <summary>
        /// Computes the cumulative distribution (CDF) of the distribution at x, i.e. P(X ≤ x).
        /// </summary>
        /// <param name="x">The location at which to compute the cumulative distribution function.</param>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <returns>the cumulative distribution at location <paramref name="x"/>.</returns>
        /// <seealso cref="CumulativeDistribution"/>
        /// <remarks>MATLAB: logncdf</remarks>
        public static double CDF(double mu, double sigma, double x)
        {
            if (sigma < 0.0) throw new ArgumentOutOfRangeException("sigma", Resources.InvalidDistributionParameters);

            return x < 0.0 ? 0.0
                : 0.5*(1.0 + SpecialFunctions.Erf((Math.Log(x) - mu)/(sigma*Constants.Sqrt2)));
        }

        /// <summary>
        /// Computes the inverse of the cumulative distribution function (InvCDF) for the distribution
        /// at the given probability. This is also known as the quantile or percent point function.
        /// </summary>
        /// <param name="p">The location at which to compute the inverse cumulative density.</param>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <returns>the inverse cumulative density at <paramref name="p"/>.</returns>
        /// <seealso cref="InverseCumulativeDistribution"/>
        /// <remarks>MATLAB: logninv</remarks>
        public static double InvCDF(double mu, double sigma, double p)
        {
            if (sigma < 0.0) throw new ArgumentOutOfRangeException("sigma", Resources.InvalidDistributionParameters);

            return p <= 0.0 ? 0.0 : p >= 1.0 ? double.PositiveInfinity
                : Math.Exp(mu - sigma*Constants.Sqrt2*SpecialFunctions.ErfcInv(2.0*p));
        }

        /// <summary>
        /// Generates a sample from the log-normal distribution using the <i>Box-Muller</i> algorithm.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <returns>a sample from the distribution.</returns>
        public static double Sample(System.Random rnd, double mu, double sigma)
        {
            return Math.Exp(Normal.Sample(rnd, mu, sigma));
        }

        /// <summary>
        /// Generates a sequence of samples from the log-normal distribution using the <i>Box-Muller</i> algorithm.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="mu">The log-scale (μ) of the distribution.</param>
        /// <param name="sigma">The shape (σ) of the distribution. Range: σ ≥ 0.</param>
        /// <returns>a sequence of samples from the distribution.</returns>
        public static IEnumerable<double> Samples(System.Random rnd, double mu, double sigma)
        {
            return Normal.Samples(rnd, mu, sigma).Select(Math.Exp);
        }
    }
}
