﻿// <copyright file="Categorical.cs" company="Math.NET">
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
using MathNet.Numerics.Statistics;

namespace MathNet.Numerics.Distributions
{
    /// <summary>
    /// Discrete Univariate Categorical distribution.
    /// For details about this distribution, see 
    /// <a href="http://en.wikipedia.org/wiki/Categorical_distribution">Wikipedia - Categorical distribution</a>. This
    /// distribution is sometimes called the Discrete distribution.
    /// </summary>
    /// <remarks><para>The distribution is parameterized by a vector of ratios: in other words, the parameter
    /// does not have to be normalized and sum to 1. The reason is that some vectors can't be exactly normalized
    /// to sum to 1 in floating point representation.</para>
    /// <para>The distribution will use the <see cref="System.Random"/> by default. 
    /// Users can set the random number generator by using the <see cref="RandomSource"/> property.</para>
    /// <para>The statistics classes will check all the incoming parameters whether they are in the allowed
    /// range. This might involve heavy computation. Optionally, by setting Control.CheckDistributionParameters
    /// to <c>false</c>, all parameter checks can be turned off.</para></remarks>
    public class Categorical : IDiscreteDistribution
    {
        System.Random _random;

        double[] _pmfNormalized;
        double[] _cdfUnnormalized;

        /// <summary>
        /// Initializes a new instance of the Categorical class.
        /// </summary>
        /// <param name="probabilityMass">An array of nonnegative ratios: this array does not need to be normalized 
        /// as this is often impossible using floating point arithmetic.</param>
        /// <exception cref="ArgumentException">If any of the probabilities are negative or do not sum to one.</exception>
        public Categorical(double[] probabilityMass)
        {
            _random = new System.Random(Random.RandomSeed.Guid());
            SetParameters(probabilityMass);
        }

        /// <summary>
        /// Initializes a new instance of the Categorical class.
        /// </summary>
        /// <param name="probabilityMass">An array of nonnegative ratios: this array does not need to be normalized 
        /// as this is often impossible using floating point arithmetic.</param>
        /// <param name="randomSource">The random number generator which is used to draw random samples.</param>
        /// <exception cref="ArgumentException">If any of the probabilities are negative or do not sum to one.</exception>
        public Categorical(double[] probabilityMass, System.Random randomSource)
        {
            _random = randomSource ?? new System.Random(Random.RandomSeed.Guid());
            SetParameters(probabilityMass);
        }

        /// <summary>
        /// Initializes a new instance of the Categorical class from a <paramref name="histogram"/>. The distribution 
        /// will not be automatically updated when the histogram changes. The categorical distribution will have
        /// one value for each bucket and a probability for that value proportional to the bucket count.
        /// </summary>
        /// <param name="histogram">The histogram from which to create the categorical variable.</param>
        public Categorical(Histogram histogram)
        {
            if (histogram == null)
            {
                throw new ArgumentNullException("histogram");
            }

            // The probability distribution vector.
            var p = new double[histogram.BucketCount];

            // Fill in the distribution vector.
            for (var i = 0; i < histogram.BucketCount; i++)
            {
                p[i] = histogram[i].Count;
            }

            _random = new System.Random(Random.RandomSeed.Guid());
            SetParameters(p);
        }

        /// <summary>
        /// A string representation of the distribution.
        /// </summary>
        /// <returns>a string representation of the distribution.</returns>
        public override string ToString()
        {
            return "Categorical(Dimension = " + _pmfNormalized.Length + ")";
        }

        /// <summary>
        /// Checks whether the parameters of the distribution are valid. 
        /// </summary>
        /// <param name="p">An array of nonnegative ratios: this array does not need to be normalized as this is often impossible using floating point arithmetic.</param>
        /// <returns>If any of the probabilities are negative returns <c>false</c>, or if the sum of parameters is 0.0; otherwise <c>true</c></returns>
        static bool IsValidProbabilityMass(double[] p)
        {
            var sum = 0.0;
            for (int i = 0; i < p.Length; i++)
            {
                double t = p[i];
                if (t < 0.0 || Double.IsNaN(t))
                {
                    return false;
                }

                sum += t;
            }

            return sum > 0.0;
        }

        /// <summary>
        /// Checks whether the parameters of the distribution are valid. 
        /// </summary>
        /// <param name="cdf">An array of nonnegative ratios: this array does not need to be normalized as this is often impossible using floating point arithmetic.</param>
        /// <returns>If any of the probabilities are negative returns <c>false</c>, or if the sum of parameters is 0.0; otherwise <c>true</c></returns>
        static bool IsValidCumulativeDistribution(double[] cdf)
        {
            var last = 0.0;
            for (int i = 0; i < cdf.Length; i++)
            {
                double t = cdf[i];
                if (t < 0.0 || Double.IsNaN(t) || t < last)
                {
                    return false;
                }

                last = t;
            }

            return last > 0.0;
        }

        /// <summary>
        /// Sets the parameters of the distribution after checking their validity.
        /// </summary>
        /// <param name="p">An array of nonnegative ratios: this array does not need to be normalized 
        /// as this is often impossible using floating point arithmetic.</param>
        /// <exception cref="ArgumentOutOfRangeException">When the parameters are out of range.</exception>
        void SetParameters(double[] p)
        {
            if (Control.CheckDistributionParameters && !IsValidProbabilityMass(p))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            // Extract unnormalized cumulative distribution
            _cdfUnnormalized = new double[p.Length];
            _cdfUnnormalized[0] = p[0];
            for (int i = 1; i < p.Length; i++)
            {
                _cdfUnnormalized[i] = _cdfUnnormalized[i - 1] + p[i];
            }

            // Extract normalized probability mass
            var sum = _cdfUnnormalized[_cdfUnnormalized.Length - 1];
            _pmfNormalized = new double[p.Length];
            for (int i = 0; i < p.Length; i++)
            {
                _pmfNormalized[i] = p[i]/sum;
            }
        }

        /// <summary>
        /// Gets or sets the probability mass vector (non-negative ratios) of the multinomial.
        /// </summary>
        /// <remarks>Sometimes the normalized probability vector cannot be represented exactly in a floating point representation.</remarks>
        public double[] P
        {
            get { return (double[])_pmfNormalized.Clone(); }
            set { SetParameters(value); }
        }

        /// <summary>
        /// Gets or sets the random number generator which is used to draw random samples.
        /// </summary>
        public System.Random RandomSource
        {
            get { return _random; }
            set { _random = value ?? new System.Random(Random.RandomSeed.Guid()); }
        }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        public double Mean
        {
            get { return _pmfNormalized.Mean(); }
        }

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
        public double StdDev
        {
            get { return _pmfNormalized.StandardDeviation(); }
        }

        /// <summary>
        /// Gets the variance of the distribution.
        /// </summary>
        public double Variance
        {
            get { return _pmfNormalized.Variance(); }
        }

        /// <summary>
        /// Gets the entropy of the distribution.
        /// </summary>
        public double Entropy
        {
            get { return _pmfNormalized.Sum(p => p*Math.Log(p)); }
        }

        /// <summary>
        /// Gets the skewness of the distribution.
        /// </summary>
        /// <remarks>Throws a <see cref="NotSupportedException"/>.</remarks>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }

        /// <summary>
        /// Gets the smallest element in the domain of the distributions which can be represented by an integer.
        /// </summary>
        public int Minimum
        {
            get { return 0; }
        }

        /// <summary>
        /// Gets the largest element in the domain of the distributions which can be represented by an integer.
        /// </summary>
        public int Maximum
        {
            get { return _pmfNormalized.Length - 1; }
        }

        /// <summary>
        /// Gets he mode of the distribution.
        /// </summary>
        /// <remarks>Throws a <see cref="NotSupportedException"/>.</remarks>
        public int Mode
        {
            get { throw new NotSupportedException(); }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        public int Median
        {
            get { return (int) _pmfNormalized.Median(); }
        }

        /// <summary>
        /// Computes the probability mass (PMF) at k, i.e. P(X = k).
        /// </summary>
        /// <param name="k">The location in the domain where we want to evaluate the probability mass function.</param>
        /// <returns>the probability mass at location <paramref name="k"/>.</returns>
        public double Probability(int k)
        {
            if (k < 0)
            {
                return 0.0;
            }

            if (k >= _pmfNormalized.Length)
            {
                return 0.0;
            }

            return _pmfNormalized[k];
        }

        /// <summary>
        /// Computes the log probability mass (lnPMF) at k, i.e. ln(P(X = k)).
        /// </summary>
        /// <param name="k">The location in the domain where we want to evaluate the log probability mass function.</param>
        /// <returns>the log probability mass at location <paramref name="k"/>.</returns>
        public double ProbabilityLn(int k)
        {
            if (k < 0)
            {
                return 0.0;
            }

            if (k >= _pmfNormalized.Length)
            {
                return 0.0;
            }

            return Math.Log(_pmfNormalized[k]);
        }

        /// <summary>
        /// Computes the cumulative distribution (CDF) of the distribution at x, i.e. P(X ≤ x).
        /// </summary>
        /// <param name="x">The location at which to compute the cumulative distribution function.</param>
        /// <returns>the cumulative distribution at location <paramref name="x"/>.</returns>
        public double CumulativeDistribution(double x)
        {
            if (x < 0.0)
            {
                return 0.0;
            }

            if (x >= _cdfUnnormalized.Length)
            {
                return 1.0;
            }

            return _cdfUnnormalized[(int) Math.Floor(x)]/_cdfUnnormalized[_cdfUnnormalized.Length - 1];
        }

        /// <summary>
        /// Computes the inverse of the cumulative distribution function (InvCDF) for the distribution
        /// at the given probability.
        /// </summary>
        /// <param name="probability">A real number between 0 and 1.</param>
        /// <returns>An integer between 0 and the size of the categorical (exclusive), that corresponds to the inverse CDF for the given probability.</returns>
        public int InverseCumulativeDistribution(double probability)
        {
            if (probability < 0.0 || probability > 1.0 || Double.IsNaN(probability))
            {
                throw new ArgumentOutOfRangeException("probability");
            }

            var denormalizedProbability = probability * _cdfUnnormalized[_cdfUnnormalized.Length - 1];
            int idx = Array.BinarySearch(_cdfUnnormalized, denormalizedProbability);
            if (idx < 0)
            {
                idx = ~idx;
            }

            return idx;
        }

        /// <summary>
        /// Computes the inverse of the cumulative distribution function (InvCDF) for the distribution
        /// at the given probability.
        /// </summary>
        /// <param name="cdfUnnormalized">An array corresponding to a CDF for a categorical distribution. Not assumed to be normalized.</param>
        /// <param name="probability">A real number between 0 and 1.</param>
        /// <returns>An integer between 0 and the size of the categorical (exclusive), that corresponds to the inverse CDF for the given probability.</returns>
        public static int InverseCumulativeDistribution(double[] cdfUnnormalized, double probability)
        {
            if (Control.CheckDistributionParameters && !IsValidCumulativeDistribution(cdfUnnormalized))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            if (probability < 0.0 || probability > 1.0 || Double.IsNaN(probability))
            {
                throw new ArgumentOutOfRangeException("probability");
            }

            var denormalizedProbability = probability * cdfUnnormalized[cdfUnnormalized.Length - 1];
            int idx = Array.BinarySearch(cdfUnnormalized, denormalizedProbability);
            if (idx < 0)
            {
                idx = ~idx;
            }

            return idx;
        }

        /// <summary>
        /// Computes the cumulative distribution function. This method performs no parameter checking.
        /// If the probability mass was normalized, the resulting cumulative distribution is normalized as well (up to numerical errors).
        /// </summary>
        /// <param name="pmfUnnormalized">An array of nonnegative ratios: this array does not need to be normalized 
        /// as this is often impossible using floating point arithmetic.</param>
        /// <returns>An array representing the unnormalized cumulative distribution function.</returns>
        internal static double[] ProbabilityMassToCumulativeDistribution(double[] pmfUnnormalized)
        {
            var cdfUnnormalized = new double[pmfUnnormalized.Length];
            cdfUnnormalized[0] = pmfUnnormalized[0];
            for (int i = 1; i < pmfUnnormalized.Length; i++)
            {
                cdfUnnormalized[i] = cdfUnnormalized[i - 1] + pmfUnnormalized[i];
            }

            return cdfUnnormalized;
        }

        /// <summary>
        /// Returns one trials from the categorical distribution.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="cdfUnnormalized">The (unnormalized) cumulative distribution of the probability distribution.</param>
        /// <returns>One sample from the categorical distribution implied by <paramref name="cdfUnnormalized"/>.</returns>
        internal static int SampleUnchecked(System.Random rnd, double[] cdfUnnormalized)
        {
            // TODO : use binary search to speed up this procedure.
            var u = rnd.NextDouble()*cdfUnnormalized[cdfUnnormalized.Length - 1];

            var idx = 0;
            while (u > cdfUnnormalized[idx])
            {
                idx++;
            }

            return idx;
        }

        /// <summary>
        /// Samples a Binomially distributed random variable.
        /// </summary>
        /// <returns>The number of successful trials.</returns>
        public int Sample()
        {
            return SampleUnchecked(_random, _cdfUnnormalized);
        }

        /// <summary>
        /// Samples an array of Bernoulli distributed random variables.
        /// </summary>
        /// <returns>a sequence of successful trial counts.</returns>
        public IEnumerable<int> Samples()
        {
            while (true)
            {
                yield return SampleUnchecked(_random, _cdfUnnormalized);
            }
        }

        /// <summary>
        /// Samples one categorical distributed random variable; also known as the Discrete distribution.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="cdfUnnormalized">An array of the cumulative distribution. Not assumed to be normalized.</param>
        /// <returns>One random integer between 0 and the size of the categorical (exclusive).</returns>
        public static int SampleWithCumulativeDistribution(System.Random rnd, double[] cdfUnnormalized)
        {
            if (Control.CheckDistributionParameters && !IsValidCumulativeDistribution(cdfUnnormalized))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            return SampleUnchecked(rnd, cdfUnnormalized);
        }

        /// <summary>
        /// Samples one categorical distributed random variable; also known as the Discrete distribution.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="pmfUnnormalized">An array of nonnegative ratios. Not assumed to be normalized.</param>
        /// <returns>One random integer between 0 and the size of the categorical (exclusive).</returns>
        public static int SampleWithProbabilityMass(System.Random rnd, double[] pmfUnnormalized)
        {
            if (Control.CheckDistributionParameters && !IsValidProbabilityMass(pmfUnnormalized))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            var cdf = ProbabilityMassToCumulativeDistribution(pmfUnnormalized);
            return SampleUnchecked(rnd, cdf);
        }

        /// <summary>
        /// Samples a categorically distributed random variable.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="cdfUnnormalized">An array of the cumulative distribution. Not assumed to be normalized.</param>
        /// <returns>random integers between 0 and the size of the categorical (exclusive).</returns>
        public static IEnumerable<int> SamplesWithCumulativeDistribution(System.Random rnd, double[] cdfUnnormalized)
        {
            if (Control.CheckDistributionParameters && !IsValidCumulativeDistribution(cdfUnnormalized))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            while (true)
            {
                yield return SampleUnchecked(rnd, cdfUnnormalized);
            }
        }

        /// <summary>
        /// Samples a categorically distributed random variable.
        /// </summary>
        /// <param name="rnd">The random number generator to use.</param>
        /// <param name="pmfUnnormalized">An array of nonnegative ratios. Not assumed to be normalized.</param>
        /// <returns>random integers between 0 and the size of the categorical (exclusive).</returns>
        public static IEnumerable<int> SamplesWithProbabilityMass(System.Random rnd, double[] pmfUnnormalized)
        {
            if (Control.CheckDistributionParameters && !IsValidProbabilityMass(pmfUnnormalized))
            {
                throw new ArgumentOutOfRangeException(Resources.InvalidDistributionParameters);
            }

            var cdf = ProbabilityMassToCumulativeDistribution(pmfUnnormalized);
            while (true)
            {
                yield return SampleUnchecked(rnd, cdf);
            }
        }
    }
}
