// <copyright file="ConstrainedHMCTests.cs" company="Fukui Shogo">
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
// Copyright (c) 2009-2016 Math.NET
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
using System.Text;
using System.Threading.Tasks;
using System.Globalization;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Appendix.Statistics.Mcmc;
using NUnit.Framework;

namespace AppendixUnitTest.MCMCTests
{
    /// <summary>
    /// Tests for the ConstrainedHMC class.
    /// </summary>
    [TestFixture, Category("Statistics")]
    public class ConstrainedHMCTest
    {
        readonly Normal _normal = new Normal(0, 1);

        /// <summary>
        /// Testing the constructor to make sure that RandomSource is
        /// assigned.
        /// </summary>
        [Test]
        public void RandomSourceTest()
        {
            var constrainedHMC = new ConstrainedHMC(new double[] { 0 }, (double[] x) => _normal.DensityLn(x[0]), 10, 0.1);
            Assert.IsNotNull(constrainedHMC.RandomSource);

            constrainedHMC.RandomSource = new System.Random(0);
            Assert.IsNotNull(constrainedHMC.RandomSource);
        }

        /// <summary>
        /// Test the range of FrogLeapSteps. Sets FrogLeapSteps
        /// to negative or zero throws <c>ArgumentOutOfRangeException</c>.
        /// </summary>
        [Test]
        public void FrogLeapTest()
        {
            Assert.Throws<ArgumentOutOfRangeException>(()
                => new ConstrainedHMC(new double[] { 0 }, (double[] x) => _normal.DensityLn(x[0]), 0, 0.1));
        }

        /// <summary>
        /// Test the range of StepSize. Sets StepSize
        /// to negative or zero throws <c>ArgumentOutOfRangeException</c>.
        /// </summary>
        [Test]
        public void StepSizeTest()
        {
            Assert.Throws<ArgumentOutOfRangeException>(()
                => new ConstrainedHMC(new double[] { 0 }, (double[] x) => _normal.DensityLn(x[0]), 1, 0));
        }

        /// <summary>
        /// Test the range of BurnInterval. Sets BurnInterval
        /// to negative throws <c>ArgumentOutOfRangeException</c>.
        /// </summary>
        [Test]
        public void BurnIntervalTest()
        {
            Assert.Throws<ArgumentOutOfRangeException>(()
                => new ConstrainedHMC(new double[] { 0 }, (double[] x) => _normal.DensityLn(x[0]), 10, 0.1, -1));
        }


        /// <summary>
        /// Test the range of MomentumStdDev. Sets MomentumStdDev
        /// to negative throws <c>ArgumentNullException</c>.
        /// </summary>
        [Test]
        public void MomentumStdDevNegativeTest()
        {
            Assert.Throws<ArgumentOutOfRangeException>(()
                => new ConstrainedHMC(new double[] { 0 }, (double[] x) => _normal.DensityLn(x[0]), 10, 0.1, 0, new double[] { -1 }));
        }

        /// <summary>
        /// Test the length of MomentumStdDev. Sets MomentumStdDev
        /// to a length different from the length of samples throws <c>ArgumentOutOfRangeException</c>.
        /// </summary>
        [Test]
        public void MomentumStdDevLengthTest()
        {
            Assert.Throws<ArgumentOutOfRangeException>(()
                => new ConstrainedHMC(new double[] { 0 }, (double[] x) => _normal.DensityLn(x[0]), 10, 0.1, 0, new double[2]));
        }

        [TestCase(new[] { 1.0, 2.0 }, new[] { 0.1 }, new[] { 10.0, 20.0 })]
        [TestCase(new[] { 1.0, 2.0 }, new[] { 0.1, 0.2 }, new[] { 10.0 })]
        [TestCase(new[] { 1.0 }, new[] { 0.1, 0.2 }, new[] { 10.0, 20.0 })]
        public void SampleValuesInfimumAndSupremumLengthTest(double[] x0, double[] xInfs, double[] xSups)
        {
            var pdfLn = new MathNet.Numerics.Statistics.Mcmc.DensityLn<double[]>(x => LogDen(x, new[] { 1.0, 1.0 }, new[] { 0.0, 0.0 }, 0.1));
            Assert.Throws<ArgumentException>(()
                => new ConstrainedHMC(x0, pdfLn, 10, 0.1, 1000, xInfs, xSups));
        }

        [TestCase(new[] { 100.0, 0.2 }, new[] { 10.0, 20.0 })]
        [TestCase(new[] { 0.1, 30.0 }, new[] { 10.0, 20.0 })]
        [TestCase(new[] { 10.0, 20.0 }, new[] { 10.0, 20.0 })]
        public void InfimumsAreGreaterThanSupremumTest(double[] xInfs, double[] xSups)
        {
            var pdfLn = new MathNet.Numerics.Statistics.Mcmc.DensityLn<double[]>(x => LogDen(x, new[] { 1.0, 1.0 }, new[] { 0.0, 0.0 }, 0.1));
            Assert.Throws<ArgumentException>(()
                => new ConstrainedHMC(new[] { 1.0, 2.0 }, pdfLn, 10, 0.1, 1000, xInfs, xSups));
        }

        [TestCase(new[] { Double.NaN, 0.2 }, new[] { 10.0, 20.0 })]
        [TestCase(new[] { 0.1, 0.2 }, new[] { 10.0, Double.NaN })]
        [TestCase(new[] { 10.0, Double.NaN }, new[] { Double.NaN, 20.0 })]
        public void InfimumsAndSupremumNaNTest(double[] xInfs, double[] xSups)
        {
            var pdfLn = new MathNet.Numerics.Statistics.Mcmc.DensityLn<double[]>(x => LogDen(x, new[] { 1.0, 1.0 }, new[] { 0.0, 0.0 }, 0.1));
            Assert.Throws<ArgumentException>(()
                => new ConstrainedHMC(new[] { 1.0, 2.0 }, pdfLn, 10, 0.1, 1000, xInfs, xSups));
        }

        /// <summary>
        /// Test the sampler using  a bivariate normal distribution with randomly selected mean and standard deviation.
        /// It is a statistical test and may not pass all the time. Note that Sdv and rho have to be between 0 and 1.
        /// </summary>
        [TestCase(new[] { 0.8, 0.2 }, new[] { 3.2, -4.6 }, 0.77, 1000)]
        [TestCase(new[] { 0.5, 0.1 }, new[] { -2.2, -1.3 }, 0.29, 1000)]
        [TestCase(new[] { 0.45, 0.78 }, new[] { 1.34, -3.3 }, 0.58, 1000)]
        public void NotConstrainedSampleTest(double[] sdv, double[] mean, double rho, int seed)
        {
            var pdfLn = new MathNet.Numerics.Statistics.Mcmc.DensityLn<double[]>(x => LogDen(x, sdv, mean, rho));
            var hybrid = new ConstrainedHMC(new double[] { 0, 0 }, pdfLn, 10, 0.1, 1000, new double[] { 1, 1 }, new System.Random(seed))
            {
                BurnInterval = 0
            };
            
            const int sampleSize = 10000;
            double[][] sample = hybrid.Sample(sampleSize);
            double[][] newSamples = ArrangeSamples(sampleSize, sample);

            var convergence = new double[2];
            var sampleMean = new double[2];
            var sampleSdv = new double[2];

            for (int i = 0; i < 2; i++)
            {
                // ReSharper disable once AccessToModifiedClosure
                convergence[i] = 1 / Math.Sqrt(MathNet.Numerics.Statistics.Mcmc.MCMCDiagnostics.EffectiveSize(sample, x => x[i]));
                var stats = new MathNet.Numerics.Statistics.DescriptiveStatistics(newSamples[i]);
                sampleMean[i] = stats.Mean;
                sampleSdv[i] = stats.StandardDeviation;

            }

            double sampleRho = MathNet.Numerics.Statistics.Correlation.Pearson(newSamples[0], newSamples[1]);

            for (int i = 0; i < 2; i++)
            {
                string index = i.ToString(CultureInfo.InvariantCulture);
                Assert.AreEqual(mean[i], sampleMean[i], 10 * convergence[i], index + "Mean");
                Assert.AreEqual(sampleSdv[i] * sampleSdv[i], sdv[i] * sdv[i], 10 * convergence[i], index + "Standard Deviation");
            }

            double convergenceRho = 1 / Math.Sqrt(MathNet.Numerics.Statistics.Mcmc.MCMCDiagnostics.EffectiveSize(sample, x => (x[0] - sampleMean[0]) * (x[1] - sampleMean[1])));
            Assert.AreEqual(sampleRho * sampleSdv[0] * sampleSdv[1], rho * sdv[0] * sdv[1], 10 * convergenceRho, "Rho");
        }

        /// <summary>
        /// Test the sampler using  a bivariate normal distribution whose domain are truncated.
        /// It is a statistical test and may not pass all the time. Note that Sdv and rho have to be between 0 and 1.
        /// </summary>
        [Test]
        public void SampleFromDomainTruncatedDensityTest()
        {
            var sdv = new[] { 1.0, 1.0 };
            var mean = new[] { 0.0, 0.0 };
            var rho = 0.0;
            int seed = 1000;
            var xInfs = new[] { -2.0, -2.0 };
            var xSups = new[] { 2.0, 2.0 };

            var pdfLn = new MathNet.Numerics.Statistics.Mcmc.DensityLn<double[]>(x => LogDen(x, sdv, mean, rho));
            var hybrid = new ConstrainedHMC(new double[] { 0, 0 }, pdfLn, 10, 0.1, 1000, new double[] { 1, 1 }, new System.Random(seed), xInfs, xSups)
            {
                BurnInterval = 0
            };

            const int sampleSize = 10000;
            double[][] sample = hybrid.Sample(sampleSize);
            double[][] newSamples = ArrangeSamples(sampleSize, sample);

            var convergence = new double[2];
            var sampleMean = new double[2];
            var sampleSdv = new double[2];

            for (int i = 0; i < 2; i++)
            {
                // ReSharper disable once AccessToModifiedClosure
                convergence[i] = 1 / Math.Sqrt(MathNet.Numerics.Statistics.Mcmc.MCMCDiagnostics.EffectiveSize(sample, x => x[i]));
                var stats = new MathNet.Numerics.Statistics.DescriptiveStatistics(newSamples[i]);
                sampleMean[i] = stats.Mean;
                sampleSdv[i] = stats.StandardDeviation;

            }

            double sampleRho = MathNet.Numerics.Statistics.Correlation.Pearson(newSamples[0], newSamples[1]);

            for (int i = 0; i < 2; i++)
            {
                string index = i.ToString(CultureInfo.InvariantCulture);
                Assert.AreEqual(mean[i], sampleMean[i], 10 * convergence[i], index + "Mean");
                Assert.AreEqual(newSamples[i].Min(), xInfs[i], 0.1, index + "Min");
                Assert.AreEqual(newSamples[i].Max(), xSups[i], 0.1, index + "Max");
            }

            double convergenceRho = 1 / Math.Sqrt(MathNet.Numerics.Statistics.Mcmc.MCMCDiagnostics.EffectiveSize(sample, x => (x[0] - sampleMean[0]) * (x[1] - sampleMean[1])));
            Assert.AreEqual(sampleRho * sampleSdv[0] * sampleSdv[1], rho * sdv[0] * sdv[1], 10 * convergenceRho, "Rho");
        }

        /// <summary>
        /// Test the sampler using  a bivariate normal distribution whose domain are truncated.
        /// It is a statistical test and may not pass all the time. Note that Sdv and rho have to be between 0 and 1.
        /// </summary>
        [Test]
        public void SampleFromDensityWithInvalidResultTest()
        {
            var sdv = new[] { 1.0, 1.0 };
            var mean = new[] { 0.0, 0.0 };
            var rho = 0.0;
            int seed = 1000;
            var xInfs = new[] { -2.0, -2.0 };
            var xSups = new[] { 2.0, 2.0 };

            var pdfLn = new MathNet.Numerics.Statistics.Mcmc.DensityLn<double[]>(x => LogDenInvalidResult(x, sdv, mean, rho, xInfs, xSups));
            var hybrid = new ConstrainedHMC(new double[] { 0, 0 }, pdfLn, 10, 0.1, 1000, new double[] { 1, 1 }, new System.Random(seed))
            {
                BurnInterval = 0
            };

            const int sampleSize = 10000;
            double[][] sample = hybrid.Sample(sampleSize);
            double[][] newSamples = ArrangeSamples(sampleSize, sample);

            var convergence = new double[2];
            var sampleMean = new double[2];
            var sampleSdv = new double[2];

            for (int i = 0; i < 2; i++)
            {
                // ReSharper disable once AccessToModifiedClosure
                convergence[i] = 1 / Math.Sqrt(MathNet.Numerics.Statistics.Mcmc.MCMCDiagnostics.EffectiveSize(sample, x => x[i]));
                var stats = new MathNet.Numerics.Statistics.DescriptiveStatistics(newSamples[i]);
                sampleMean[i] = stats.Mean;
                sampleSdv[i] = stats.StandardDeviation;

            }

            double sampleRho = MathNet.Numerics.Statistics.Correlation.Pearson(newSamples[0], newSamples[1]);

            for (int i = 0; i < 2; i++)
            {
                string index = i.ToString(CultureInfo.InvariantCulture);
                Assert.AreEqual(mean[i], sampleMean[i], 10 * convergence[i], index + "Mean");
                Assert.AreEqual(newSamples[i].Min(), xInfs[i], 0.1, index + "Min");
                Assert.AreEqual(newSamples[i].Max(), xSups[i], 0.1, index + "Max");
            }

            double convergenceRho = 1 / Math.Sqrt(MathNet.Numerics.Statistics.Mcmc.MCMCDiagnostics.EffectiveSize(sample, x => (x[0] - sampleMean[0]) * (x[1] - sampleMean[1])));
            Assert.AreEqual(sampleRho * sampleSdv[0] * sampleSdv[1], rho * sdv[0] * sdv[1], 10 * convergenceRho, "Rho");
        }

        /// <summary>
        /// The log density of the bivariate normal distribution used in SampleTest.
        /// </summary>
        /// <param name="x">Location to evaluate the density.</param>
        /// <param name="sdv">Standard deviation.</param>
        /// <param name="mean">Mean.</param>
        /// <param name="rho">Correlation of the two variables.</param>
        /// <returns>Value of the log density.</returns>
        double LogDen(double[] x, double[] sdv, double[] mean, double rho)
        {
            if (x.Length != 2 || sdv.Length != 2 || mean.Length != 2)
                throw new ArgumentException("LogDen must take a 2 dimensional array");

            double xDiv = x[0] - mean[0];
            double yDiv = x[1] - mean[1];
            double xVar = sdv[0] * sdv[0];
            double yVar = sdv[1] * sdv[1];

            return (-(0.5 / (1 - rho * rho)) * ((xDiv * xDiv) / xVar + (yDiv * yDiv) / yVar - 2 * rho * xDiv * yDiv / (sdv[0] * sdv[1])));
        }

        /// <summary>
        /// The log density of the bivariate normal distribution used in SampleTest.
        /// </summary>
        /// <param name="x">Location to evaluate the density.</param>
        /// <param name="sdv">Standard deviation.</param>
        /// <param name="mean">Mean.</param>
        /// <param name="rho">Correlation of the two variables.</param>
        /// <param name="xInfs">The infimums of the values of sample.</param>
        /// <param name="xSups">The supremums of the values of sample.</param>
        /// <returns>If x[0] is out of the bound, then return Double.NaN. Else if x[1] is equal to or lower than the infimum, then return negative infinity. Else if x[1] is equal to or greater than the supremum, then return positive infinity. Otherwise, return the value of the log density of bivariate normal distribution.</returns>
        double LogDenInvalidResult(double[] x, double[] sdv, double[] mean, double rho, double[] xInfs, double[] xSups)
        {
            if (x.Length != 2 || sdv.Length != 2 || mean.Length != 2)
                throw new ArgumentException("LogDen must take a 2 dimensional array");

            double xDiv = x[0] - mean[0];
            double yDiv = x[1] - mean[1];
            double xVar = sdv[0] * sdv[0];
            double yVar = sdv[1] * sdv[1];

            if (x[0] <= xInfs[0] || x[0] >= xSups[0])
                return Double.NaN;
            else if (x[1] <= xInfs[1])
                return Double.NegativeInfinity;
            else if (x[1] >= xSups[1])
                return Double.PositiveInfinity;
            else
                return (-(0.5 / (1 - rho * rho)) * ((xDiv * xDiv) / xVar + (yDiv * yDiv) / yVar - 2 * rho * xDiv * yDiv / (sdv[0] * sdv[1])));
        }

        /// <summary>
        /// Method to rearrange the array of samples in to separated arrays.
        /// </summary>
        /// <param name="sampleSize">Size of the sample.</param>
        /// <param name="sample">Sample from the HybridMC.</param>
        /// <returns>An array whose first entry is the samples in the first variable and
        /// second entry is the samples in the second variable.</returns>
        double[][] ArrangeSamples(int sampleSize, double[][] sample)
        {
            var xSample = new double[sampleSize];
            var ySample = new double[sampleSize];

            for (int i = 0; i < sampleSize; i++)
            {
                xSample[i] = sample[i][0];
                ySample[i] = sample[i][1];
            }

            return new[] { xSample, ySample };
        }
    }
}
