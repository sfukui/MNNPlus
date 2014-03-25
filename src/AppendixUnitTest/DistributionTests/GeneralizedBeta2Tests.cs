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
using System.Linq;
using MathNet.Numerics.Distributions;
using NUnit.Framework;

namespace AppendixUnitTest.DistributionTests
{
    /// <summary>
    /// Tests of generalized beta disturion of second kind.
    /// </summary>
    [TestFixture, Category("Distributions")]
    public class GeneralizedBeta2Tests
    {
        double[] m_GoodParameters, m_BadParameters, m_ParametersTooShort, m_ParametersTooLong;

        /// <summary>
        /// Set up parameters
        /// </summary>
        [SetUp]
        public void SetUp()
        {
            m_GoodParameters = new double[4] { 1.0, 700.0, 2.0, 3.0 };
            m_BadParameters = new double[4] { 0.0, -0.1, Double.NegativeInfinity, Double.NaN };
            m_ParametersTooShort = new double[3] { 1.0, 700.0, 2.0 };
            m_ParametersTooLong = new double[5] { 1.0, 700.0, 2.0, 3.0, 4.0 };
        }

        /// <summary>
        /// Can create GB2.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        [TestCase(0.01, 0.01, 0.01, 0.01)]
        [TestCase(1.0, 1.0, 1.0, 1.0)]
        [TestCase(10.0, 10.0, 10.0, 10.0)]
        [TestCase(1.0, 700.0, 2.0, 3.0)]
        [TestCase(0.01, 1.0, 10.0, 0.1)]
        [TestCase(Double.PositiveInfinity, Double.PositiveInfinity, Double.PositiveInfinity, Double.PositiveInfinity)]
        public void CanCreateGB2(double alpha, double beta, double p, double q)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(alpha, gb2.Alpha);
            Assert.AreEqual(beta, gb2.Beta);
            Assert.AreEqual(p, gb2.P);
            Assert.AreEqual(q, gb2.Q);
        }

        /// <summary>
        /// Can create GB2 with a double array "m_GoodParameters".
        /// </summary>
        [Test]
        public void CanCreateGB2WithDoubleArray()
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            for(int i = 0 ; i < m_GoodParameters.Length ; i++)
            {
                Assert.AreEqual(m_GoodParameters[i], gb2.Parameters[i]);
            }
        }

        /// <summary>
        /// GB2 creation fails with bad parameter value.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        [TestCase(0.0, 0.01, 0.01, 0.01)]
        [TestCase(0.01, 0.0, 0.01, 0.01)]
        [TestCase(0.01, 0.01, 0.0, 0.01)]
        [TestCase(0.01, 0.01, 0.01, 0.0)]
        [TestCase(-1.0, 1.0, 1.0, 1.0)]
        [TestCase(1.0, -1.0, 1.0, 1.0)]
        [TestCase(1.0, 1.0, -1.0, 1.0)]
        [TestCase(1.0, 1.0, 1.0, -1.0)]
        [TestCase(Double.NaN, 10.0, 10.0, 10.0)]
        [TestCase(10.0, Double.NaN, 10.0, 10.0)]
        [TestCase(10.0, 10.0, Double.NaN, 10.0)]
        [TestCase(10.0, 10.0, 10.0, Double.NaN)]
        [TestCase(Double.NegativeInfinity, 0.1, 1.0, 10.0)]
        [TestCase(0.01, Double.NegativeInfinity, 1.0, 10.0)]
        [TestCase(0.01, 0.1, Double.NegativeInfinity, 10.0)]
        [TestCase(0.01, 0.1, 1.0, Double.NegativeInfinity)]
        public void GB2CreateFailWithBadParameters(double alpha, double beta, double p, double q)
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new GeneralizedBeta2(alpha, beta, p, q));
        }

        /// <summary>
        /// GB2 creation fails with an array containing bad values.
        /// </summary>
        [Test]
        public void GB2CreateFailWithBadParametersArray()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new GeneralizedBeta2(m_BadParameters));
        }

        /// <summary>
        /// GB2 creation fails with short array of parameters.
        /// </summary>
        [Test]
        public void GB2CreateFailWithTooShortParametersArray()
        {
            Assert.Throws<ArgumentException>(() => new GeneralizedBeta2(m_ParametersTooShort));
        }

        /// <summary>
        /// GB2 creation fails with long array of parameters.
        /// </summary>
        [Test]
        public void GB2CreateFailWithTooLongParametersArray()
        {
            Assert.Throws<ArgumentException>(() => new GeneralizedBeta2(m_ParametersTooLong));
        }

        /// <summary>
        /// Validate ToString
        /// </summary>
        [Test]
        public void ValidateToString()
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.AreEqual("GB2(alpha = 1, beta = 700, p = 2, q = 3)", gb2.ToString());
        }

        /// <summary>
        /// Can set alpha.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        [TestCase(0.01)]
        [TestCase(0.1)]
        [TestCase(1.0)]
        [TestCase(10.0)]
        [TestCase(Double.PositiveInfinity)]
        public void CanSetAlpha(double alpha)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            gb2.Alpha = alpha;
            Assert.AreEqual(alpha, gb2.Alpha);
        }

        /// <summary>
        /// Set alpha fails with invalid values.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        [TestCase(0.0)]
        [TestCase(-1.0)]
        [TestCase(Double.NaN)]
        [TestCase(Double.NegativeInfinity)]
        public void SetAlphaFailWithInvalidValue(double alpha)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.Throws<ArgumentOutOfRangeException>(() => gb2.Alpha = alpha);
        }

        /// <summary>
        /// Can set beta.
        /// </summary>
        /// <param name="beta">Value of beta.</param>
        [TestCase(0.01)]
        [TestCase(0.1)]
        [TestCase(1.0)]
        [TestCase(10.0)]
        [TestCase(Double.PositiveInfinity)]
        public void CanSetBeta(double beta)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            gb2.Beta = beta;
            Assert.AreEqual(beta, gb2.Beta);
        }

        /// <summary>
        /// Set beta fails with invalid values.
        /// </summary>
        /// <param name="beta">Value of beta.</param>
        [TestCase(0.0)]
        [TestCase(-1.0)]
        [TestCase(Double.NaN)]
        [TestCase(Double.NegativeInfinity)]
        public void SetBetaFailWithInvalidValue(double beta)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.Throws<ArgumentOutOfRangeException>(() => gb2.Beta = beta);
        }

        /// <summary>
        /// Can set p.
        /// </summary>
        /// <param name="p">Value of p.</param>
        [TestCase(0.01)]
        [TestCase(0.1)]
        [TestCase(1.0)]
        [TestCase(10.0)]
        [TestCase(Double.PositiveInfinity)]
        public void CanSetP(double p)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            gb2.P = p;
            Assert.AreEqual(p, gb2.P);
        }

        /// <summary>
        /// Set p fails with invalid values.
        /// </summary>
        /// <param name="p">Value of p.</param>
        [TestCase(0.0)]
        [TestCase(-1.0)]
        [TestCase(Double.NaN)]
        [TestCase(Double.NegativeInfinity)]
        public void SetPFailWithInvalidValue(double p)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.Throws<ArgumentOutOfRangeException>(() => gb2.P = p);
        }

        /// <summary>
        /// Can set q.
        /// </summary>
        /// <param name="q">Value of q.</param>
        [TestCase(0.01)]
        [TestCase(0.1)]
        [TestCase(1.0)]
        [TestCase(10.0)]
        [TestCase(Double.PositiveInfinity)]
        public void CanSetQ(double q)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            gb2.Q = q;
            Assert.AreEqual(q, gb2.Q);
        }

        /// <summary>
        /// Set q fails with invalid values.
        /// </summary>
        /// <param name="q">Value of q.</param>
        [TestCase(0.0)]
        [TestCase(-1.0)]
        [TestCase(Double.NaN)]
        [TestCase(Double.NegativeInfinity)]
        public void SetQFailWithInvalidValue(double q)
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.Throws<ArgumentOutOfRangeException>(() => gb2.Q = q);
        }

        /// <summary>
        /// Validate mean.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        /// <param name="xmax">Upper limit in numerical integration.</param>
        /// <param name="mean">Expected mean.</param>
        [TestCase(1.0, 10.0, 4.0, 2.0, 200.0, 31.34105)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 5000.0, 646.5517)]
        public void ValidateMean(double alpha, double beta, double p, double q, double xmax, double mean)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(mean, gb2.Mean(xmax), 1e-3);
        }

        /// <summary>
        /// Validate variance.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        /// <param name="xmax">Upper limit in numerical integration.</param>
        /// <param name="variance">Expected variance.</param>
        [TestCase(1.0, 10.0, 4.0, 2.0, 200.0, 974.1118)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 5000.0, 449517.2)]
        public void ValidateVariance(double alpha, double beta, double p, double q, double xmax, double variance)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(variance, gb2.Variance(xmax), 1.0);
        }

        /// <summary>
        /// Validate standard deviation.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        /// <param name="xmax">Upper limit in numerical integration.</param>
        /// <param name="stddev">Expected standard deviation.</param>
        [TestCase(1.0, 10.0, 4.0, 2.0, 200.0, 31.21076)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 5000.0, 670.4605)]
        public void ValidateStdDev(double alpha, double beta, double p, double q, double xmax, double stddev)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(stddev, gb2.StdDev(xmax), 1e-3);
        }

        /// <summary>
        /// Validate statistics with parameter array.
        /// </summary>
        /// <param name="xmax">Upper limit in numerical integration.</param>
        /// <param name="mean">Expected mean.</param>
        /// <param name="variance">Expected variance.</param>
        /// <param name="stddev">Expected standard deviation.</param>
        [Test]
        public void ValidateStatisticsWithParameterArray()
        {
            double xmax = 5000.0, mean = 646.5517, variance = 449517.2, stddev = 670.4605; 
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.AreEqual(mean, gb2.Mean(xmax), 1e-3);
            Assert.AreEqual(variance, gb2.Variance(xmax), 1.0);
            Assert.AreEqual(stddev, gb2.StdDev(xmax), 1e-3);
        }

        /// <summary>
        /// Validate minimum.
        /// </summary>
        [Test]
        public void ValidateMinimum()
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.AreEqual(0.0, gb2.Minimum);
        }

        /// <summary>
        /// Validate maximum.
        /// </summary>
        [Test]
        public void ValideteMaximum()
        {
            var gb2 = new GeneralizedBeta2(m_GoodParameters);
            Assert.AreEqual(Double.PositiveInfinity, gb2.Maximum);
        }

        /// <summary>
        /// Validate density.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        /// <param name="x">Input x value.</param>
        /// <param name="expected">Expected value of density.</param>
        [TestCase(1.0, 10.0, 4.0, 2.0, 0.0, 0.0)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 1.0, 0.001128948)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 10.0, 0.03125)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 30.0, 0.01318359)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 80.0, 0.001926837)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 0.0, 0.0)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 10.0, 0.0002281306)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 250.0, 0.001329834)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 1000.0, 0.0002898883)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 1500.0, 0.000119799)]
        public void ValidateDensity(double alpha, double beta, double p, double q, double x, double expected)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(expected, gb2.Density(x), 1e-4);
        }

        /// <summary>
        /// Validatge log of density.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        /// <param name="x">Input x value.</param>
        /// <param name="expected">Expected value of log of density.</param>
        [TestCase(1.0, 10.0, 4.0, 2.0, 0.0, Double.NegativeInfinity)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 1.0, -6.786469)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 10.0, -3.465736)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 30.0, -4.328782)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 80.0, -6.251876)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 0.0, Double.NegativeInfinity)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 10.0, -8.385592)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 250.0, -6.622701)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 1000.0, -8.146015)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 1500.0, -9.029695)]
        public void ValidateDensityLn(double alpha, double beta, double p, double q, double x, double expected)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(expected, gb2.DensityLn(x), 1e-5);
        }

        /// <summary>
        /// Validate cumulative density.
        /// </summary>
        /// <param name="alpha">Value of alpha.</param>
        /// <param name="beta">Value of beta.</param>
        /// <param name="p">Value of p.</param>
        /// <param name="q">Value of q.</param>
        /// <param name="x">Input x value.</param>
        /// <param name="expected">Expected value of cumulative density.</param>
        [TestCase(1.0, 10.0, 4.0, 2.0, 0.0, 0.0)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 1.0, 0.0003166699)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 10.0, 0.1875)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 30.0, 0.6328125)]
        [TestCase(1.0, 10.0, 4.0, 2.0, 80.0, 0.9017596)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 0.0, 0.0)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 10.0, 0.001168006)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 250.0, 0.2841062)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 1000.0, 0.8069827)]
        [TestCase(1.0, 700.0, 2.0, 3.0, 1500.0, 0.9018979)]
        public void ValidateCumulativeDensity(double alpha, double beta, double p, double q, double x, double expected)
        {
            var gb2 = new GeneralizedBeta2(alpha, beta, p, q);
            Assert.AreEqual(expected, gb2.CumulativeDistribution(x), 1e-3);
        }
    }
}
