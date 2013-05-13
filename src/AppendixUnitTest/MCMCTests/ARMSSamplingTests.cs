using System;
using System.Linq;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.Statistics.Mcmc;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.MCMCTests
{
    [TestFixture]
    public class ARMSSamplingTests
    {
        /// <summary>
        /// Rate of accept range for parameter comparison.
        /// </summary>
        private double acceptRangeRate = 0.01;

        /// <summary>
        /// Sample sizes in ARMS.
        /// </summary>
        private int Iteration_NormalTest = 50000;
        private int Iteration_GB2Test = 50000;

        /// <summary>
        /// Log of normal density function.
        /// </summary>
        /// <param name="parameters">Vector value of parameters.</param>
        /// <returns>Natural log of density.</returns>
        private double lnPDF_Normal(double x, Vector<double> parameters)
        {
            return (-0.5) * (System.Math.Log(2.0 * System.Math.PI) + (System.Math.Log(parameters[1]))) -
                    (System.Math.Pow((x - parameters[0]), 2.0) / (2.0 * parameters[1]));
        }

        /// <summary>
        /// Log of GB2 density function.
        /// </summary>
        /// <param name="x">Value of variable.</param>
        /// <param name="parameters">Vector value of parameters.</param>
        /// <returns>Natural log of density.</returns>
        private double lnPDF_GB2(double x, Vector<double> parameters)
        {
            double lnumer = Math.Log(parameters[0]) + (parameters[0] * parameters[2] - 1.0) * Math.Log(x);
            double ldenom = (parameters[0] * parameters[2]) * Math.Log(parameters[1]) +
                SpecialFunctions.BetaLn(parameters[2], parameters[3]) +
                (parameters[2] + parameters[3]) * Math.Log(1.0 + Math.Pow((x / parameters[1]), parameters[0]));
            return lnumer - ldenom;
        }

        /// <summary>
        /// Unit tests for sampling from normal distribution with ARMS.
        /// </summary>
        /// <param name="mean">Mean value.</param>
        /// <param name="variance">Variance value.</param>
        [TestCase(10.0, 0.1)]
        [TestCase(-5.0, 1.0)]
        [TestCase(0.0, 10.0)]
        [TestCase(10.0, 100.0)]
        public void CompareSampledStats_Normal(double mean, double variance)
        {
            DenseVector parameters = new DenseVector(new double[2] { mean, variance });
            Func<double, double> lnpdf = (double x) => { return lnPDF_Normal(x, parameters); };
            
            double xmin = mean - 4.0 * Math.Sqrt(variance), xmax = mean + 4.0 * Math.Sqrt(variance);
            var arms = new AdaptiveRejectionMetropolisSampler(lnpdf, xmin, xmax);
            double[] sample = arms.Sample(Iteration_NormalTest);
 
            double mean_sim = sample.Average(), variance_sim = 0.0;
            int i = 0;
            while(i < sample.Length)
                variance_sim += Math.Pow((sample[i++] - mean_sim), 2.0);
            variance_sim /= (sample.Length - 1);

            double delta_mean = mean == 0.0 ? Math.Abs(acceptRangeRate) : Math.Abs(acceptRangeRate * mean);
            double delta_sd = Math.Abs(acceptRangeRate * Math.Sqrt(variance));

            Assert.AreEqual(mean, mean_sim, delta_mean);
            Assert.AreEqual(Math.Sqrt(variance), Math.Sqrt(variance_sim), delta_sd);
        }

        // Additional Test
        // let normARMSampler = AdaptiveRejectionMetropolisSampler((normalPdfLn normTheta), -7.0, 9.0, -2.0, 4.0)

        /// <summary>
        /// Unit test for sampling from generalized beta distribution of second kind with ARMS.
        /// The parameters are: (a, b, p, q) = (1.65, 700.0, 2.34, 3.24).
        /// </summary>
        [Test]
        public void CompareSampledStats_GB2()
        {
            DenseVector parameters = new DenseVector(new double[4] { 1.65, 700.0, 2.34, 3.24 });
            Func<double, double> lnpdf = (double x) => { return lnPDF_GB2(x, parameters); };

            double xmin = 0.0, xmax = 20000.0, x1 = 150.0, xn = 2000.0;
            var arms = new AdaptiveRejectionMetropolisSampler(lnpdf, xmin, xmax, x1, xn);
            double[] sample = arms.Sample(Iteration_GB2Test);

            double mean = Integrate.OnClosedInterval((double x) => { return x * Math.Exp(lnpdf(x)); }, 0.0, 20000.0);
            double variance = Integrate.OnClosedInterval((double x) => { return Math.Pow((x - mean), 2.0) * Math.Exp(lnpdf(x)); }, 0.0, 20000.0);

            double mean_sim = sample.Average(), variance_sim = 0.0;
            int i = 0;
            while(i < sample.Length)
                variance_sim += Math.Pow((sample[i++] - mean_sim), 2.0);
            variance_sim /= (sample.Length - 1);

            double delta_mean = Math.Abs(acceptRangeRate * mean), delta_sd = Math.Abs(acceptRangeRate * Math.Sqrt(variance));

            Assert.AreEqual(mean, mean_sim, delta_mean);
            Assert.AreEqual(Math.Sqrt(variance), Math.Sqrt(variance_sim), delta_sd);
        }

        // Additional Test
        //let gbeta2ARMSampler = new AdaptiveRejectionMetropolisSampler((gbeta2PdfLn gbeta2Theta), 0.0, 40000.0, 200.0, 1500.0)

        /// <summary>
        /// Unit tests for comparing two samples whose samplers have same seed.
        /// </summary>
        /// <param name="burnIn">Burn-in times.</param>
        /// <param name="seed">Seed of samplers.</param>
        [TestCase(0,1)]
        [TestCase(1000,100)]
        public void CompareSample_GB2_SameSeed(int burnIn, int seed)
        {
            DenseVector parameters = new DenseVector(new double[4] { 1.65, 700.0, 2.34, 3.24 });
            Func<double, double> lnpdf = (double x) => { return lnPDF_GB2(x, parameters); };

            double xmin = 0.0, xmax = 20000.0, x1 = 150.0, xn = 2000.0;
            var arms1 = new AdaptiveRejectionMetropolisSampler(lnpdf, xmin, xmax, x1, xn, burnIn, seed);
            var arms2 = new AdaptiveRejectionMetropolisSampler(lnpdf, xmin, xmax, x1, xn, burnIn, seed);
            double[] sample1 = arms1.Sample(Iteration_GB2Test);
            double[] sample2 = arms2.Sample(Iteration_GB2Test);

            CollectionAssert.AreEqual(sample1, sample2);
        }
    }
}
