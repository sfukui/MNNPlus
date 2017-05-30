using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.OptimizationTests
{
    /// <summary>
    /// "LBFGSBSimpleTests" is optimization tests for a simple function; f(x) = x * sin(x), x \in [0, 4].
    /// The function has three local minima.
    /// In this test, we validate whether the gradient is almost 0. 
    /// </summary>
    [TestFixture]
    public class LBFGSBSimpleTests
    {
        /// <summary>
        /// Test function for L-BFGS-B optimization.
        /// The function has three local minima when x is in [0.0, 4.0].
        /// </summary>
        /// <param name="x">Vector value of input variable.</param>
        /// <returns></returns>
        private double TargetFunction(double[] x)
        {
            return 0.1 * x[0] * Math.Sin(x[0]);
        }

        /// <summary>
        /// Differentiation of target function.
        /// If an x is at one of local minima, this function value become 0.0.
        /// </summary>
        private double DerivTargetFunction(double[] x)
        {
            return Math.Sin(x[0]) + x[0] * Math.Cos(x[0]);
        }

        /// <summary>
        /// Bounds of input variable.
        /// </summary>
        private Tuple<double, double>[] m_bounds = { new Tuple<double, double>(0.0, 4.0 * Math.PI) };

        /// <summary>
        /// Rate of accept range for parameter comparison.
        /// </summary>
        private double acceptRangeRate = 1e-3;

        /// <summary>
        /// Maximum iterations and tolerances of L-BFGS-B method.
        /// </summary>
        private int lbfgsbIteration = 100;
        private double lbfgsbTolerance = 1e-4;

        /// <summary>
        /// An object of L-BFGS-B optimization class.
        /// </summary>
        private LBFGSB m_lbfgsb;

        [SetUp]
        public void SetUp()
        {
            m_lbfgsb = new LBFGSB(TargetFunction, m_bounds, lbfgsbIteration, lbfgsbTolerance);
        }

        /// <summary>
        /// Can minimize function by L-BFGS-B method.
        /// </summary>
        /// <param name="xInit">Initial value of x[0].</param>
        [TestCase(0.0)]
        [TestCase(Math.PI)]
        [TestCase(2.0 * Math.PI)]
        [TestCase(3.0 * Math.PI)]
        [TestCase(4.0 * Math.PI)]
        public void LBFGSBOptimization(double xInit)
        {
            double[] init = new double[1] { xInit };
            double expectedDValue = 0.0;
            var result = m_lbfgsb.Minimize(init);

            Assert.AreEqual(LBFGSBResultStatus.Converged, result.Status);
            double actualDValue = DerivTargetFunction(result.Values);
            Assert.AreEqual(expectedDValue, actualDValue, acceptRangeRate);
        }
    }

    /// <summary>
    /// "LBFGSBGB2EstimationTests" is a test of density estimation of generalized beta of second kind.  
    /// </summary>
    [TestFixture]
    public class LBFGSBGB2EstimationTests
    {
        /// <summary>
        /// Rate of accept range for parameter comparison.
        /// </summary>
        private double acceptRangeRate = 0.1;

        /// <summary>
        /// Dataset sampled from generarized beta distribution of second kind(GB2).
        /// The parameters (\alpha, \beta, p, q) = (1.65, 700.0, 2.34, 3.24).
        /// </summary>
        private System.Collections.Generic.List<double> TestData = new System.Collections.Generic.List<double>();

        /// <summary>
        /// Log of GB2 density function.
        /// </summary>
        /// <param name="x">Value of variable.</param>
        /// <param name="parameters">Vector value of parameters.</param>
        /// <returns>Natural log of density.</returns>
        private double LnPDF_GB2(double x, Vector<double> parameters)
        {
            double lnumer = Math.Log(parameters[0]) + (parameters[0] * parameters[2] - 1.0) * Math.Log(x);
            double ldenom = (parameters[0] * parameters[2]) * Math.Log(parameters[1]) +
                SpecialFunctions.BetaLn(parameters[2], parameters[3]) +
                (parameters[2] + parameters[3]) * Math.Log(1.0 + Math.Pow((x / parameters[1]), parameters[0]));
            return lnumer - ldenom;
        }

        /// <summary>
        /// GB2 density function.
        /// </summary>
        /// <param name="x">Value of variable.</param>
        /// <param name="parameters">Vector value of parameters.</param>
        /// <returns>Value of density.</returns>
        private double PDF_GB2(double x, Vector<double> parameters)
        {
            return Math.Exp(LnPDF_GB2(x, parameters));
        }

        /// <summary>
        /// Likelihood function of GB2.
        /// </summary>
        /// <param name="parameters">Vector value of parameters.</param>
        /// <returns>Value of log-likelihood.</returns>
        private double Likelihood(double[] parameters)
        {
            int i = 0;
            double res = 0.0;
            while (i < TestData.Count)
            {
                var paramVec = DenseVector.OfArray(parameters);
                res += LnPDF_GB2(TestData[i], paramVec);
                i++;
            }

            return res;
        }

        /// <summary>
        /// Bounds of input variable.
        /// </summary>
        private Tuple<double, double>[] m_Bounds = { new Tuple<double, double>(1e-5, Double.PositiveInfinity),
                                                     new Tuple<double, double>(1e-5, Double.PositiveInfinity),
                                                     new Tuple<double, double>(1e-5, Double.PositiveInfinity),
                                                     new Tuple<double, double>(1e-5, Double.PositiveInfinity)};

        /// <summary>
        /// Maximum iterations and tolerances of Nelder-Mead and L-BFGS-B method.
        /// </summary>
        private int nmIteration = 100;
        private double nmTolerance = 1e-4;
        private int lbfgsbIteration = 100;
        private double lbfgsbTolerance = 1e-2;

        /// <summary>
        /// Target fuction.
        /// </summary>
        private System.Func<double[], double> m_TargetFunc;

        /// <summary>
        /// Objects of optimization class.
        /// </summary>
        private NelderMead m_nm;
        private LBFGSB m_lbfgsb;

        /// <summary>
        /// Set-Up section.
        /// Reading data from CSV file.
        /// </summary>
        [SetUp]
        public void SetUp()
        {
            var currentPath = System.AppDomain.CurrentDomain.BaseDirectory;
            System.Environment.CurrentDirectory = currentPath;

            using (var reader = new System.IO.StreamReader(@"..\..\TestData\GB2Sample.csv"))
            {
                while (reader.EndOfStream != true)
                    TestData.Add(Double.Parse(reader.ReadLine()));
            }

            double scale = 1.0e-3;
            m_TargetFunc =
                    (parameters) => { return (-1.0) * scale * Likelihood(parameters); };

            m_nm = new NelderMead(m_TargetFunc, nmIteration, nmTolerance);
            m_lbfgsb = new LBFGSB(m_TargetFunc, m_Bounds, lbfgsbIteration, lbfgsbTolerance);
        }

        /// <summary>
        /// Maximum likelihood method to estimate GB2 distribution from sample in "GB2Sample.csv".
        /// </summary>
        [Test]
        public void GB2MLEstimationTest()
        {
            var initParams = new double[4] { 1.65, 700.0, 2.34, 3.24 };
            var expectedParams = new double[4] { 1.65, 700.0, 2.34, 3.24 };

            var nmResult = m_nm.Minimize(initParams);
            var lbfgsbResult = m_lbfgsb.Minimize(nmResult.Parameters);

            int i = 0;
            double delta = 0.0;
            while (i < lbfgsbResult.Values.Length)
            {
                delta = Math.Abs(acceptRangeRate * expectedParams[i]);
                Assert.AreEqual(expectedParams[i], lbfgsbResult.Values[i], delta);
                i++;
            }
        }
    }

}