using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Appendix.Optimization;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.OptimizationTests
{
    [TestFixture]
    public class LinearModelEstimationTests
    {
        /// <summary>
        /// Rate of accept range for parameter comparison.
        /// </summary>
        private double acceptRangeRate = 0.01;

        /// <summary>
        /// Two series of data; x and y.
        /// x is a list from 1 to 10.
        /// y is created from the following linear model.
        /// y = 10.0 + 5.0 * x + u
        /// u is distributed according to normal distribution with mean=0 and variance=4.0.
        /// </summary>
        private System.Collections.Generic.List<double> xLM = new System.Collections.Generic.List<double>();
        private System.Collections.Generic.List<double> yLM = new System.Collections.Generic.List<double>();

        /// <summary>
        /// Likelihood function.
        /// </summary>
        /// <param name="parameters">Parameter values.</param>
        /// <returns>Value of log-likelihood.</returns>
        private double likehood(double[] parameters)
        {
            int i = 0;
            double res = 0.0;
            while (i < yLM.Count)
            {
                res += (-0.5) * (System.Math.Log(2.0 * System.Math.PI) + (System.Math.Log(parameters[2]))) -
                    (System.Math.Pow((yLM[i] - parameters[0] - parameters[1] * xLM[i]), 2.0) / (2.0 * parameters[2]));
                i++;
            }

            return res;
        }

        /// <summary>
        /// Set-up section.
        /// Reading data from CSV file.
        /// </summary>
        [SetUp]
        public void SetUp()
        {
            var currentPath = System.AppDomain.CurrentDomain.BaseDirectory;
            System.Environment.CurrentDirectory = currentPath;

            using (var reader = new System.IO.StreamReader(@"..\..\TestData\XYData.csv"))
            {
                while(reader.EndOfStream != true)
                {
                    string oneline = reader.ReadLine();
                    string[] xystr = oneline.Split(new char[1] {','});
                    xLM.Add(Double.Parse(xystr[0]));
                    yLM.Add(double.Parse(xystr[1]));
                }
            }
        }

        /// <summary>
        /// Maximum likelihood estimation of linear model.
        /// The estimate of parameters will be (8.46, 5.45, 6.75).
        /// </summary>
        [Test]
        public void LinearModelMLEstimationTest()
        {
            int nmIter = 100, bfgsIter = 100;
            double nmToler = 1e-3, bfgsToler = 1e-1;

            Func<double[], double> targetFunction = (parameters) => { return (-1.0) * likehood(parameters); };

            NelderMead nm = new NelderMead(targetFunction, nmIter, nmToler);
            BFGS bfgs = new BFGS(targetFunction, bfgsIter, bfgsToler);

            double[] initParameters = new double[3] { 10.0, 5.0, 4.0 };
            double[] expectedParameters = new double[3] { 8.46, 5.45, 6.75 };

            var nmResult = nm.Minimize(initParameters);
            var bfgsResult = bfgs.Minimize(nmResult.Parameters);

            Assert.AreEqual(BFGSResultStatus.Converged, bfgsResult.Status);

            int i = 0;
            double delta = 0.0;
            while (i < bfgsResult.Parameters.Length)
            {
                delta = acceptRangeRate * expectedParameters[i];
                Assert.AreEqual(expectedParameters[i], bfgsResult.Parameters[i], delta);
                i++;
            }
        }
    }

    [TestFixture]
    public class GB2EstimationTests
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
        private double lnPDF_GB2(double x, Vector<double> parameters)
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
        private double pdf_GB2(double x, Vector<double> parameters)
        {
            return Math.Exp(lnPDF_GB2(x, parameters));
        }

        /// <summary>
        /// Likelihood function of GB2.
        /// </summary>
        /// <param name="parameters">Vector value of parameters.</param>
        /// <returns>Value of log-likelihood.</returns>
        private double likelihood(double[] parameters)
        {
            int i = 0;
            double res = 0.0;
            while (i < TestData.Count)
            {
                var paramVec = DenseVector.OfArray(parameters);
                res += lnPDF_GB2(TestData[i], paramVec);
                i++;
            }

            return res;
        }

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
                while(reader.EndOfStream != true)
                    TestData.Add(Double.Parse(reader.ReadLine()));
            }
        }

        /// <summary>
        /// Maximum likelihood method to estimate GB2 distribution from sample in "GB2Sample.csv".
        /// </summary>
        [Test]
        public void GB2MLEstimationTest()
        {
            int nmIter = 100, bfgsIter = 100;
            double nmToler = 1e-3, bfgsToler = 1e-1;

            Func<double[], double> targetfunction = (parameters) => { return (-1.0) * (1e-3) * likelihood(parameters); };
            var nm = new NelderMead(targetfunction, nmIter, nmToler);
            var bfgs = new BFGS(targetfunction, bfgsIter, bfgsToler);

            var nderiv = new MathNet.Numerics.Appendix.Parallel.Differentiation.ParallelNumericalJacobian();
            bfgs.DerivationMethod = (x) => { return nderiv.Evaluate(targetfunction, x); };

            var initParams = new double[4] { 1.65, 700.0, 2.34, 3.24 };
            var expectedParams = new double[4] { 1.65, 700.0, 2.34, 3.24 };

            var nmResult = nm.Minimize(initParams);
            var bfgsResult = bfgs.Minimize(nmResult.Parameters);

            int i = 0;
            double delta = 0.0;
            while(i < bfgsResult.Parameters.Length)
            {
                delta = Math.Abs(acceptRangeRate * expectedParams[i]);
                Assert.AreEqual(expectedParams[i], bfgsResult.Parameters[i], delta);
                i++;
            }

            while (i < bfgs.LatestXVector.Count)
            {
                Assert.AreEqual(bfgsResult.Parameters[i], bfgs.LatestXVector[i]);
                i++;
            }
        }
    }
}
