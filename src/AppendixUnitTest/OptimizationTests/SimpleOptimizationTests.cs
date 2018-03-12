using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Appendix.Optimization;
using NUnit.Framework;

namespace AppendixUnitTest.OptimizationTests
{
    /// <summary>
    /// SimpleOptimizationTests の概要の説明
    /// </summary>
    [TestFixture]
    public class SimpleOptimizationTests
    {
        /// <summary>
        /// Test function for Nelder-Meand and BFGS optimization.
        /// The function is minimized at (x[0], x[1]) = (1.0, 1.0).
        /// </summary>
        /// <param name="x">Vector value of independent variable.</param>
        /// <returns></returns>
        private double targetFunction(double[] x)
        {
            return (x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 1.0) * (x[1] - 1.0);
        }

        /// <summary>
        /// Rate of accept range for parameter comparison.
        /// </summary>
        private double acceptRangeRate = 0.01;

        /// <summary>
        /// Maximum iterations and tolerances of Nelder-Mead method and BFGS method.
        /// </summary>
        private int nmIteration = 100;
        private int bfgsIteration = 100;
        private double nmTolerance = 1e-5;
        private double bfgsTolerance = 1e-3;


        /// <summary>
        /// Can minimize f(x[0], x[1]) = (x[0] - 1)^2 + x([1] - 1)^2 by Nelder-Mead method.
        /// </summary>
        /// <param name="x0Init">Initial value of x[0].</param>
        /// <param name="x1Init">Initial value of x[1].</param>
        [TestCase(1.0, 1.0)]
        [TestCase(0.0, 0.0)]
        [TestCase(-1.0, 0.0)]
        [TestCase(0.0, -1.0)]
        public void NelderMeadOptimization(double x0Init, double x1Init)
        {
            NelderMead nm = new NelderMead(targetFunction, nmIteration, nmTolerance);

            double[] initParams = new double[2] { x0Init, x1Init };
            double[] expectedParams = new double[2] { 1.0, 1.0 };
            var intermediateresult = nm.Minimize(initParams);
            var result = nm.Minimize(intermediateresult.Parameters);

            Assert.AreEqual(true, result.Converged);

            int i = 0;
            while (i < result.Parameters.Length)
            {
                double delta = acceptRangeRate * expectedParams[i];
                Assert.AreEqual(expectedParams[i], result.Parameters[i], delta);
                i++;
            }
        }


        /// <summary>
        /// Can minimize f(x[0], x[1]) = (x[0] - 1)^2 + x([1] - 1)^2 by BFGS method.
        /// </summary>
        /// <param name="x0Init">Initial value of x[0].</param>
        /// <param name="x1Init">Initial value of x[1].</param>
        [TestCase(1.0, 1.0)]
        [TestCase(0.0, 0.0)]
        [TestCase(-1.0, 0.0)]
        [TestCase(0.0, -1.0)]
        public void BFGSOptimization(double x0Init, double x1Init)
        {
            BFGS bfgs = new BFGS(targetFunction, bfgsIteration, bfgsTolerance);

            double[] initParams = new double[2] { x0Init, x1Init };
            double[] expectedParams = new double[2] { 1.0, 1.0 };
            var result = bfgs.Minimize(initParams);

            Assert.AreEqual(BFGSResultStatus.Converged, result.Status);

            int i = 0;
            while (i < result.Parameters.Length)
            {
                double delta = acceptRangeRate * expectedParams[i];
                Assert.AreEqual(expectedParams[i], result.Parameters[i], delta);
                i++;
            }
        }
    }
}
