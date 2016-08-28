using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.OptimizationTests
{
    [TestFixture]
    public class SimpleTests
    {
        /// <summary>
        /// Test function for L-BFGS-B optimization.
        /// The function is minimized at (x[0], x[1]) = (1.0, 1.0).
        /// </summary>
        /// <param name="x">Vector value of independent variable.</param>
        /// <returns></returns>
        private double targetFunction(double[] x)
        {
            return 0.1 * x[0] * Math.Sin(x[0]);
        }

        /// <summary>
        /// Rate of accept range for parameter comparison.
        /// </summary>
        private double acceptRangeRate = 0.01;

        /// <summary>
        /// Maximum iterations and tolerances of L-BFGS-B method.
        /// </summary>
        private int lbfgsbIteration = 100;
        private double lbfgsbTolerance = 1e-4;

        /// <summary>
        /// Can minimize function by L-BFGS-B method.
        /// </summary>
        /// <param name="x0Init">Initial value of x[0].</param>
        /// <param name="x1Init">Initial value of x[1].</param>
        [TestCase(Math.PI)]
        [TestCase(0.0)]
        [TestCase(9.0)]
        [TestCase(4.0 * Math.PI)]
        public void LBFGSBOptimization(double xInit)
        {
            Tuple<double, double>[] bounds = { new Tuple<double, double>(0.0, 4.0 * Math.PI) };
            LBFGSB lbfgsb = new LBFGSB(targetFunction, bounds, lbfgsbIteration, lbfgsbTolerance);

            double[] initParams = new double[1] { xInit };
            double[] expectedParams = new double[1] { 1.0 };
            var result = lbfgsb.Minimize(initParams);

            Assert.AreEqual(LBFGSBResultStatus.Converged, result.Status);

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