using System;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics.Mcmc;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.MCMCTests
{
    [TestFixture]
    public class ARMSSamplingTests
    {
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

        [Test]
        public void CompareSampledStats_Normal()
        {

        }

        [Test]
        public void CompareSampledStats_GB2()
        {
        }
    }
}
