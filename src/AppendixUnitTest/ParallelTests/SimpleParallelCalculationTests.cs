﻿using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Parallel.Differentiation;
using MathNet.Numerics.Parallel.Integration;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.ParallelTests
{
    [TestFixture]
    public class SImpleParallelCalculationTests
    {
        /// <summary>
        /// Test function for parallel numerical differentiation.
        /// </summary>
        /// <param name="x">Vector value of independent variable.</param>
        /// <returns>Function value.</returns>
        private double TargetFunction1(Vector<double> x)
        {
            return (x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 1.0) * (x[1] - 1.0);
        }

        /// <summary>
        /// Test function for parallel numerical integration.
        /// </summary>
        /// <param name="x">Value of independent variable.</param>
        /// <returns>Function value.</returns>
        private double TargetFunction2(double x)
        {
            return (x + 2.0) * (x + 1.0); 
        }

        /// <summary>
        /// Analytical indefinite integral of TargetFunction2.
        /// </summary>
        /// <param name="x">Value of independent variable.</param>
        /// <returns>Value of indefinete integral.</returns>
        private double ValidIndefiniteIntegral(double x)
        {
            return (((1.0 / 3.0) * x + (3.0 / 2.0)) * x + 2.0) * x;
        }

        private static double delta = 1e-4; 

        /// <summary>
        /// Can calculate partial derivatives with 2-points method.
        /// </summary>
        /// <param name="x0">Value of x[0].</param>
        /// <param name="x1">Value of x[1].</param>
        /// <param name="d0">Expected partial derivative of \delta f / \delta x[0].</param>
        /// <param name="d1">Expected partial derivative of \delta f / \delta x[1].</param>
        [TestCase(1.0, 1.0, 0.0, 0.0)]
        [TestCase(0.0, 0.0, -2.0, -2.0)]
        [TestCase(-1.0, 0.0, -4.0, -2.0)]
        [TestCase(0.0, -1.0, -2.0, -4.0)]
        [TestCase(10.0, -10.0, 18.0, -22.0)]
        public void NumericalDifferentiation_Coarse1(double x0, double x1, double d0, double d1)
        {
            Vector<double> x = new DenseVector(new double[2] { x0, x1 });
            Vector<double> d = new DenseVector(new double[2] { d0, d1 });
            Vector<double> res = Parallel.Differentiation.Differentiation.Derivative(TargetFunction1, x);

            for(int i = 0; i < x.Count ; ++i)
            {
                Assert.AreEqual(d[i], res[i], delta);
            }
        }

        /// <summary>
        /// Can calculate partial derivatives with 6-points method.
        /// </summary>
        /// <param name="x0">Value of x[0].</param>
        /// <param name="x1">Value of x[1].</param>
        /// <param name="d0">Expected partial derivative of \delta f / \delta x[0].</param>
        /// <param name="d1">Expected partial derivative of \delta f / \delta x[1].</param>
        [TestCase(1.0, 1.0, 0.0, 0.0)]
        [TestCase(0.0, 0.0, -2.0, -2.0)]
        [TestCase(-1.0, 0.0, -4.0, -2.0)]
        [TestCase(0.0, -1.0, -2.0, -4.0)]
        [TestCase(10.0, -10.0, 18.0, -22.0)]
        public void NumericalDifferentiation_Fine1(double x0, double x1, double d0, double d1)
        {
            Vector<double> x = new DenseVector(new double[2] { x0, x1 });
            Vector<double> d = new DenseVector(new double[2] { d0, d1 });
            Vector<double> res = Parallel.Differentiation.Differentiation.Derivative(TargetFunction1, x, true);

            for (int i = 0; i < x.Count; ++i)
            {
                Assert.AreEqual(d[i], res[i], delta);
            }
        }

        /// <summary>
        /// Can calculate numerical integral.
        /// </summary>
        /// <param name="l">Lower bound of independent variable.</param>
        /// <param name="u">Upper bound of independent variable.</param>
        [TestCase(0.0, 1.0)]
        [TestCase(-1.0, 1.0)]
        [TestCase(-10.0, 10.0)]
        [TestCase(-2.0, 4.5)]
        public void NumericalIntegral1(double l, double u)
        {
            double actural = Parallel.Integration.Integrate.OnClosedInterval(TargetFunction2, l, u);
            double expected = ValidIndefiniteIntegral(u) - ValidIndefiniteIntegral(l);

            Assert.AreEqual(expected, actural, delta);
        }
    }
}
