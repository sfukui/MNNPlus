using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Appendix.Optimization;
using MathNet.Numerics.Appendix;

namespace MathNet.Numerics.Appendix.Parallel.Integration
{
    /// <summary>
    /// Numeric Integration (Quadrature).
    /// The codes are copied from "MathNet.Numerics.Parallel.Integration.Integration.cs"
    /// </summary>
    public static class Integrate
    {
        /// <summary>
        /// Approximation of the definite integral of an analytic smooth function on a closed interval.
        /// </summary>
        /// <param name="f">The analytic smooth function to integrate.</param>
        /// <param name="intervalBegin">Where the interval starts, inclusive and finite.</param>
        /// <param name="intervalEnd">Where the interval stops, inclusive and finite.</param>
        /// <param name="targetAbsoluteError">The expected relative accuracy of the approximation.</param>
        /// <returns>Approximation of the finite integral in the given interval.</returns>
        public static double OnClosedInterval(Func<double, double> f, double intervalBegin, double intervalEnd, double targetAbsoluteError)
        {
            return DoubleExponentialTransformation.Integrate(f, intervalBegin, intervalEnd, targetAbsoluteError);
        }

        /// <summary>
        /// Approximation of the definite integral of an analytic smooth function on a closed interval.
        /// </summary>
        /// <param name="f">The analytic smooth function to integrate.</param>
        /// <param name="intervalBegin">Where the interval starts, inclusive and finite.</param>
        /// <param name="intervalEnd">Where the interval stops, inclusive and finite.</param>
        /// <returns>Approximation of the finite integral in the given interval.</returns>
        public static double OnClosedInterval(Func<double, double> f, double intervalBegin, double intervalEnd)
        {
            return DoubleExponentialTransformation.Integrate(f, intervalBegin, intervalEnd, 1e-8);
        }
    }
}
