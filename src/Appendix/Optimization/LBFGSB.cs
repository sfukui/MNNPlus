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

namespace MathNet.Numerics.Appendix.Optimization
{
    /// <summary>
    /// Defines a function with bound(s) of the variable(s).
    /// </summary>
    public class FuncWithBounds
    {
        private FuncWithBoundsFSharp m_FuncWithBounds;

        /// <summary>
        /// Initializes a new instance of BoundedFunction class.
        /// </summary>
        /// <param name="f">A function.</param>
        /// <param name="bounds">Lower and upper bounds of the variables.</param>
        public FuncWithBounds(System.Func<double[], double> f, Tuple<double, double>[] bounds)
        {
            m_FuncWithBounds = new FuncWithBoundsFSharp(f, bounds);
        }

        /// <summary>
        /// Gets the number of variables of the function, 
        /// </summary>
        public int VariableNumber
        {
            get { return m_FuncWithBounds.VariableNumber; }
        }

        /// <summary>
        /// Gets the bounds of the variables.
        /// </summary>
        public Tuple<double, double>[] Bounds
        {
            get { return m_FuncWithBounds.Bounds; }
        }

        /// <summary>
        /// Gets the original function.
        /// </summary>
        public System.Func<double[], double> Func
        {
            get { return m_FuncWithBounds.Func; }
        }

        internal FuncWithBoundsFSharp BoundedFunc
        {
            get { return m_FuncWithBounds; }
        }

        /// <summary>
        /// Evaluates function value without bound checking of the variables.
        /// If one or more values of variables are outside of the bounds, the result of this function is undefined; This may return NaN or throw some exception, or cause the other result. 
        /// </summary>
        /// <param name="values">The values of variables.</param>
        /// <returns>The function value calculated with "values" argument.</returns>
        public double EvaluateRaw(double[] values)
        {
            return m_FuncWithBounds.EvaluateRaw(values);
        }

        /// <summary>
        /// Evaluates function value with bound checking of the variables.
        /// If one or more values of variables are outside of the bounds, this throws System.ArgumentException. 
        /// </summary>
        /// <param name="values">The values of variables.</param>
        /// <returns>The function value calculated with "values" argument.</returns>
        public double Evaluate(double[] values)
        {
            return m_FuncWithBounds.Evaluate(values);
        }
    }

    /// <summary>
    /// Summarizes result of L-BFGS-B optimization.
    /// </summary>
    public class LBFGSBResult
    {
        /// <summary>
        /// L-BFGS-B optimization result status.
        /// </summary>
        /// <remarks> The detail is <see cref="LBFGSBResult"/>. </remarks>
        public LBFGSBResultStatus Status { get; private set; }

        /// <summary>
        /// The values of variables after L-BFGS-B minimization.
        /// </summary>
        public double[] Values { get; private set; }

        /// <summary>
        /// The function value after L-BFGS-B minimization.
        /// </summary>
        public double? FunctionValue { get; private set; }

        /// <summary>
        /// The inverted weight matrix after L-BFGS-B minimization.
        /// </summary>
        public double[,] InverseBFGSMatrix { get; private set; }

        /// <summary>
        /// Initializes new instance of LBFGSBResult class. 
        /// </summary>
        /// <param name="result"></param>
        public LBFGSBResult(LBFGSBResultFSharp result)
        {
            Status = (LBFGSBResultStatus)(result.Status);
            Values = result.Values;
            FunctionValue = result.FunctionValue;
            InverseBFGSMatrix = result.InverseBFGSMatrix;
        }
    }

    /// <summary>
    /// Defines L-BFGS-B optimization result status.
    /// </summary>
    public enum LBFGSBResultStatus
    {
        InProcess = 0,
        Converged = 1,
        NotConverged = 2,
        FunctionValueInvalid = 3,
        GradientInvalid = 4,
        BFGSMatrixInvalid = 5,
        LineSearchFailure = 6,
        InverseBFGSMatrixInvalid = 7,
        NoGeneralizedCauchyPoint = 8,
        BlockOfBFGSMatrixInvalid = 9,
        BlockOfInverseBFGSMatrixInvalid = 10,
        CorrectionHistoryInvalid = 11,
        ConvergedAtCorner = 12,
    }

    public class LBFGSB
    {
        private LBFGSBFSharp m_LBFGSBFS;

        public LBFGSB(FuncWithBounds boundedFunction, int iteration, double tolerance, int approxDimension)
        {
            m_LBFGSBFS = new LBFGSBFSharp(boundedFunction.BoundedFunc, iteration, tolerance, approxDimension);
        }

        public LBFGSB(System.Func<double[], double> f, Tuple<double, double>[] bounds, int iteration, double tolerance, int approxDimension)
        {
            m_LBFGSBFS = new LBFGSBFSharp(f, bounds, iteration, tolerance, approxDimension);
        }

        public LBFGSB(FuncWithBounds boundedFunction, int iteration, double tolerance)
        {
            m_LBFGSBFS = new LBFGSBFSharp(boundedFunction.BoundedFunc, iteration, tolerance);
        }

        public LBFGSB(System.Func<double[], double> f, Tuple<double, double>[] bounds, int iteration, double tolerance)
        {
            m_LBFGSBFS = new LBFGSBFSharp(f, bounds, iteration, tolerance);
        }

        public LBFGSBResult Minimize(double[] initVal)
        {
            var resFS = m_LBFGSBFS.Minimize(initVal);
            return new LBFGSBResult(m_LBFGSBFS.FSResultToCSResult(resFS));
        }

        public System.Func<double[], double[]> DerivationMethod
        {
            get { return m_LBFGSBFS.DerivationMethod; }
            set { m_LBFGSBFS.DerivationMethod = value; }
        }
    }
}
