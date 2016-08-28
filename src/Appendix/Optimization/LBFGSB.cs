using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using MathNet.Numerics.Appendix;

namespace MathNet.Numerics.Optimization
{
    public class BoundedFunction
    {
        private BoundedFunctionFSharp m_BoundedFunction;

        public BoundedFunction(System.Func<double[], double> f, Tuple<double, double>[] bounds)
        {
            m_BoundedFunction = new BoundedFunctionFSharp(f, bounds);
        }

        public int ParameterNum
        {
            get { return m_BoundedFunction.ParameterNum; }
        }

        public Tuple<double, double>[] Bounds
        {
            get { return m_BoundedFunction.Bounds; }
        }

        public System.Func<double[], double> Func
        {
            get { return m_BoundedFunction.Func; }
        }

        internal BoundedFunctionFSharp BoundedFunc
        {
            get { return m_BoundedFunction; }
        }

        public double EvaluateRaw(double[] values)
        {
            return m_BoundedFunction.EvaluateRaw(values);
        }

        public double Evaluate(double[] values)
        {
            return m_BoundedFunction.Evaluate(values);
        }
    }

    public class LBFGSBResult
    {
        public LBFGSBResultStatus Status { get; private set; }
        public double[] Parameters { get; private set; }
        public double? FunctionValue { get; private set; }
        public double[,] InvertedWeightMatrix { get; private set; }

        public LBFGSBResult(LBFGSBResultFSharp result)
        {
            Status = (LBFGSBResultStatus)(result.Status);
            Parameters = result.Parameters;
            FunctionValue = result.FunctionValue;
            InvertedWeightMatrix = result.InvertedWeightMatrix;
        }
    }

    public enum LBFGSBResultStatus
    {
        Converged = 0,
        NotConverged = 1,
        FunctionValueInvalid = 2,
        GradientInvalid = 3,
        WeightMatrixInvalid = 4,
        LineSearchFailure = 5,
        InvertedWeightMatrixInvalid = 6,
        NoGeneralizedCauchyPoint = 7,
    }

    public class LBFGSB
    {
        private LBFGSBFSharp m_LBFGSBFS;

        public LBFGSB(BoundedFunction boundedFunction, int iteration, double tolerance, int approxDimension)
        {
            m_LBFGSBFS = new LBFGSBFSharp(boundedFunction.BoundedFunc, iteration, tolerance, approxDimension);
        }

        public LBFGSB(System.Func<double[], double> f, Tuple<double, double>[] bounds, int iteration, double tolerance, int approxDimension)
        {
            m_LBFGSBFS = new LBFGSBFSharp(f, bounds, iteration, tolerance, approxDimension);
        }

        public LBFGSB(BoundedFunction boundedFunction, int iteration, double tolerance)
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