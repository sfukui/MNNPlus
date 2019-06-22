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
    //public class LineSearch
    //{
    //    private LineSearchFSharp m_LineSearchFS;

    //    public LineSearch(System.Func<double[], double> f, double xInit, double xMax, int trialMax)
    //    {
    //        var fConv = new Converter<double[], double>(f);
    //        var fFS = Microsoft.FSharp.Core.FSharpFunc<double[], double>.FromConverter(fConv);

    //        m_LineSearchFS = new LineSearchFSharp(fFS, xInit, xMax, trialMax);
    //    }

    //    public double C1
    //    {
    //        get { return m_LineSearchFS.C1; }
    //        set { m_LineSearchFS.C1 = value; }
    //    }

    //    public double C2
    //    {
    //        get { return m_LineSearchFS.C2; }
    //        set { m_LineSearchFS.C2 = value; }
    //    }

    //    public double MinimumDifference
    //    {
    //        get { return m_LineSearchFS.MinimumDifference; }
    //        set { m_LineSearchFS.MinimumDifference = value; }
    //    }

    //    public double MaxStepSearchMultiplier
    //    {
    //        get { return m_LineSearchFS.MaxStepSearchMultiplier; }
    //        set { m_LineSearchFS.MaxStepSearchMultiplier = value; }
    //    }

    //    public int MaxStepSearchTrial
    //    {
    //        get { return m_LineSearchFS.MaxStepSearchTrial; }
    //        set { m_LineSearchFS.MaxStepSearchTrial = value; }
    //    }

    //    public int InterpolationSearchTrial
    //    {
    //        get { return m_LineSearchFS.InterpolationSearchTrial; }
    //        set { m_LineSearchFS.InterpolationSearchTrial = value; }
    //    }

    //    public double NextStepMultiplier
    //    {
    //        get { return m_LineSearchFS.NextStepMultiplier; }
    //        set { m_LineSearchFS.NextStepMultiplier = value; }
    //    }

    //    public System.Func<System.Func<double, double>, double, double> NumericalDerivation
    //    {
    //        get { return m_LineSearchFS.NumericalDerivation;  }
    //        set { m_LineSearchFS.NumericalDerivation = value; }
    //    }

    //    public double Search(double[] v, double[] d)
    //    {
    //        return m_LineSearchFS.Search(v, d);
    //    }
    //}

    public class NelderMeadResult
    {
        public double[] Parameters { get; private set; }
        public double FunctionValue { get; private set; }
        public bool Converged { get; private set; }

        public NelderMeadResult(NelderMeadResultFSharp result)
        {
            Parameters = result.Parameters;
            FunctionValue = result.FunctionValue;
            Converged = result.Converged;
        }
    }

    public class NelderMead
    {
        private NelderMeadFSharp m_NelderMeadFS;
        
        public NelderMead(System.Func<double[], double> f, int iteration, double tolerance)
        {
            //var fConv = new Converter<double[], double>(f);
            //var fFS = Microsoft.FSharp.Core.FSharpFunc<double[], double>.FromConverter(fConv);

            m_NelderMeadFS = new NelderMeadFSharp(f, iteration, tolerance);
        }

        public int Iteration
        {
            get { return m_NelderMeadFS.Iteration; }
            set { m_NelderMeadFS.Iteration = value; }
        }

        public double Tolerance
        {
            get { return m_NelderMeadFS.Tolerance; }
            set { m_NelderMeadFS.Tolerance = value; }
        }

        public double ZDelta
        {
            get { return m_NelderMeadFS.ZDelta; }
            set { m_NelderMeadFS.ZDelta = value; }
        }

        public double Delta
        {
            get { return m_NelderMeadFS.Delta; }
            set { m_NelderMeadFS.Delta = value; }
        }

        public double Rho
        {
            get { return m_NelderMeadFS.Rho; }
            set { m_NelderMeadFS.Delta = value; }
        }

        public double Chi
        {
            get { return m_NelderMeadFS.Chi; }
            set { m_NelderMeadFS.Chi = value; }
        }

        public double Psi
        {
            get { return m_NelderMeadFS.Psi; }
            set { m_NelderMeadFS.Psi = value; }
        }

        public double Sigma
        {
            get { return m_NelderMeadFS.Sigma; }
            set { m_NelderMeadFS.Sigma = value; }
        }

        public NelderMeadResult Minimize(double[] initVal)
        {
            var resFS = m_NelderMeadFS.Minimize(initVal);
            return new NelderMeadResult(m_NelderMeadFS.ResultConvertToType(resFS));
        }
    }

    public class BFGSResult
    {
        public BFGSResultStatus Status { get; private set; }
        public double[] Parameters { get; private set; }
        public double? FunctionValue { get; private set; }
        public double[,] InvertedWeightMatrix { get; private set; }

        public BFGSResult(BFGSResultFSharp result)
        {
            Status = (BFGSResultStatus)(result.Status);
            Parameters = result.Parameters;
            FunctionValue = result.FunctionValue;
            InvertedWeightMatrix = result.InvertedWeightMatrix;
        }
    }

    public enum BFGSResultStatus
    {
        InProcess = 0,
        Converged = 1,
        NotConverged = 2,
        FunctionValueInvalid = 3,
        GradientInvalid = 4,
        BFGSMatrixInvalid = 5,
        LineSearchFailure = 6,
        InverseBFGSMatrixInvalid = 7,
    }

    public class BFGS
    {
        private BFGSFSharp m_BFGSFS;

        public BFGS(System.Func<double[], double> f, int iteration, double tolerance)
        {
            m_BFGSFS = new BFGSFSharp(f, iteration, tolerance);
        }

        public int Iteration
        {
            get { return m_BFGSFS.Iteration; }
            set { m_BFGSFS.Iteration = value; }
        }

        public double Tolerance
        {
            get { return m_BFGSFS.Tolerance; }
            set { m_BFGSFS.Tolerance = value; }
        }

        public BFGSResult Minimize(double[] initVal)
        {
            var resFS = m_BFGSFS.Minimize(initVal);
            return new BFGSResult(m_BFGSFS.FSResultToCSResult(resFS));
        }

        public System.Func<double[], double[]> DerivationMethod
        {
            get { return m_BFGSFS.DerivationMethod; }
            set { m_BFGSFS.DerivationMethod = value; }
        }

        public LineSearch LineSearch
        {
            get { return m_BFGSFS.LineSearch; }
            set { m_BFGSFS.LineSearch = value; }
        }

        public double? LatestStepSize
        {
            get
            {
                return Appendix.Common.CommonTools.FSOptionToCSNullable(m_BFGSFS.LatestStepSize);
            }
        }

        public Vector<double> LatestXVector
        {
            get
            {
                return Appendix.Common.CommonTools.FSOptionVectorToCSVector(m_BFGSFS.LatestXVector);
            }
        }

        public Vector<double> LatestGradientVector
        {
            get
            {
                return Appendix.Common.CommonTools.FSOptionVectorToCSVector(m_BFGSFS.LatestGradientVector);
            }
        }

        public Matrix<double> LatestInvertedWeightMatrix
        {
            get
            {
                return Appendix.Common.CommonTools.FSOptionMatrixToCSMatrix(m_BFGSFS.LatestInvertedWeightMatrix);

            }
        }
    }
}
