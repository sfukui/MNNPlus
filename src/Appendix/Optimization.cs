using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using MathNet.Numerics.Appendix;

namespace MathNet.Numerics.Optimization
{
    public class LineSearch
    {
        private LineSearchFSharp m_LineSearchFS;

        public LineSearch(System.Func<Vector<double>, double> f, double xInit, double xMax)
        {
            var fConv = new Converter<Vector<double>, double>(f);
            var fFS = Microsoft.FSharp.Core.FSharpFunc<Vector<double>, double>.FromConverter(fConv);

            m_LineSearchFS = new LineSearchFSharp(fFS, xInit, xMax);
        }

        public double Search(Vector<double> v, Vector<double> d)
        {
            return m_LineSearchFS.Search(v, d);
        }
    }

    public class NelderMeadResult
    {
        public Vector<double> Parameters { get; private set; }
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
        
        public NelderMead(System.Func<Vector<double>, double> f, int iteration, double tolerance)
        {
            var fConv = new Converter<Vector<double>,double>(f);
            var fFS = Microsoft.FSharp.Core.FSharpFunc<Vector<double>, double>.FromConverter(fConv);

            m_NelderMeadFS = new NelderMeadFSharp(fFS, iteration, tolerance);
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

        public NelderMeadResult Minimize(Vector<double> initVal)
        {
            var resFS = m_NelderMeadFS.Minimize(initVal);
            return new NelderMeadResult(m_NelderMeadFS.ResultConvertToType(resFS));
        }
    }

    public class QuasiNewtonMethodResult
    {
        public QuasiNewtonMethodResultStatus Status { get; private set; }
        public Vector<double> Parameters { get; private set; }
        public double? FunctionValue { get; private set; }
        public Matrix<double> InvertedWeightMatrix { get; private set; }

        public QuasiNewtonMethodResult(QuasiNewtonMethodResultFSharp result)
        {
            Status = (QuasiNewtonMethodResultStatus)(result.Status);
            Parameters = result.Parameters;
            FunctionValue = result.FunctionValue;
            InvertedWeightMatrix = result.InvertedWeightMatrix;
        }
    }

    public enum QuasiNewtonMethodResultStatus
    {
        Converged = 0,
        NotConverged = 1,
        FunctionValueInvalid = 2,
        GradientInvalid = 3,
        WeightMatrixInvalid = 4,
        LineSearchFailure = 5,
    }

    public class BFGS
    {
        private BFGSFSharp m_BFGSFS;

        public BFGS(System.Func<Vector<double>, double> f, int iteration, double tolerance)
        {
            var fConv = new Converter<Vector<double>,double>(f);
            var fFS = Microsoft.FSharp.Core.FSharpFunc<Vector<double>, double>.FromConverter(fConv);

            m_BFGSFS = new BFGSFSharp(fFS, iteration, tolerance);
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

        public double InitialStepSize
        {
            get { return m_BFGSFS.InitialStepSize; }
            set { m_BFGSFS.InitialStepSize = value; }
        }

        public double MaxStepSize
        {
            get { return m_BFGSFS.MaxStepSize; }
            set { m_BFGSFS.MaxStepSize = value; }
        }

        public QuasiNewtonMethodResult Minimize(Vector<double> initVal)
        {
            var resFS = m_BFGSFS.Minimize(initVal);
            return new QuasiNewtonMethodResult(m_BFGSFS.FSResultToCSResult(resFS));
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

        public Matrix<double> LatestWeightMatrix
        {
            get
            {
                return Appendix.Common.CommonTools.FSOptionMatrixToCSMatrix(m_BFGSFS.LatestWeightMatrix);

            }
        }
    }
}
