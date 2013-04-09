using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;

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
            return m_NelderMeadFS.FSResultToCSResult(resFS);
        }
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
            return m_BFGSFS.FSResultToCSResult(resFS);
        }
    }
}
