using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics.Mcmc;

namespace MathNet.Numerics.Statistics.Mcmc
{
    public class AdaptiveRejectionMetropolisSampler
    {
        private AdaptiveRejectionMetropolisSamplerFSharp m_ARMSFS;

        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, MathNet.Numerics.Random.AbstractRandomNumberGenerator sampler)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn, sampler);
        }

        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn);
        }

        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, int seed)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn, seed);
        }

        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax, int seed)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, seed);
        }

        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax);
        }

        private Microsoft.FSharp.Core.FSharpFunc<double, double> CSFuncToFSFunc(System.Func<double, double> func)
        {
            var fConv = new Converter<double, double>(func);
            return Microsoft.FSharp.Core.FSharpFunc<double, double>.FromConverter(fConv);
        }

        public double[] Sample(double x0, int iteration)
        {
            return m_ARMSFS.Sample(x0, iteration).ToArray<double>();
        }

        public double[] Sample(int iteration)
        {
            return m_ARMSFS.Sample(iteration).ToArray<double>();
        }
    }
}
