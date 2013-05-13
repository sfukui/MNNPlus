using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics.Mcmc;

namespace MathNet.Numerics.Statistics.Mcmc
{
    public class AdaptiveRejectionMetropolisSampler
    {
        private AdaptiveRejectionMetropolisSamplerFSharp m_ARMSFS;

        public AdaptiveRejectionMetropolisSampler(Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, int burnIn, MathNet.Numerics.Random.AbstractRandomNumberGenerator sampler)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn,
                burnIn, sampler);
        }

        public AdaptiveRejectionMetropolisSampler(Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, int burnIn)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn, burnIn);
        }

        public AdaptiveRejectionMetropolisSampler(Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, int burnIn, int seed)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn, burnIn, seed);
        }
        
        public AdaptiveRejectionMetropolisSampler(Func<double, double> lnPdf, double xMin, double xMax, int burnIn, int seed)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, burnIn, seed);
        }

        public AdaptiveRejectionMetropolisSampler(Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax, x1, xn);
        }

        public AdaptiveRejectionMetropolisSampler(Func<double, double> lnPdf, double xMin, double xMax)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(CSFuncToFSFunc(lnPdf), xMin, xMax);
        }

        private Microsoft.FSharp.Core.FSharpFunc<double, double> CSFuncToFSFunc(Func<double, double> func)
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
