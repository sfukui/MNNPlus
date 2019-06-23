using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Appendix.Statistics.Mcmc;

namespace MathNet.Numerics.Appendix.Statistics.Mcmc
{

    /// <summary>
    /// Defines adaptive rejection Metropolis sampler class.
    /// </summary>
    public class AdaptiveRejectionMetropolisSampler
    {
        private AdaptiveRejectionMetropolisSamplerFSharp m_ARMSFS;

        /// <summary>
        /// Initializes a new instance of AdaptiveRejectionMetropolisSampler class. 
        /// The initial abscissas consist of three values; x1, (x1 + xn) / 2, and xn.
        /// </summary>
        /// <param name="lnPdf">A log of probability density function(PDF).</param>
        /// <param name="xMin">The minimum value of domain.</param>
        /// <param name="xMax">The maximum value of domain.</param>
        /// <param name="x1">The minimum value of initial abscissas.</param>
        /// <param name="xn">The maximum value of initial abscissas.</param>
        /// <param name="generator">The randome value sampler.</param>
        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, MathNet.Numerics.Random.RandomSource generator)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(lnPdf, xMin, xMax, x1, xn, generator);
        }

        /// <summary>
        /// Initializes a new instance of AdaptiveRejectionMetropolisSampler class.
        /// The initial abscissas consist of three values; x1, (x1 + xn) / 2, and xn.
        /// The random number generator is Mersenne twister.
        /// </summary>
        /// <param name="lnPdf">A log of probability density function(PDF).</param>
        /// <param name="xMin">The minimum value of domain.</param>
        /// <param name="xMax">The maximum value of domain.</param>
        /// <param name="x1">The minimum value of initial abscissas.</param>
        /// <param name="xn">The maximum value of initial abscissas.</param>
        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax,
        double x1, double xn)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(lnPdf, xMin, xMax, x1, xn);
        }

        /// <summary>
        /// Initializes a new instance of AdaptiveRejectionMetropolisSampler class. 
        /// The initial abscissas consist of three values; x1, (x1 + xn) / 2, and xn.
        /// The random number generator is Mersenne twister.
        /// </summary>
        /// <param name="lnPdf">A log of probability density function(PDF).</param>
        /// <param name="xMin">The minimum value of domain.</param>
        /// <param name="xMax">The maximum value of domain.</param>
        /// <param name="x1">The minimum value of initial abscissas.</param>
        /// <param name="xn">The maximum value of initial abscissas.</param>
        /// <param name="seed">The seed of Mersenne twister.</param>
        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax,
            double x1, double xn, int seed)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(lnPdf, xMin, xMax, x1, xn, seed);
        }

        /// <summary>
        /// Initializes a new instance of AdaptiveRejectionMetropolisSampler class. The random number generator is Mersenne twister. 
        /// The minimum value of initial abscissas(x1) is $\max(\mu - 2 \times \sigma, (xMin + \mu) / 2)$, the maximum value(xn) is $\min(\mu + 2 \times \sigma, (\mu + xMax) / 2)$,
        /// and the central value is (x1 + xn) / 2.
        /// $\mu$ is the mean of the random variable and $\sigma$ is the the standard deviation of the random variable, 
        /// The random number generator is Mersenne twister.
        /// </summary>
        /// <param name="lnPdf">A log of probability density function(PDF).</param>
        /// <param name="xMin">The minimum value of domain.</param>
        /// <param name="xMax">The maximum value of domain.</param>
        /// <param name="seed">The seed of Mersenne twister.</param>
        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax, int seed)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(lnPdf, xMin, xMax, seed);
        }

        /// <summary>
        /// Initializes a new instance of AdaptiveRejectionMetropolisSampler class. The random number generator is Mersenne twister. 
        /// The minimum value of initial abscissas(x1) is $\max(\mu - 2 \times \sigma, (xMin + \mu) / 2)$, the maximum value(xn) is $\min(\mu + 2 \times \sigma, (\mu + xMax) / 2)$,
        /// and the central value is (x1 + xn) / 2.
        /// $\mu$ is the mean of the random variable and $\sigma$ is the the standard deviation of the random variable, 
        /// The random number generator is Mersenne twister.
        /// </summary>
        /// <param name="lnPdf">A log of probability density function(PDF).</param>
        /// <param name="xMin">The minimum value of domain.</param>
        /// <param name="xMax">The maximum value of domain.</param>
        public AdaptiveRejectionMetropolisSampler(System.Func<double, double> lnPdf, double xMin, double xMax)
        {
            m_ARMSFS = new AdaptiveRejectionMetropolisSamplerFSharp(lnPdf, xMin, xMax);
        }

        /// <summary>
        /// Samples random numberss with adaptive rejection Metropolis sampler.
        /// </summary>
        /// <param name="x0">The initial value of the random variable.</param>
        /// <param name="iteration">The number of random number sampling.</param>
        /// <returns></returns>
        public double[] Sample(double x0, int iteration)
        {
            return m_ARMSFS.Sample(x0, iteration).ToArray<double>();
        }

        /// <summary>
        /// Samples random numberss with adaptive rejection Metropolis sampler.
        /// The initial value of the random variable is the mean calculated from the probability density function, 
        /// </summary>
        /// <param name="iteration">The number of random number sampling.</param>
        /// <returns></returns>
        public double[] Sample(int iteration)
        {
            return m_ARMSFS.Sample(iteration).ToArray<double>();
        }
    }
}
