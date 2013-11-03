﻿// <copyright file="Mrg32k3a.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
//
// Copyright (c) 2009-2010 Math.NET
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using System;

namespace MathNet.Numerics.Random
{
    /// <summary>
    /// A 32-bit combined multiple recursive generator with 2 components of order 3.
    /// </summary>
    ///<remarks>Based off of P. L'Ecuyer, "Combined Multiple Recursive Random Number Generators," Operations Research, 44, 5 (1996), 816--822. </remarks>
    public class Mrg32k3a : AbstractRandomNumberGenerator
    {
        private const double A12 = 1403580;
        private const double A13 = 810728;
        private const double A21 = 527612;
        private const double A23 = 1370589;
        private const double Modulus1 = 4294967087;
        private const double Modulus2 = 4294944443;

        private const double Reciprocal = 1.0 / Modulus1;
        private double _xn1 = 1;
        private double _xn2 = 1;
        private double _xn3;
        private double _yn1 = 1;
        private double _yn2 = 1;
        private double _yn3 = 1;

        /// <summary>
        /// Initializes a new instance of the <see cref="Mcg31m1"/> class using
        /// the current time as the seed.
        /// </summary>
        /// <remarks>If the seed value is zero, it is set to one. Uses the
        /// value of <see cref="Control.ThreadSafeRandomNumberGenerators"/> to
        /// set whether the instance is thread safe.</remarks>
        public Mrg32k3a() : this((int)DateTime.Now.Ticks)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Mcg31m1"/> class using
        /// the current time as the seed.
        /// </summary>
        /// <param name="threadSafe">if set to <c>true</c> , the class is thread safe.</param>
        public Mrg32k3a(bool threadSafe) : this((int)DateTime.Now.Ticks, threadSafe)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Mrg32k3a"/> class.
        /// </summary>
        /// <param name="seed">The seed value.</param>
        /// <remarks>If the seed value is zero, it is set to one. Uses the
        /// value of <see cref="Control.ThreadSafeRandomNumberGenerators"/> to
        /// set whether the instance is thread safe.</remarks>
        public Mrg32k3a(int seed) : this(seed, Control.ThreadSafeRandomNumberGenerators)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Mrg32k3a"/> class.
        /// </summary>
        /// <param name="seed">The seed value.</param>
        /// <param name="threadSafe">if set to <c>true</c>, the class is thread safe.</param>
        public Mrg32k3a(int seed, bool threadSafe) : base(threadSafe)
        {
            if (seed == 0)
            {
                seed = 1;
            }
            _xn3 = (uint)seed;
        }


        /// <summary>
        /// Returns a random number between 0.0 and 1.0.
        /// </summary>
        /// <returns>
        /// A double-precision floating point number greater than or equal to 0.0, and less than 1.0.
        /// </returns>
        protected override double DoSample()
        {
            double xn = A12*_xn2 - A13*_xn3;
            double k = (long) (xn/Modulus1);
            xn -= k*Modulus1;
            if (xn < 0)
            {
                xn += Modulus1;
            }

            double yn = A21*_yn1 - A23*_yn3;
            k = (long) (yn/Modulus2);
            yn -= k*Modulus2;
            if (yn < 0)
            {
                yn += Modulus2;
            }
            _xn3 = _xn2;
            _xn2 = _xn1;
            _xn1 = xn;
            _yn3 = _yn2;
            _yn2 = _yn1;
            _yn1 = yn;

            if (xn <= yn)
            {
                return (xn - yn + Modulus1)*Reciprocal;
            }
            return (xn - yn)*Reciprocal;
        }
    }
}