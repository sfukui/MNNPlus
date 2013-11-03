﻿// <copyright file="WH2006.cs" company="Math.NET">
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
    /// Wichmann-Hill’s 2006 combined multiplicative congruential generator. 
    /// </summary>
    /// <remarks>See: Wichmann, B. A. &amp; Hill, I. D. (2006), "Generating good pseudo-random numbers".
    /// Computational Statistics &amp; Data Analysis 51:3 (2006) 1614-1622
    /// </remarks>
    public class WH2006 : AbstractRandomNumberGenerator
    {
        private const uint Modw = 2147483123;
        private const double ModwRecip = 1.0/Modw;
        private const uint Modx = 2147483579;
        private const double ModxRecip = 1.0/Modx;
        private const uint Mody = 2147483543;
        private const double ModyRecip = 1.0/Mody;
        private const uint Modz = 2147483423;
        private const double ModzRecip = 1.0/Modz;
        private ulong _wn = 1;
        private ulong _xn;
        private ulong _yn = 1;
        private ulong _zn = 1;

        /// <summary>
        /// Initializes a new instance of the <see cref="WH2006"/> class using
        /// the current time as the seed.
        /// </summary>
        public WH2006() : this((int) DateTime.Now.Ticks)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="WH2006"/> class using
        /// the current time as the seed.
        /// </summary>
        /// <param name="threadSafe">if set to <c>true</c> , the class is thread safe.</param>
        public WH2006(bool threadSafe)
            : this((int) DateTime.Now.Ticks, threadSafe)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="WH2006"/> class.
        /// </summary>
        /// <param name="seed">The seed value.</param>
        /// <remarks>If the seed value is zero, it is set to one. Uses the
        /// value of <see cref="Control.ThreadSafeRandomNumberGenerators"/> to
        /// set whether the instance is thread safe.</remarks>
        public WH2006(int seed) : this(seed, Control.ThreadSafeRandomNumberGenerators)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="WH2006"/> class.
        /// </summary>
        /// <param name="seed">The seed value.</param>
        /// <remarks>The seed is set to 1, if the zero is used as the seed.</remarks>
        /// <param name="threadSafe">if set to <c>true</c> , the class is thread safe.</param>
        public WH2006(int seed, bool threadSafe)
            : base(threadSafe)
        {
            if (seed == 0)
            {
                seed = 1;
            }
            _xn = (uint) seed%Modx;
        }

        /// <summary>
        /// Returns a random number between 0.0 and 1.0.
        /// </summary>
        /// <returns>
        /// A double-precision floating point number greater than or equal to 0.0, and less than 1.0.
        /// </returns>
        protected override double DoSample()
        {
            _xn = 11600*_xn%Modx;
            _yn = 47003*_yn%Mody;
            _zn = 23000*_zn%Modz;
            _wn = 33000*_wn%Modw;

            double u = _xn*ModxRecip + _yn*ModyRecip + _zn*ModzRecip + _wn*ModwRecip;
            u -= (int) u;
            return u;
        }
    }
}