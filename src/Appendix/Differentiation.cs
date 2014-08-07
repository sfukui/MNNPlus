﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using MathNet.Numerics.Appendix;

namespace MathNet.Numerics.Differentiation
{
    public static class Differentiation
    {
        public static Vector<double> Derivative(System.Func<Vector<double>, double> f, Vector<double> xs,
            bool fine_on = false, bool parallel_on = false)
        {
            var fConv = new Converter<Vector<double>, double>(f);
            var fFS = Microsoft.FSharp.Core.FSharpFunc<Vector<double>, double>.FromConverter(fConv);

            var fine_on0 = new Microsoft.FSharp.Core.FSharpOption<bool>(fine_on);
            var parallel_on0 = new Microsoft.FSharp.Core.FSharpOption<bool>(parallel_on);
            return DifferentiationFSharp.Derivative(fFS, xs, fine_on0, parallel_on0);
        }
    }
}
