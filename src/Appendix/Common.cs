using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace MathNet.Numerics.Appendix.Common
{
    public static class CommonTools
    {
        public static double? FSOptionToCSNullable(Microsoft.FSharp.Core.FSharpOption<double> optionValue)
        {
            double? res = null;
            if (Microsoft.FSharp.Core.FSharpOption<double>.get_IsSome(optionValue)) res = optionValue.Value;
            return res;
        }

        public static Vector<double> FSOptionVectorToCSVector(Microsoft.FSharp.Core.FSharpOption<Vector<double>> optionValue)
        {
            if (Microsoft.FSharp.Core.FSharpOption<Vector<double>>.get_IsNone(optionValue)) return null;
            else return optionValue.Value;
        }

        public static Vector<double> FSOptionVectorToCSVector(Microsoft.FSharp.Core.FSharpOption<DenseVector> optionValue)
        {
            if (Microsoft.FSharp.Core.FSharpOption<DenseVector>.get_IsNone(optionValue)) return null;
            else return optionValue.Value;
        }

        public static Matrix<double> FSOptionMatrixToCSMatrix(Microsoft.FSharp.Core.FSharpOption<Matrix<double>> optionValue)
        {
            if (Microsoft.FSharp.Core.FSharpOption<Matrix<double>>.get_IsNone(optionValue)) return null;
            else return optionValue.Value;
        }
    }
}
