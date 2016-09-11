using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace MathNet.Numerics.Appendix.Common
{
    /// <summary>
    /// Common functions bridge a difference between F# and C#.   
    /// </summary>
    public static class CommonTools
    {
        /// <summary>
        /// Converts F# option&gt;double&lt; type to double? type.  
        /// </summary>
        /// <param name="optionValue">A value of option&gt;double&lt; type.</param>
        /// <returns>A nullable double(double?) value.</returns>
        public static double? FSOptionToCSNullable(Microsoft.FSharp.Core.FSharpOption<double> optionValue)
        {
            double? res = null;
            if (Microsoft.FSharp.Core.FSharpOption<double>.get_IsSome(optionValue)) res = optionValue.Value;
            return res;
        }

        /// <summary>
        /// Converts F# option&gt;Vector&gt;double&lt;&lt; type to Vector&gt;double&lt; type.  
        /// </summary>
        /// <param name="optionValue">A value of option&gt;Vector&gt;double&lt;&lt; type</param>
        /// <returns>A Vector&gt;double&lt; value.</returns>
        public static Vector<double> FSOptionVectorToCSVector(Microsoft.FSharp.Core.FSharpOption<Vector<double>> optionValue)
        {
            if (Microsoft.FSharp.Core.FSharpOption<Vector<double>>.get_IsNone(optionValue)) return null;
            else return optionValue.Value;
        }

       /// <summary>
        /// Converts F# option&gt;DenseVector&lt; type to Vector&gt;double&lt; type.  
        /// </summary>
        /// <param name="optionValue">A value of option&gt;DenseVector&lt; type</param>
        /// <returns>A Vector&gt;double&lt; value.</returns>
        public static Vector<double> FSOptionVectorToCSVector(Microsoft.FSharp.Core.FSharpOption<DenseVector> optionValue)
        {
            if (Microsoft.FSharp.Core.FSharpOption<DenseVector>.get_IsNone(optionValue)) return null;
            else return optionValue.Value;
        }

        /// <summary>
        /// Converts F# option&gt;Matrix&gt;double&lt;&lt; type to Matrix&gt;double&lt; type.  
        /// </summary>
        /// <param name="optionValue">A value of option&gt;Matrix&gt;double&lt;&lt; type</param>
        /// <returns>A Matrix&gt;double&lt; value.</returns>
        public static Matrix<double> FSOptionMatrixToCSMatrix(Microsoft.FSharp.Core.FSharpOption<Matrix<double>> optionValue)
        {
            if (Microsoft.FSharp.Core.FSharpOption<Matrix<double>>.get_IsNone(optionValue)) return null;
            else return optionValue.Value;
        }
    }
}
