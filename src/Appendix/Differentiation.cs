using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.Numerics.Differentiation
{
    public static class Differentiation
    {
        public static double InitialDenominator
        {
            get { return DifferentiationFSharp.InitialDenominator; }
            set { DifferentiationFSharp.InitialDenominator = value; }
        }

        public static double DenominatorMultiplier
        {
            get { return DifferentiationFSharp.DenominatorMultiplier; }
            set { DifferentiationFSharp.DenominatorMultiplier = value; }
        }

        public static int ExtrapolationTime
        {
            get { return DifferentiationFSharp.ExtrapolationLength; }
            set { DifferentiationFSharp.ExtrapolationLength = value; }
        }

        public static double ZeroValue
        {
            get { return DifferentiationFSharp.ZeroValue; }
            set { DifferentiationFSharp.ZeroValue = value; }
        }

        public static int SearchTimeOfMaxInitalDenominator
        {
            get { return DifferentiationFSharp.SearchTimeOfMaxInitalDenominator; }
            set { DifferentiationFSharp.SearchTimeOfMaxInitalDenominator = value; }
        }

        public static int CriterionTimeToZeroValue
        {
            get { return DifferentiationFSharp.CriterionTimeToZeroValue; }
            set { DifferentiationFSharp.CriterionTimeToZeroValue = value; }
        }


    }
}
