using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathNet.Numerics.Differentiation.Obsolete
{
    public static class DifferentiationObslolete
    {
        public static double InitialDenominator
        {
            get { return DifferentiationObsoleteFSharp.InitialDenominator; }
            set { DifferentiationObsoleteFSharp.InitialDenominator = value; }
        }

        public static double DenominatorMultiplier
        {
            get { return DifferentiationObsoleteFSharp.DenominatorMultiplier; }
            set { DifferentiationObsoleteFSharp.DenominatorMultiplier = value; }
        }

        public static int ExtrapolationTime
        {
            get { return DifferentiationObsoleteFSharp.ExtrapolationLength; }
            set { DifferentiationObsoleteFSharp.ExtrapolationLength = value; }
        }

        public static double ZeroValue
        {
            get { return DifferentiationObsoleteFSharp.ZeroValue; }
            set { DifferentiationObsoleteFSharp.ZeroValue = value; }
        }

        public static int SearchTimeOfMaxInitalDenominator
        {
            get { return DifferentiationObsoleteFSharp.SearchTimeOfMaxInitalDenominator; }
            set { DifferentiationObsoleteFSharp.SearchTimeOfMaxInitalDenominator = value; }
        }

        public static int CriterionTimeToZeroValue
        {
            get { return DifferentiationObsoleteFSharp.CriterionTimeToZeroValue; }
            set { DifferentiationObsoleteFSharp.CriterionTimeToZeroValue = value; }
        }


    }
}
