// Math.NET Numerics Appendix Library License (MIT/X11)
// ===================================================
// 
// Copyright (c) 2013 Fukui Shogo
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

namespace MathNet.Numerics.Differentiation

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double
open System
open System.Threading.Tasks

/// A type which imprements numerical differentiations.
[<CompiledName "DifferentiationFSharp">]
type Differentiation() =
    static let m_DeltaCoarse = MathNet.Numerics.Precision.DoublePrecision**(1.0/3.0)
    static let m_DeltaFine = MathNet.Numerics.Precision.DoublePrecision**(1.0/7.0)

    static member private ShiftX((xs: Vector<float>), (dx: float), index, op) =
        DenseVector.init xs.Count (fun i -> if i = index then
                                                (op xs.[i] dx)
                                            else
                                                xs.[i])

    static member private PartialDerivative_Coarse ((f: Vector<float> -> float), (xs: Vector<float>), index) =
        let xpdx = Differentiation.ShiftX(xs, m_DeltaCoarse, index, (+))
        let xmdx = Differentiation.ShiftX(xs, m_DeltaCoarse, index, (-))
        
        (f(xpdx) - f(xmdx)) / (2.0 * m_DeltaCoarse) 

    static member private PartialDerivative_Fine ((f: Vector<float> -> float), (xs: Vector<float>), index) =
        let deriv_temps = DenseVector.init 3 (fun i ->
            let m = (double (i + 1)) 
            let xpdx = Differentiation.ShiftX(xs, m * m_DeltaFine, index, (+))
            let xmdx = Differentiation.ShiftX(xs, m * m_DeltaFine, index, (-))
            (f(xpdx) - f(xmdx)) / (2.0 * m) )

        (15.0 * deriv_temps.[0] - 6.0 * deriv_temps.[1] + deriv_temps.[2]) / (10.0 * m_DeltaFine)

    static member public Derivative((f: Vector<float> -> float), (xs: Vector<float>), ?fine_on0 : bool, ?parallel_on0 : bool) =
        let fine_on = defaultArg fine_on0 false
        let parallel_on = defaultArg parallel_on0 false

        let getpd i =  if fine_on then
                           Differentiation.PartialDerivative_Fine(f, xs, i)
                       else
                           Differentiation.PartialDerivative_Coarse(f, xs, i)

        let res_array = if parallel_on then Array.Parallel.init xs.Count getpd
                        else Array.init xs.Count getpd

        DenseVector.init xs.Count (fun i -> res_array.[i])
