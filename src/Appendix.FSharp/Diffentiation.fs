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
open MathNet.Numerics.LinearAlgebra.Generic
open MathNet.Numerics.LinearAlgebra.Double
open System
open System.Threading.Tasks

type GradientResult =
    | Result of float
    | PositiveInfinity
    | NegativeInfinity
    | NaN
    | NoDenominator

type SearchDenomResult =
    | Success of GradientResult array
    | Failure

type Differentiation() =
    static let mutable m_InitialDenominator = MathNet.Numerics.Precision.DoubleMachinePrecision**(1.0/3.0)
    static let mutable m_ExtrapolationTime = 4
    static let mutable m_ZeroValue = 100.0 * System.Double.Epsilon
    static let mutable m_SearchInitalDenominatorMaxNumber = 10
    static let mutable m_TimeToZeroValue = 3

    static member InitialDenominator
        with get() = m_InitialDenominator
        and set(value) = m_InitialDenominator <- value

    static member ExtrapolationTime
        with get() = m_ExtrapolationTime
        and set(value) = m_ExtrapolationTime <- value

    static member ZeroValue
        with get() = m_ZeroValue
        and set(value) = m_ZeroValue <- value

    static member SearchInitalDenominatorMaxNumber
        with get() = m_SearchInitalDenominatorMaxNumber
        and set(value) = m_SearchInitalDenominatorMaxNumber <- value

    static member TimeToZeroValue
        with get() = m_TimeToZeroValue
        and set(value) = m_TimeToZeroValue <- value

    static member private OneGradient ((f: Vector<float> -> float), (xs: Vector<float>), h, index) =
        let largef = Vector.mapi (fun j x -> if j = index then x + h else x) xs |> f

        if System.Double.IsNaN(largef) then NaN
        else
            let res = 0.5 *
                      ( largef - (Vector.mapi (fun j x -> if j = index then x - h else x) xs |> f) ) / h

            if System.Double.IsNaN(res) then NaN
            else if System.Double.IsPositiveInfinity(res) then PositiveInfinity
            else if System.Double.IsNegativeInfinity(res) then NegativeInfinity
            else Result(res)

    static member private SearchInitial ((f: Vector<float> -> float), (xs: Vector<float>), index) = 
        let rec search h oldg i zeronum =
            if i > Differentiation.SearchInitalDenominatorMaxNumber then Failure
            else if zeronum > Differentiation.TimeToZeroValue then Array.Parallel.init Differentiation.ExtrapolationTime (fun i -> Result(0.0)) |> Success
            else let isinvalid = (fun x -> match x with
                                           | PositiveInfinity | NegativeInfinity | NaN -> true
                                           | _ -> false)
                 let iszero = (fun x -> match x with
                                        | Result(v) -> abs(v) < Differentiation.ZeroValue
                                        | _ -> true)

                 let nexth = (h * 10.0**(float Differentiation.ExtrapolationTime - 1.0))
                 let largestgrad = Differentiation.OneGradient(f, xs, nexth, index)
                 if (isinvalid oldg) then search nexth largestgrad (i+1) zeronum
                 else
                     match largestgrad with
                     | Result(res) -> let gs = Array.Parallel.init Differentiation.ExtrapolationTime (fun j -> if j = 0 then oldg
                                                                                                               else if j = (Differentiation.ExtrapolationTime - 1) then largestgrad
                                                                                                               else Differentiation.OneGradient(f, xs, (h * (10.0**(float j))), index)
                                                                                                     )
                                      // if the following conditions are true, then denominator search is retried.
                                      // * One or more gradient candidates(gs) are NaN or infinite.
                                      // * All candidates are almost zero.                     
                                      if Array.exists isinvalid gs then 
                                          search nexth largestgrad (i+1) zeronum
                                      else if Array.forall iszero gs then
                                          search nexth largestgrad (i+1) (zeronum+1)
                                      else Success(gs)
                     | _ -> search nexth largestgrad (i+1) zeronum
                
        search Differentiation.InitialDenominator (Differentiation.OneGradient(f, xs, Differentiation.InitialDenominator, index)) 0 0

    static member private OneGradient ((f: Vector<float> -> float), (xs: Vector<float>), index) =
        let rec richardsonExtrapolation (fds: float array) trial =
            let getOneExtrapolation (smallerhg, greaterhg) =
                smallerhg + (smallerhg - greaterhg) / (4.0**(float trial) - 1.0)

            let hgtuples = Array.zip (Array.sub fds 0 (fds.Length-1)) (Array.sub fds 1 (fds.Length-1)) 
            let newds = Array.Parallel.map getOneExtrapolation hgtuples
            if Array.length newds = 1 then newds.[0] else richardsonExtrapolation newds (trial+1)
        
        let getGradientResultValue g =
                match g with
                | Result(x) -> x
                | _ -> failwith "Improbable input value."

        let searchresult = Differentiation.SearchInitial (f, xs, index)
        
        let res = match searchresult with
                  | Failure -> NoDenominator
                  | Success(initDs) -> let initDsVal = Array.Parallel.map getGradientResultValue initDs
                                       Result(richardsonExtrapolation initDsVal 1)

        res
    
    static member private GradientResultToDouble gr =
        match gr with
        | Result(x) -> x
        | PositiveInfinity -> Double.PositiveInfinity
        | NegativeInfinity -> Double.NegativeInfinity
        | _ -> Double.NaN

    static member Gradient ((f: Vector<float> -> float), xs: Vector<float>, hs: Vector<float>) =
        let tres = Array.Parallel.init xs.Count (fun i -> Differentiation.OneGradient(f, xs, hs.[i], i))
        DenseVector.init tres.Length (fun i -> Differentiation.GradientResultToDouble tres.[i])
        
    static member Gradient ((f: Vector<float> -> float), xs: Vector<float>) =
        let tres = Array.Parallel.init xs.Count (fun i -> Differentiation.OneGradient(f, xs, i))
        DenseVector.init tres.Length (fun i -> Differentiation.GradientResultToDouble tres.[i])

(*    
    static member Hessian ((f: Vector<float> -> float), xs: Vector<float>, hs: Vector<float>) =
        if xs.Count <> hs.Count then invalidArg "h" "The length of denominator should be equal to that of x."

        let res = new DenseMatrix(xs.Count)
        let xhs i c xs = Vector.mapi (fun j x -> if i = j then (c x hs.[j]) else x) xs

        let aac = new DenseMatrix(xs.Count)
        Matrix.mapi (fun i j _ -> if i = j then let xa = xhs i (+) xs
                                                let xb = xhs i (-) xs
                                                ((f xa) + (f xb) - 2.0 * (f xs)) / hs.[i]**2.0
                                  else if i > j then aac.Item(j,i)
                                  else let xa = xs |> xhs i (+) |> xhs j (+)
                                       let xb = xs |> xhs i (-) |> xhs j (-)
                                       let xa2 = xa |> xhs j (fun x h -> x - 2.0 * h)
                                       let xb2 = xb |> xhs j (fun x h -> x + 2.0 * h)
                                       let res = 0.25 * (f(xa) + f(xb) - f(xa2) - f(xb2)) / (hs.[i] * hs.[j])
                                       do aac.Item(i,j) <- res
                                       res ) res

    static member Hessian ((f: Vector<float> -> float), xs: Vector<float>) =
        let hs = [ for index in 0..(xs.Count - 1) -> ((Differentiation.SearchInitialH(f, xs, (Differentiation.realH xs.[index]), index)) |> fst)**0.5] |> DenseVector.ofList
        Differentiation.Hessian(f, xs, hs)
*)


                
