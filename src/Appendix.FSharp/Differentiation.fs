﻿// Math.NET Numerics Appendix Library License (MIT/X11)
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

[<CompiledName "GradientResultFSharp">]
type GradientResult =
    | Result of float
    | PositiveInfinity
    | NegativeInfinity
    | NaN
    | NoDenominator

[<CompiledName "SearchDenomResultFSharp">]
type SearchDenomResult =
    | Success of GradientResult array
    | Failure

[<CompiledName "DifferentiationFSharp">]
type Differentiation() =
    static let mutable m_InitialDenominator = MathNet.Numerics.Precision.DoubleMachinePrecision**(1.0/3.0)
    static let mutable m_ExtrapolationLength = 4
    static let mutable m_ZeroValue = 100.0 * System.Double.Epsilon
    static let mutable m_SearchTimeOfMaxInitalDenominator = 5
    static let mutable m_CriterionTimeToZeroValue = 3

    static member InitialDenominator
        with get() = m_InitialDenominator
        and set(value) = m_InitialDenominator <- value

    static member ExtrapolationLength
        with get() = m_ExtrapolationLength
        and set(value) = m_ExtrapolationLength <- value

    static member ZeroValue
        with get() = m_ZeroValue
        and set(value) = m_ZeroValue <- value

    static member SearchTimeOfMaxInitalDenominator
        with get() = m_SearchTimeOfMaxInitalDenominator
        and set(value) = m_SearchTimeOfMaxInitalDenominator <- value

    static member CriterionTimeToZeroValue
        with get() = m_CriterionTimeToZeroValue
        and set(value) = m_CriterionTimeToZeroValue <- value

    static member private oneGradient ((f: Vector<float> -> float), (xs: Vector<float>), h, index) =
        let largeF = Vector.mapi (fun j x -> if j = index then x + h else x) xs |> f

        if System.Double.IsNaN(largeF) then NaN
        else
            let res = 0.5 *
                      ( largeF - (Vector.mapi (fun j x -> if j = index then x - h else x) xs |> f) ) / h

            if System.Double.IsNaN(res) then NaN
            else if System.Double.IsPositiveInfinity(res) then PositiveInfinity
            else if System.Double.IsNegativeInfinity(res) then NegativeInfinity
            else Result(res)

    static member private searchInitial ((f: Vector<float> -> float), (xs: Vector<float>), index) = 
        let getGradientCandidates (initialGradients : (float * GradientResult) array) (initialDenominator : float) (length : int) =
            let newGradients = Array.Parallel.init (length - initialGradients.Length) (fun i -> let h = initialDenominator * 10.0**(float i)
                                                                                                (h, Differentiation.oneGradient(f, xs, h, index)) ) |> Array.rev
            Array.append newGradients initialGradients

        let isAnyInvalid gradients = Array.tryFindIndex (fun (h, res) -> match res with
                                                                         | PositiveInfinity | NegativeInfinity | NaN -> true
                                                                         | _ -> false) gradients
        let isAllZero gradients = Array.forall (fun (h, res) -> match res with
                                                                | Result(v) -> abs(v) < Differentiation.ZeroValue
                                                                | _ -> false) gradients 

        let rec search gradientsCandidates i zeroNum =
            if i > Differentiation.SearchTimeOfMaxInitalDenominator then Failure
            else if zeroNum > Differentiation.CriterionTimeToZeroValue then Array.Parallel.init Differentiation.ExtrapolationLength (fun i -> Result(0.0)) |> Success
            else match (isAnyInvalid gradientsCandidates) with
                 | Some(j) -> let nextInitGradients = Array.sub gradientsCandidates 0 j
                              let nextInitDenominator = (fst gradientsCandidates.[0]) * 10.0
                              let nextCandidates = getGradientCandidates nextInitGradients nextInitDenominator Differentiation.ExtrapolationLength
                              search nextCandidates (i+1) zeroNum
                 | None -> if isAllZero gradientsCandidates then let nextInitDenominator = (fst gradientsCandidates.[0]) * 10.0
                                                                 let nextCandidates = getGradientCandidates Array.empty<(float * GradientResult)> nextInitDenominator Differentiation.ExtrapolationLength
                                                                 search nextCandidates (i + 1) (zeroNum + 1)
                           else let (hs, grads) = Array.unzip gradientsCandidates
                                Success(Array.rev grads)
        
        let initialGradients = getGradientCandidates Array.empty<(float * GradientResult)> Differentiation.InitialDenominator Differentiation.ExtrapolationLength       
        search initialGradients 0 0

    static member private oneGradient ((f: Vector<float> -> float), (xs: Vector<float>), index) =
        let rec richardsonExtrapolation (fds: float array) trial =
            let getOneExtrapolation (smallerHG, greaterHG) =
                smallerHG + (smallerHG - greaterHG) / (10.0**(float trial) - 1.0)

            let hGTuples = Array.zip (Array.sub fds 0 (fds.Length-1)) (Array.sub fds 1 (fds.Length-1)) 
            let newDs = Array.Parallel.map getOneExtrapolation hGTuples
            if Array.length newDs = 1 then newDs.[0] else richardsonExtrapolation newDs (trial+1)
        
        let getGradientResultValue g =
                match g with
                | Result(x) -> x
                | _ -> failwith "Improbable input value."

        let searchResult = Differentiation.searchInitial (f, xs, index)
        
        let res = match searchResult with
                  | Failure -> NoDenominator
                  | Success(initDs) -> let initDsVal = Array.Parallel.map getGradientResultValue initDs
                                       Result(richardsonExtrapolation initDsVal 1)

        res
    
    static member private gradientResultToDouble gr =
        match gr with
        | Result(x) -> x
        | PositiveInfinity -> Double.PositiveInfinity
        | NegativeInfinity -> Double.NegativeInfinity
        | _ -> Double.NaN

    static member Gradient ((f: Vector<float> -> float), xs: Vector<float>, hs: Vector<float>) =
        let tRes = Array.Parallel.init xs.Count (fun i -> Differentiation.oneGradient(f, xs, hs.[i], i))
        DenseVector.init tRes.Length (fun i -> Differentiation.gradientResultToDouble tRes.[i])
        
    static member Gradient ((f: Vector<float> -> float), xs: Vector<float>) =
        let tRes = Array.Parallel.init xs.Count (fun i -> Differentiation.oneGradient(f, xs, i))
        DenseVector.init tRes.Length (fun i -> Differentiation.gradientResultToDouble tRes.[i])

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


                