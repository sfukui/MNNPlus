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

namespace MathNet.Numerics.Differentiation.Obsolete

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra
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

/// A type which imprements numerical differentiations.
[<CompiledName "DifferentiationObsoleteFSharp">]
type Differentiation() =
    static let mutable m_InitialDenominator = MathNet.Numerics.Precision.DoublePrecision**(1.0/3.0)
    static let mutable m_DenominatorMultiplier = 2.0
    static let mutable m_NumberOfCandidates = 8
    static let mutable m_ExtrapolationLength = 4
    static let mutable m_ZeroValue = 100.0 * System.Double.Epsilon
    static let mutable m_SearchTimeOfMaxInitalDenominator = 5
    static let mutable m_CriterionTimeToZeroValue = 2

    static member InitialDenominator
        with get() = m_InitialDenominator
        and set(value) = m_InitialDenominator <- value

    static member DenominatorMultiplier
        with get() = m_DenominatorMultiplier
        and set(value) = m_DenominatorMultiplier <- value

    static member NumberOfCandidates
        with get() = m_NumberOfCandidates
        and set(value) = m_NumberOfCandidates <- value

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
        let xDiffVec (operator : float -> float -> float) (xVec : Vector<float>) =
            Vector.mapi (fun j x -> if j = index then (operator x h) else x) xVec
        let dxUpperLower = [| (xDiffVec (+) xs) ; (xDiffVec (-) xs) |]
        let fUpperLower = Array.map f dxUpperLower

        let res = 0.5 * (fUpperLower.[0] - fUpperLower.[1]) / h

        if System.Double.IsNaN(res) then NaN
        else if System.Double.IsPositiveInfinity(res) then PositiveInfinity
        else if System.Double.IsNegativeInfinity(res) then NegativeInfinity
        else Result(res)

    static member private searchInitial ((f: Vector<float> -> float), (xs: Vector<float>), index) = 
        let getGradientCandidates (initialGradients : (float * GradientResult) array) (initialDenominator : float) (length : int) =
            let newGradients = Array.init (length - initialGradients.Length) (fun i -> let h = initialDenominator * Differentiation.DenominatorMultiplier**(float i)
                                                                                       (h, Differentiation.oneGradient(f, xs, h, index)) ) |> Array.rev
            Array.append newGradients initialGradients

        let isInvalid grad = match (snd grad) with
                             | PositiveInfinity | NegativeInfinity | NaN -> true
                             | _ -> false
        
        let isZero grad = match (snd grad) with
                          | Result(v) -> abs(v) < Differentiation.ZeroValue
                          | _ -> false

        let rec search gradientsCandidates i zeroNum allValidCands =
            if i > Differentiation.SearchTimeOfMaxInitalDenominator then
                if Array.isEmpty(allValidCands) then Failure
                else Success( (Array.unzip allValidCands) |> snd |> Array.rev )
            else if zeroNum > Differentiation.CriterionTimeToZeroValue then Array.init Differentiation.NumberOfCandidates (fun i -> Result(0.0)) |> Success
            else if (Array.forall isZero gradientsCandidates) then let nextInitDenominator = (fst gradientsCandidates.[0]) * Differentiation.DenominatorMultiplier
                                                                   let nextCandidates = getGradientCandidates Array.empty<(float * GradientResult)> nextInitDenominator Differentiation.NumberOfCandidates
                                                                   search nextCandidates (i + 1) (zeroNum + 1) allValidCands
            else match Array.tryFindIndex (fun x -> (isInvalid x) || (isZero x)) gradientsCandidates with
                 | Some(j) -> let nextInitGradients = Array.sub gradientsCandidates 0 j
                              let nextInitDenominator = (fst gradientsCandidates.[0]) * Differentiation.DenominatorMultiplier
                              let nextCandidates = getGradientCandidates nextInitGradients nextInitDenominator Differentiation.NumberOfCandidates
                              let nextAllValidCands = if Array.exists isInvalid gradientsCandidates then allValidCands else gradientsCandidates
                              search nextCandidates (i+1) zeroNum nextAllValidCands
                 | None -> let (hs, grads) = Array.unzip gradientsCandidates
                           Success(Array.rev grads)
        
        let initialGradients = getGradientCandidates Array.empty<(float * GradientResult)> Differentiation.InitialDenominator Differentiation.NumberOfCandidates      
        search initialGradients 0 0 Array.empty<(float * GradientResult)>

    static member private oneGradient ((f: Vector<float> -> float), (xs: Vector<float>), index) =
        let getGradientResultValue g =
            match g with
            | Result(x) -> x
            | _ -> failwith "Improbable input value."
        
        let getSmallestRangeSet gs =
            let range (gsSub : float array) = ( (Array.max gsSub) - (Array.min gsSub) ) |> abs
            let gsValue = Array.map getGradientResultValue gs
            let gsSubsets = Array.init (Differentiation.NumberOfCandidates - Differentiation.ExtrapolationLength + 1) (fun i -> Array.sub gsValue i Differentiation.ExtrapolationLength)
            let gsSubsetRanges = Array.map range gsSubsets
            let minError = Array.min gsSubsetRanges
            gsSubsets.[(Array.findIndex (fun x -> x = minError) gsSubsetRanges)]

        let rec richardsonExtrapolation (fds: float array) trial =
            let getOneExtrapolation (smallerHG, greaterHG) =
                smallerHG + (smallerHG - greaterHG) / (Differentiation.DenominatorMultiplier**(2.0 * (float trial)) - 1.0)

            let hGTuples = Array.zip (Array.sub fds 0 (fds.Length-1)) (Array.sub fds 1 (fds.Length-1)) 
            let newDs = Array.map getOneExtrapolation hGTuples
            if Array.length newDs = 1 then newDs.[0] else richardsonExtrapolation newDs (trial+1)
        
        let searchResult = Differentiation.searchInitial (f, xs, index)
        
        let res = match searchResult with
                  | Failure -> NoDenominator
                  | Success(initDs) -> let gsSets = getSmallestRangeSet initDs 
                                       Result(richardsonExtrapolation gsSets 1)

        res
    
    static member private gradientResultToDouble gr =
        match gr with
        | Result(x) -> x
        | PositiveInfinity -> Double.PositiveInfinity
        | NegativeInfinity -> Double.NegativeInfinity
        | _ -> Double.NaN

    static member Gradient ((f: Vector<float> -> float), xs: Vector<float>, hs: Vector<float>) =
        let tRes = Array.init xs.Count (fun i -> Differentiation.oneGradient(f, xs, hs.[i], i))
//        let tRes = Array.Parallel.init xs.Count (fun i -> Differentiation.oneGradient(f, xs, hs.[i], i))
        DenseVector.init tRes.Length (fun i -> Differentiation.gradientResultToDouble tRes.[i])
        
    static member Gradient ((f: Vector<float> -> float), xs: Vector<float>) =
        let tRes = Array.init xs.Count (fun i -> Differentiation.oneGradient(f, xs, i))
//        let tRes = Array.Parallel.init xs.Count (fun i -> Differentiation.oneGradient(f, xs, i))
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


                
