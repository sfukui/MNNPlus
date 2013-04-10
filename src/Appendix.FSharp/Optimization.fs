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

namespace MathNet.Numerics.Optimization

open MathNet.Numerics
open MathNet.Numerics.Differentiation
open MathNet.Numerics.Interpolation.Algorithms
open MathNet.Numerics.LinearAlgebra.Generic
open MathNet.Numerics.LinearAlgebra.Double

[<CompiledName "LineSearchFSharp">]
type LineSearch (f: (Vector<float> -> float), xInit, xMax) =
    let c1, c2 = 10.0**(-4.0), 0.9
    let dMin = 1e-16
    let rescueStepSize = 1.0
    let searchMaxStepMult = 0.5
    let searchMaxStepMax = 100

    let cubicInterpolation((a_old, f_old, df_old), (a_cur, f_cur, df_cur)) =
        if abs( (a_old - a_cur) / a_cur) < dMin then Some(a_cur)
        else
            let d1 = df_old + df_cur - 3.0 * ((f_old - f_cur) / (a_old - a_cur))
            let d2Temp = (d1**2.0 - (df_old * df_cur)) |> sqrt
            let d2Sign = (a_cur - a_old) |> sign |> float
            let d2 = d2Sign * d2Temp
            let tRes = a_cur - (a_cur - a_old) * ((df_cur + d2 - d1) / (df_cur - df_old + 2.0 * d2))
            let (aMin, aMax) = ((min a_old a_cur), (max a_old a_cur))

            if System.Double.IsNaN(tRes) then None
            else if tRes < aMin then Some(aMin)
            else if tRes > aMax then Some(aMax)
            else Some(tRes)

    let initialMaxStep (v: Vector<float>) (d: Vector<float>) =
        if Vector.exists System.Double.IsNaN v || Vector.exists System.Double.IsInfinity v then None
        else if Vector.exists System.Double.IsNaN d || Vector.exists System.Double.IsInfinity d then None
        else
            let rec searchMaxStep (cMax: float, searchNum: int) =
                if searchNum >= searchMaxStepMax then None
                else
                    let y = f (v + cMax * d)
                    if System.Double.IsInfinity(y) || System.Double.IsNaN(y) then let nextxmax = cMax * searchMaxStepMult
                                                                                  searchMaxStep (nextxmax, searchNum + 1)
                    else Some(cMax)

            searchMaxStep(xMax, 0)
            
    member this.Search (v: Vector<float>) (d: Vector<float>) =
        let actualMaxStep = initialMaxStep v d
        match actualMaxStep with
        | None -> System.Double.NaN
        | Some(maxStep) -> 
            let actualInitStep = if maxStep < xInit then 0.5 * maxStep else xInit
            let phi = (fun (a: Vector<float>) -> f (v + a.[0] * d))
            let dphi = (fun a -> let res = Differentiation.Gradient(phi, a)
                                 res.[0])
            let phi_0, dphi_0 = DenseVector(1, 0.0) |> (fun x -> (phi x, dphi x))
            //let phi_i, dphi_i = DenseVector(1, actualInitStep) |> (fun x -> (phi x, dphi x))
            let phiMax, dphiMax = DenseVector(1, maxStep) |> (fun x -> (phi x, dphi x))

            let rec zoom ((a_l, phi_l, dphi_l) as low) ((a_h, phi_h, dphi_h) as high) =
                if a_h < dMin then rescueStepSize
                else
                    let int_a = (low, high) |> cubicInterpolation

                    match int_a with
                    | None -> rescueStepSize
                    | Some a_j
                        -> let zoomNotForSamePoint newlow newhigh =
                               if (newlow = low) && (newhigh = high) then a_j
                               else zoom newlow newhigh

                           let phi_j, dphi_j = [a_j] |> DenseVector.ofList |> (fun x -> (phi x, dphi x))
            
                           match phi_j with
                           | _ when phi_j > (phi_0 + c1 * a_j * dphi_0) || (phi_j >= phi_l) -> zoomNotForSamePoint (a_l, phi_l, dphi_l) (a_j, phi_j, dphi_j)
                           | _ -> match dphi_j with
                                  | _ when (abs dphi_j) <= -c2 * dphi_0 -> a_j
                                  | _ when dphi_j * (a_h - a_l) >= 0.0 -> zoomNotForSamePoint (a_j, phi_j, dphi_j) (a_l, phi_l, dphi_l)
                                  | _ -> zoomNotForSamePoint (a_j, phi_j, dphi_j) (a_h, phi_h, dphi_h)

            let rec search a_cur (a_old, phi_old, dphi_old as old) =
                let a_curVec = [a_cur] |> DenseVector.ofList
                let phi_cur, dphi_cur = (phi a_curVec, dphi a_curVec)

                match phi_cur with
                | _ when (phi_cur > phi_0 + c1 * a_cur * dphi_0) || (phi_cur >= phi_old && a_old <> actualInitStep)
                    -> zoom (a_old, phi_old, dphi_old) (a_cur, phi_cur, dphi_cur)
                | _ -> match dphi_cur with
                       | _ when (abs dphi_cur) <= -c2 * dphi_0 -> a_cur
                       | _ when dphi_cur >= 0.0 -> zoom (a_cur, phi_cur, dphi_cur) (a_old, phi_old, dphi_old)
                       | _ -> let a_new = cubicInterpolation ((a_cur, phi_cur, dphi_cur), (maxStep, phiMax, dphiMax))
                              match a_new with
                              | None -> rescueStepSize
                              | Some a_new ->
                                  if a_new >= maxStep then maxStep
                                  else if a_new <= a_cur then a_cur
                                  else (search a_new (a_cur, phi_cur, dphi_cur))

            match cubicInterpolation((0.0, phi_0, dphi_0), (maxStep, phiMax, dphiMax)) with
            | None -> rescueStepSize
            | Some x -> search actualInitStep (0.0, phi_0, dphi_0)

// "NelderMeadResult" class is for C#.
type NelderMeadResultFSharp = { Parameters: Vector<float>; FunctionValue: float; Converged: bool }

[<CompiledName "NelderMeadFSharp">]
type NelderMead (f:(Vector<float> -> float), iteration: int, tolerance: float) =
    let defaultIteration = 100
    let defaultTolerance = 1e-3

    let mutable m_Iteration = iteration
    let mutable m_Tolerance = tolerance
    let mutable m_ZDelta = 0.00025
    let mutable m_Delta = 0.05
    let mutable m_Rho = 1.0
    let mutable m_Chi = 2.0
    let mutable m_Psi = 0.5
    let mutable m_Sigma = 0.5

    member this.Iteration with get() = m_Iteration and set v = m_Iteration <- if v <= 0 then defaultIteration else v
    member this.Tolerance with get() = m_Tolerance and set v = m_Tolerance <- if v <= 0.0 then defaultTolerance else v
    member this.ZDelta with get() = m_ZDelta and set v = m_ZDelta <- v
    member this.Delta with get() = m_Delta and set v = m_Delta <- v
    member this.Rho with get() = m_Rho and set v = m_Rho <- v
    member this.Chi with get() = m_Chi and set v = m_Chi <- v
    member this.Psi with get() = m_Psi and set v = m_Psi <- v
    member this.Sigma with get() = m_Sigma and set v = m_Sigma <- v
                
    member this.Minimize(initval: Vector<float>) =
        // Creating Simplex
        let rec loopCS (vec: Vector<float>) acc simplex =
            if acc = vec.Count then simplex
            else vec |> Vector.mapi (fun i x -> if i = acc then match x with
                                                                | 0.0 -> this.ZDelta
                                                                | x -> (1.0 + this.Delta) * x
                                                else x) 
                     |> (fun s -> loopCS vec (acc+1) (s :: simplex))
        let ss = initval :: (loopCS initval 0 []) |> List.map (fun x -> (x, (f x)))
        
        let rec loopIT (oldSimplex: (Vector<float> * float) list) (newSimplex: (Vector<float> * float) list) (count: int) =
            let ascendingSimplex = List.sortBy (fun elem -> (snd elem)) newSimplex
            let descendingSimplex = List.rev ascendingSimplex

            let getVectorFX (s: (Vector<float> * float) list) =
                s |> List.unzip |> snd |> DenseVector.ofList

            let createNextNewSimplex (newv: (Vector<float> * float)) = 
                let otherThanLargest = List.tail descendingSimplex
                newv :: otherThanLargest

            // Checking Convergence
            let f_L2 =
                if count > 0 then
                    (getVectorFX oldSimplex) - (getVectorFX newSimplex) |> Vector.fold (fun a x -> a + x*x) 0.0 |> System.Math.Sqrt
                else
                    0.0
            if f_L2 < this.Tolerance && count > 0 then (List.head newSimplex) |> (fun x -> ((fst x), (snd x), true)) 
            else if count > this.Iteration && count > 0 then (List.head newSimplex) |> (fun x -> ((fst x), (snd x), false))          
            else
                let (best, worst, sndWorst) =
                    (List.head ascendingSimplex, List.head descendingSimplex, (descendingSimplex |> List.tail |> List.head))
                let centroid = List.reduce (+) (List.tail descendingSimplex |> List.unzip |> fst)
                               |> (fun x -> x / (float (descendingSimplex.Length) - 1.0))

                // Reflect
                let r = centroid + this.Rho * (centroid - (fst worst)) |> (fun x -> (x, f x))
            
                if (snd r) < (snd best) then       // Expand
                    let e = (fst r) + this.Chi * ((fst r) - centroid) |> (fun x -> (x, f x))
                    loopIT ascendingSimplex (createNextNewSimplex e) (count + 1)
                else if (snd sndWorst) <= (snd r) then     // Contract
                    let c =
                        if (snd r) < (snd worst) then
                            centroid + this.Psi * ((fst r) - centroid) |> (fun x -> (x, f x))    // Outside Contract
                        else
                            centroid + this.Psi * ((fst worst) - centroid) |> (fun x -> (x, f x))    // Inside Contract
                    if (snd c) <= (snd worst) then
                        loopIT newSimplex (createNextNewSimplex c) (count + 1)
                    else                // Shrink
                        let nNewSimplex = List.map (fun (v: (Vector<float> * float)) ->
                                                        let x = (fst best) + this.Sigma * ((fst v) - (fst best))
                                                        (x, f x)) newSimplex
                        loopIT newSimplex nNewSimplex (count + 1)
                else
                    loopIT newSimplex (createNextNewSimplex r) (count + 1)

        loopIT [] ss 0

    member this.ResultConvertToType (result:(Vector<float> * float * bool)) =
        let (parameters, fValue, converged) = result
        { Parameters = parameters; FunctionValue = fValue; Converged = converged } 

[<CompiledName "QuasiNewtonMethodStatusFSharp">]
type QuasiNewtonMethodStatus<'a> =
    Converged of 'a
    | NotConverged of 'a
    | FunctionValueInvalid
    | GradientInvalid
    | WeightMatrixInvalid
    | LineSearchFailure

[<CompiledName "QuasiNewtonSearchBuilderFSharp">]
type QuasiNewtonSearchBuilder() =
    member this.Bind(result: QuasiNewtonMethodStatus<'a>, rest: 'a -> QuasiNewtonMethodStatus<'b>) =
        match result with
        Converged(x)
        | NotConverged(x) -> rest x
        | FunctionValueInvalid -> FunctionValueInvalid
        | GradientInvalid -> GradientInvalid
        | WeightMatrixInvalid -> WeightMatrixInvalid
        | LineSearchFailure -> LineSearchFailure

    member this.Return(result) = result


type QuasiNewtonMethodResultFSharp = { Status: int; Parameters: Vector<double>; FunctionValue: System.Nullable<float>; InvertedWeightMatrix: Matrix<float> }


[<CompiledName "BFGSFSharp">]
type BFGS (f:(Vector<float> -> float), iteration: int, tolerance: float) =  
    let defaultIteration = 100
    let defaultTolerance = 1e-3

    let mutable m_Iteration = iteration
    let mutable m_Tolerance = tolerance
    let mutable m_InitialStepSize = 1.0
    let mutable m_MaxStepSize = 2.0
    let mutable m_DerivationMethod = (fun x -> Differentiation.Gradient(f, x))
    let mutable m_FirstTimeStepSizeMuiltiplier = 1.0

    let isInvalidFloat (x: float) =
        if System.Double.IsInfinity x || System.Double.IsNaN x then true
        else false
         
    let qnsearch = new QuasiNewtonSearchBuilder()

    member this.Iteration with get() = iteration and set v = m_Iteration <- if v <= 0 then defaultIteration else v
    member this.Tolerance with get() = tolerance and set v = m_Tolerance <- if v <= 0.0 then defaultTolerance else v
    member this.DerivationMethod with get() = m_DerivationMethod and set v = m_DerivationMethod <- v
    member this.InitialStepSize with get() = m_InitialStepSize and set v = m_InitialStepSize <- v
    member this.MaxStepSize with get() = m_MaxStepSize and set v = m_MaxStepSize <- v
    member this.FirstTimeStepSizeMultiplier with get() = m_FirstTimeStepSizeMuiltiplier and set v = m_FirstTimeStepSizeMuiltiplier <- v

    member private this.differentiation x =
        let tRes = this.DerivationMethod x
        if Vector.exists isInvalidFloat tRes then GradientInvalid
        else NotConverged(tRes)

    member private this.BFGSWeightMatrix (w: Matrix<float>) (xds: Vector<float>) (dds: Vector<float>) =
        let tRes = w + ((1.0 / (dds * xds)) * Vector.OuterProduct(dds, dds)) - ((1.0 / (xds * w * xds)) * Vector.OuterProduct((w * xds), (w * xds)))
        if Matrix.exists isInvalidFloat tRes then WeightMatrixInvalid
        else NotConverged(tRes)

    member private this.lineSearch r g =
        let ls = LineSearch(f, this.InitialStepSize, this.MaxStepSize)
        let tRes = ls.Search r g
        if isInvalidFloat tRes then LineSearchFailure
        else NotConverged(tRes)

    member this.FSResultToCSResult(result: QuasiNewtonMethodStatus<Vector<float> * float * Matrix<float>>) =
        match result with
            Converged((x, f, w)) -> { Status = 0; Parameters = x; FunctionValue = new System.Nullable<float>(f); InvertedWeightMatrix = w }
            | NotConverged((x, f, w)) -> { Status = 1; Parameters = x; FunctionValue = new System.Nullable<float>(f); InvertedWeightMatrix = w }
            | FunctionValueInvalid -> { Status = 2; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | GradientInvalid -> { Status = 3; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | WeightMatrixInvalid -> { Status = 4; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | LineSearchFailure -> { Status = 5; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }

    member this.Minimize(initVal: Vector<float>) =
        let rec search (w: Matrix<float>) (r: Vector<float>) (g: Vector<float>) count =
            qnsearch {
                let wi = w.Inverse()

                if g.Norm(2.0) < this.Tolerance && count > 0 then return Converged(r, (f r), wi)
                else if (count >= this.Iteration) then return NotConverged(r, (f r), wi)
                else let! step = this.lineSearch r ((-1.0) * (wi * g))
                     let newR = r - step * (wi * g)
                     let! newG = this.differentiation newR
                     let! newW = this.BFGSWeightMatrix w (newR - r) (newG - g)
                     return search newW newR newG (count + 1)
            }

        qnsearch {
            let initX = initVal.Clone()
            let! initD = this.differentiation initX
            let w = let mult = if this.FirstTimeStepSizeMultiplier = 1.0 then 1.0
                               else (1.0 / this.FirstTimeStepSizeMultiplier) * initD.Norm(2.0)
                    mult * DenseMatrix.Identity(initX.Count)

            return search w initX initD 0
        }
