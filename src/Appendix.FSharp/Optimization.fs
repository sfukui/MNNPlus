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
open MathNet.Numerics.Interpolation
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double

type LineSearchInfo = { A: float; Phi: float; DPhi: float; }

[<CompiledName "LineSearchFSharp">]
type LineSearch ( f: (Vector<float> -> float), xInit: float, xMax: float, trialMax: int) =
    let (c1, c2) = (10.0**(-4.0), 0.9)
    let dMin = 1e-16
    let rescueStepSize = 1.0

    let searchMaxStepMult = 0.5
    let searchMaxStepMax = 100

    let cubicInterpolation (oldInfo: LineSearchInfo) (newInfo: LineSearchInfo) =
        if (abs (newInfo.A - oldInfo.A)) < dMin then Some(newInfo.A)
        else
            let d1 = oldInfo.DPhi + newInfo.DPhi - 3.0 * ((oldInfo.Phi - newInfo.Phi) / (oldInfo.A - newInfo.A))
            let d2Temp = (d1 * d1 - oldInfo.DPhi * newInfo.DPhi) |> sqrt
            let d2Sign = (newInfo.A - oldInfo.A) |> sign |> float
            let d2 = d2Sign * d2Temp
            let tRes = newInfo.A - (newInfo.A - oldInfo.A) * (newInfo.DPhi + d2 - d1) / (newInfo.DPhi - oldInfo.DPhi + 2.0 * d2)

            if System.Double.IsNaN(tRes) then None
            else
                if oldInfo.A <= newInfo.A then
                    if tRes < oldInfo.A then Some(oldInfo.A)
                    else if tRes > newInfo.A then Some(newInfo.A)
                    else Some(tRes)
                else
                    if tRes < newInfo.A then Some(newInfo.A)
                    else if tRes > oldInfo.A then Some(oldInfo.A)
                    else Some(tRes)

    let initialMaxStep (v: Vector<float>) (d: Vector<float>) =
        if Vector.exists System.Double.IsNaN v || Vector.exists System.Double.IsInfinity v then None
        else if Vector.exists System.Double.IsNaN d || Vector.exists System.Double.IsInfinity d then None
        else
            let rec searchMaxStep (cMax: float, searchNum: int) =
                if searchNum >= searchMaxStepMax then None
                else
                    let y = f (v + cMax * d)
                    if System.Double.IsInfinity(y) || System.Double.IsNaN(y) then
                        let nextxmax = cMax * searchMaxStepMult
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
            let dphi = (fun a -> let res = Differentiation.Derivative(phi, a)
                                 res.[0])
            let info0  = let (phi0, dphi0) = DenseVector.Create(1, 0.0) |> (fun x -> (phi x, dphi x))
                         in { LineSearchInfo.A = 0.0; LineSearchInfo.Phi = phi0; LineSearchInfo.DPhi = dphi0 }
            let infoMax = let (phiMax, dphiMax) = DenseVector.Create(1, maxStep) |> (fun x -> (phi x, dphi x))
                          in { LineSearchInfo.A = maxStep; LineSearchInfo.Phi = phiMax; LineSearchInfo.DPhi = dphiMax }

            let rec zoom (low: LineSearchInfo) (high: LineSearchInfo) (trialCount: int) =
                if (abs (low.A - high.A)) < dMin then low.A
                else
                    let interpolatedA = cubicInterpolation low high
 
                    match interpolatedA with
                    | None -> rescueStepSize
                    | Some aJ ->
                        if trialCount > trialMax then aJ
                        else
                            let pick = let (phiJ, dphiJ) = [aJ] |> DenseVector.ofList |> (fun x -> (phi x, dphi x))
                                       in { LineSearchInfo.A = aJ; LineSearchInfo.Phi = phiJ; LineSearchInfo.DPhi = dphiJ }
            
                            match pick.Phi with
                            | _ when pick.Phi > (info0.Phi + c1 * pick.A * info0.DPhi) || (pick.Phi >= low.Phi)
                                     -> zoom low pick (trialCount + 1)
                            | _ -> match pick.DPhi with
                                   | _ when (abs pick.DPhi) <= -c2 * info0.DPhi -> pick.A
                                   | _ when pick.DPhi * (high.A - low.A) >= 0.0 -> zoom pick low (trialCount + 1)
                                   | _ -> zoom pick high (trialCount + 1)

            let rec search aCur (old: LineSearchInfo) (trialCount: int) =
                if trialCount > trialMax then aCur
                else
                    let current = let (phiCur, dphiCur) = [aCur] |> DenseVector.ofList |> (fun a -> (phi a, dphi a))
                                  in { LineSearchInfo.A = aCur; LineSearchInfo.Phi = phiCur; LineSearchInfo.DPhi = dphiCur }

                    match current.Phi with
                    | _ when (current.Phi > info0.Phi + c1 * current.A * info0.DPhi) || (current.Phi >= old.Phi && old.A <> actualInitStep)
                        -> zoom old current trialCount
                    | _ -> match current.DPhi with
                           | _ when (abs current.DPhi) <= -c2 * info0.DPhi -> current.A
                           | _ when current.DPhi >= 0.0 -> zoom current old trialCount
                           | _ -> let aNewOption = cubicInterpolation current infoMax
                                  match aNewOption with
                                  | None -> rescueStepSize
                                  | Some aNewVal ->
                                      if aNewVal >= infoMax.A then infoMax.A
                                      else if aNewVal <= current.A then current.A
                                      else search aNewVal current (trialCount + 1)

            match cubicInterpolation info0 infoMax with
            | None -> rescueStepSize
            | Some x -> search actualInitStep info0 0

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
            else if System.Double.IsNaN(f_L2) || System.Double.IsInfinity(f_L2) then
                let invalidres = List.init initval.Count (fun i -> System.Double.NaN) |> DenseVector.ofList |> Vector.map (fun x -> x)
                (invalidres, f_L2, false)
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
    | InvertedWeightMatrixInvalid

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
        | InvertedWeightMatrixInvalid -> InvertedWeightMatrixInvalid

    member this.Return(result) = result


type QuasiNewtonMethodResultFSharp = { Status: int; Parameters: Vector<double>; FunctionValue: System.Nullable<float>; InvertedWeightMatrix: Matrix<float> }

[<CompiledName "TraceOutputFSharp">]
type TraceOutput<'a> =
    TextWriter of 'a
    | StdOut
    | NoTrace 

[<CompiledName "BFGSFSharp">]
type BFGS (f:(Vector<float> -> float), iteration: int, tolerance: float) =  
    let defaultIteration = 100
    let defaultTolerance = 1e-1

    let mutable m_Iteration = iteration
    let mutable m_Tolerance = tolerance

    let mutable m_InitialStepSize = 1.0
    let mutable m_MaxStepSize = 10.0

    let defaultLineSearchMaxTrial = 10
    let mutable m_LineSearchMaxTrial = defaultLineSearchMaxTrial

    let mutable m_DerivationMethod = (fun x -> Differentiation.Derivative(f, x, true, false))
    let mutable m_FirstTimeStepSizeMuiltiplier = 1.0

    let mutable m_LatestStepSize = None
    let mutable m_LatestXVector = None
    let mutable m_LatestGradientVector = None
    let mutable m_LatestInvertedWeightMatrix = None

    let mutable m_WriteTrace = NoTrace

    let isInvalidFloat (x: float) =
        if System.Double.IsInfinity x || System.Double.IsNaN x then true
        else false
         
    let qnsearch = new QuasiNewtonSearchBuilder()

    member this.Iteration with get() = m_Iteration and set v = m_Iteration <- if v <= 0 then defaultIteration else v
    member this.Tolerance with get() = m_Tolerance and set v = m_Tolerance <- if v <= 0.0 then defaultTolerance else v
    member this.DerivationMethod with get() = m_DerivationMethod and set v = m_DerivationMethod <- v
    member this.InitialStepSize with get() = m_InitialStepSize and set v = m_InitialStepSize <- v
    member this.MaxStepSize with get() = m_MaxStepSize and set v = m_MaxStepSize <- v
    member this.LineSearchMaxTrial with get() = m_LineSearchMaxTrial and set v = m_LineSearchMaxTrial <- if v <= 0 then defaultLineSearchMaxTrial else v
    member this.FirstTimeStepSizeMultiplier with get() = m_FirstTimeStepSizeMuiltiplier and set v = m_FirstTimeStepSizeMuiltiplier <- v

    member this.LatestStepSize with get() = m_LatestStepSize
    member this.LatestXVector with get() = m_LatestXVector
    member this.LatestGradientVector with get() = m_LatestGradientVector
    member this.LatestInvertedWeightMatrix with get() = m_LatestInvertedWeightMatrix

    member this.TraceToStdOut = 
        do m_WriteTrace <- StdOut

    member this.TraceToTextWriter (writer: System.IO.TextWriter) =
        do m_WriteTrace <- TextWriter(writer)

    member this.TraceNone =
        do m_WriteTrace <- NoTrace

    member private this.differentiation x =
        let tRes = this.DerivationMethod x
        if Vector.exists isInvalidFloat tRes then GradientInvalid
        else do m_LatestGradientVector <- Some(tRes)
             NotConverged(tRes)

    member private this.BFGSWeightMatrix (w: Matrix<float>) (xds: Vector<float>) (dds: Vector<float>) =
        let tRes = w + ((1.0 / (dds * xds)) * Vector.OuterProduct(dds, dds)) - ((1.0 / (xds * w * xds)) * Vector.OuterProduct((w * xds), (w * xds)))
        if Matrix.exists isInvalidFloat tRes then
            WeightMatrixInvalid
        else NotConverged(tRes)

    member private this.InverseOfBFGSWeightMatrix (wInv: Matrix<float>) (xds: Vector<float>) (dds: Vector<float>) =
        let (rho, ident) = (1.0 / (dds * xds), (DenseMatrix.identity xds.Count))
        let tRes = let (matL, matR) = (ident - rho * Vector.OuterProduct(xds, dds), ident - rho * Vector.OuterProduct(dds, xds))
                   in matL * wInv * matR + rho * Vector.OuterProduct(xds, xds)
        if Matrix.exists isInvalidFloat tRes then
            InvertedWeightMatrixInvalid
        else do m_LatestInvertedWeightMatrix <- Some(tRes)
             NotConverged(tRes)

    member private this.lineSearch r g =
        let ls = LineSearch(f, this.InitialStepSize, this.MaxStepSize, m_LineSearchMaxTrial)
        let tRes = ls.Search r g
        if isInvalidFloat tRes then LineSearchFailure
        else do m_LatestStepSize <- Some(tRes)
             NotConverged(tRes)

    member this.FSResultToCSResult(result: QuasiNewtonMethodStatus<Vector<float> * float * Matrix<float>>) =
        match result with
            Converged((x, f, w)) -> { Status = 0; Parameters = x; FunctionValue = new System.Nullable<float>(f); InvertedWeightMatrix = w }
            | NotConverged((x, f, w)) -> { Status = 1; Parameters = x; FunctionValue = new System.Nullable<float>(f); InvertedWeightMatrix = w }
            | FunctionValueInvalid -> { Status = 2; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | GradientInvalid -> { Status = 3; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | WeightMatrixInvalid -> { Status = 4; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | LineSearchFailure -> { Status = 5; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | InvertedWeightMatrixInvalid -> { Status = 6; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }

    member private this.Trace (writeline: string -> unit) (sw: System.Diagnostics.Stopwatch) =
        do writeline("---- Tracing Log of BFGS Oprimization ----")
        do writeline("Elapsed Time:")
        do writeline(sw.Elapsed.ToString())
        do writeline("Estimated Parameters:")
        match this.LatestXVector with
        | Some(x : Vector<float>) -> writeline(x.ToVectorString())
        | None -> writeline("NaN")
        do writeline("Gradients:")
        match this.LatestGradientVector with
        | Some(x : Vector<float>) -> writeline(x.ToVectorString())
        | None -> writeline("NaN")
        do writeline("Inverted Weight Matrix:")
        match this.LatestInvertedWeightMatrix with
        | Some(x : Matrix<float>) -> writeline(x.ToMatrixString(x.RowCount, x.ColumnCount, null))
        | None -> writeline("NaN")
        do writeline("Step Size:")
        match this.LatestStepSize with
        | Some(x : float) -> writeline(x.ToString())
        | None -> writeline("NaN")

    member this.Minimize(initVal: Vector<float>) =
        let sw = new System.Diagnostics.Stopwatch();
        do sw.Start()

        let rec search (winv: Matrix<float>) (r: Vector<float>) (g: Vector<float>) count = 
            match m_WriteTrace with
            | StdOut -> do this.Trace (System.Console.WriteLine) sw
            | TextWriter(writer) -> do this.Trace writer.WriteLine sw
                                    do writer.Flush()
            | NoTrace -> ()

            qnsearch {
                let cur_fval = f r

                if System.Double.IsNaN(cur_fval) || System.Double.IsInfinity(cur_fval) then return FunctionValueInvalid
                else if g.Norm(2.0) < this.Tolerance
                then do sw.Stop()
                     return Converged(r, cur_fval, winv)
                else if (count >= this.Iteration) then return NotConverged(r, cur_fval, winv)
                else let! step = this.lineSearch r ((-1.0) * (winv * g))
                     let newR = r - step * (winv * g)
                     do m_LatestXVector <- Some(newR)
                     let! newG = this.differentiation newR
                     let! newWInv = this.InverseOfBFGSWeightMatrix winv (newR - r) (newG - g)
                     return search newWInv newR newG (count + 1)
            }

        qnsearch {
            let initX = initVal.Clone()
            let! initD = this.differentiation initX
            let initW = let mult = if this.FirstTimeStepSizeMultiplier = 1.0 then 1.0
                                   else (1.0 / this.FirstTimeStepSizeMultiplier) * initD.Norm(2.0)
                        mult * DenseMatrix.identity(initX.Count)
            let initWInv = let wDiagInv = initW.Diagonal() |> Vector.map (fun x -> 1.0 / x)
                           in DenseMatrix.initDiag wDiagInv.Count wDiagInv.Count (fun i -> wDiagInv.[i])

            return search initWInv initX initD 0
        }
