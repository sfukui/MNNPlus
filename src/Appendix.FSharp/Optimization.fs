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

type internal LineSearchInfo = { A: float; Phi: float; DPhi: float; }

type internal searchInterpolationResult = | Plain of LineSearchInfo
                                          | Conditioned of LineSearchInfo
                                          | Failed

type LineSearch ( f: System.Func<float[], float>, xInit: float, xMax: float, trialMax: int) =
    let defaultNumericalDerivative = new Differentiation.NumericalDerivative()
    let defaultDerivationMethod = new System.Func<System.Func<float, float>, float, float>( (fun f a -> defaultNumericalDerivative.EvaluateDerivative(f, a, 1)) )

    let mutable (m_c1, m_c2) = (10.0**(-4.0), 0.9)  // General settings for Wolfe conditions.
    let mutable m_DMin = 1e-16                      // When a difference of two values is less than dMin, the difference is considered as zero.   
    let mutable m_MaxStepSearchMult = 0.5           // Multiplier for temporal max step size in valid max step search. 
    let mutable m_MaxStepSearchTrial = 100          // Max trial number of valid max step search.
    let mutable m_InterpolationSearchTrial = 10     // Max trial number of valid interpolated step size.
    let mutable m_NextStepMult = 0.1                // Multiplier to get next step size candidate.
    let mutable m_DerivationMethod = defaultDerivationMethod // Numerical differentiation function.

    let createFuncs (v: float[]) (d: float[]) =
        let vad a = Array.map2 (fun vi di -> vi + a * di) v d
        let phi = fun (a: float) -> let ov = a |> vad
                                    in f.Invoke(ov)
        let phiFunc = new System.Func<float, float>(phi)
        let dphi = (fun (a: float) -> m_DerivationMethod.Invoke(phiFunc, a))
        (phi, dphi)

    let cubicInterpolation (oldInfo: LineSearchInfo) (newInfo: LineSearchInfo) =
        if (abs (newInfo.A - oldInfo.A)) < m_DMin then Some(newInfo.A)
        else
            let d1 = oldInfo.DPhi + newInfo.DPhi - 3.0 * ((oldInfo.Phi - newInfo.Phi) / (oldInfo.A - newInfo.A))
            let d2Temp = (d1 * d1 - oldInfo.DPhi * newInfo.DPhi)
            if d2Temp < 0.0 then None
            else
                let d2Sign = (newInfo.A - oldInfo.A) |> sign |> float
                let d2 = d2Temp |> sqrt |> (*) d2Sign
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

    let searchDownValidStep (xCurrentMax: float) (v: float[]) (d: float[]) =
        if Array.exists System.Double.IsNaN v || Array.exists System.Double.IsInfinity v then Failed
        else if Array.exists System.Double.IsNaN d || Array.exists System.Double.IsInfinity d then Failed
        else
            let rec searchStep (cMax: float, searchCount: int) =
                if searchCount >= m_MaxStepSearchTrial then Failed
                else
                    let (phi, dphi) = let funcs = createFuncs v d
                                      in (fst funcs, snd funcs)
                    let phiCMax = phi cMax
                    let dphiCMax = dphi cMax
                    if System.Double.IsInfinity(phiCMax) || System.Double.IsNaN(phiCMax) ||
                       System.Double.IsInfinity(dphiCMax) || System.Double.IsNaN(dphiCMax) then
                        let nextxmax = cMax * m_MaxStepSearchMult
                        searchStep (nextxmax, searchCount + 1)
                    else
                        let res = { A = cMax; Phi = phiCMax; DPhi = dphiCMax }
                        if searchCount = 0 then Plain(res)
                        else Conditioned(res)

            searchStep(xCurrentMax, 0)

    let searchValidCubicInterpolation (v: float[]) (d: float[]) =
        let rec searchValidPoint (oldInfo: LineSearchInfo) (newInfo: LineSearchInfo) (count: int) =
            let tempInterpolation = cubicInterpolation oldInfo newInfo
            if (count >= m_InterpolationSearchTrial) then None
            else match tempInterpolation with
                 | Some x -> match searchDownValidStep x v d with
                             | Plain info -> Some info
                             | Conditioned info -> if oldInfo.A <= newInfo.A then searchValidPoint oldInfo info (count + 1)
                                                   else searchValidPoint info newInfo (count + 1)
                             | Failed -> None
                 | None -> None

        fun (oldInfo: LineSearchInfo) (newInfo: LineSearchInfo) ->
            searchValidPoint oldInfo newInfo 0

    member this.C1 with get() = m_c1 and set v = m_c1 <- if v <= 0.0 then 10.0**(-4.0) else v
    member this.C2 with get() = m_c2 and set v = m_c2 <- if v <= 0.0 then 0.9 else v
    member this.MinimumDifference with get() = m_DMin and set v = m_DMin <- if v <= 0.0 then 1e-16 else v
    member this.MaxStepSearchMultiplier with get() = m_MaxStepSearchMult and set v = m_MaxStepSearchMult <- if v <= 0.0 || v >= 1.0 then 0.5 else v    
    member this.MaxStepSearchTrial with get() = m_MaxStepSearchTrial and set v = m_MaxStepSearchTrial <- if v <= 0 then 100 else v
    member this.InterpolationSearchTrial with get() = m_InterpolationSearchTrial and set v = m_InterpolationSearchTrial <- if v <= 0 then 10 else v
    member this.NextStepMultiplier with get() = m_NextStepMult and set v = m_NextStepMult <- if v <= 0.0 || v >= 1.0 then 0.1 else v
    member this.NumericalDerivation with get() = m_DerivationMethod and set v = m_DerivationMethod <- v
         
    member this.Search (v: float[]) (d: float[]) =
        let actualMaxStep = searchDownValidStep xMax v d
        let validCubicInterpolationSearch = searchValidCubicInterpolation v d
        match actualMaxStep with
        | Failed -> System.Double.NaN
        | Plain infoMax | Conditioned infoMax -> 
            let actualInitialStep = let tempInit = if xInit >= infoMax.A then infoMax.A * 0.5 else xInit
                                    searchDownValidStep tempInit v d
            let maxStep = infoMax.A
            let (phi, dphi) = let funcs = createFuncs v d
                              in (fst funcs, snd funcs)
            let info0  = let (phi0, dphi0) = 0.0 |> (fun x -> (phi x, dphi x))
                         in { LineSearchInfo.A = 0.0; LineSearchInfo.Phi = phi0; LineSearchInfo.DPhi = dphi0 }

            let searchNextStepRange = m_NextStepMult * (infoMax.A - info0.A)

            let rec zoom (low: LineSearchInfo) (high: LineSearchInfo) (trialCount: int) =
                if (abs (low.A - high.A)) < m_DMin then low.A
                else
                    let interpolatedInfo = validCubicInterpolationSearch low high
 
                    match interpolatedInfo with
                    | None -> failwith "Cannot interpolate current minimum and maximum step values."
                    | Some pick ->
                        if trialCount > trialMax then pick.A
                        else    
                            match pick.Phi with
                            | _ when pick.Phi > (info0.Phi + m_c1 * pick.A * info0.DPhi) || (pick.Phi >= low.Phi)
                                     -> zoom low pick (trialCount + 1)
                            | _ -> match pick.DPhi with
                                   | _ when (abs pick.DPhi) <= -m_c2 * info0.DPhi -> pick.A
                                   | _ when pick.DPhi * (high.A - low.A) >= 0.0 -> zoom pick low (trialCount + 1)
                                   | _ -> zoom pick high (trialCount + 1)

            let rec search (currentInfo: LineSearchInfo) (oldInfo: LineSearchInfo) (trialCount: int) =
                if trialCount > trialMax then currentInfo.A
                else
                    match currentInfo.Phi with
                    | _ when (currentInfo.Phi > info0.Phi + m_c1 * currentInfo.A * info0.DPhi) || (currentInfo.Phi >= oldInfo.Phi && trialCount <> 0)
                        -> zoom oldInfo currentInfo trialCount
                    | _ -> match currentInfo.DPhi with
                           | _ when (abs currentInfo.DPhi) <= -m_c2 * info0.DPhi -> currentInfo.A
                           | _ when currentInfo.DPhi >= 0.0 -> zoom currentInfo oldInfo trialCount
                           | _ -> let newInfoInterpolation = searchDownValidStep (currentInfo.A + searchNextStepRange) v d
                                  match newInfoInterpolation with
                                  | Failed -> failwith "Cannot find next step value in line search."
                                  | Plain newInfo | Conditioned newInfo ->
                                      if newInfo.A >= infoMax.A || newInfo.A <= currentInfo.A then failwith "Next step value went over the limits."
                                      else search newInfo currentInfo (trialCount + 1)

            match actualInitialStep with
            | Plain info | Conditioned info -> search info info0 0
            | Failed -> failwith "Cannot find initial step value in line search."


// "NelderMeadResult" class is for C#.
type NelderMeadResultFSharp = { Parameters: float[]; FunctionValue: float; Converged: bool }

[<CompiledName "NelderMeadFSharp">]
type NelderMead (f:(float[] -> float), iteration: int, tolerance: float) =
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
                
    member this.Minimize(initval: float[]) =
        // Creating Simplex
        let rec loopCS (vec: float[]) acc simplex =
            if acc = vec.Length then simplex
            else vec |> Array.mapi (fun i x -> if i = acc then match x with
                                                               | 0.0 -> this.ZDelta
                                                               | x -> (1.0 + this.Delta) * x
                                                else x) 
                     |> (fun s -> loopCS vec (acc+1) (s :: simplex))
        let ss = initval :: (loopCS initval 0 []) |> List.map (fun x -> ((DenseVector.OfArray x) :> Vector<float>, (f x)))
        
        let rec loopIT (oldSimplex: (Vector<float> * float) list) (newSimplex: (Vector<float> * float) list) (count: int) =
            let ascendingSimplex = List.sortBy (fun elem -> (snd elem)) newSimplex
            let descendingSimplex = List.rev ascendingSimplex

            let getVectorFX (s: (Vector<float> * float) list) =
                s |> List.unzip |> snd |> DenseVector.ofList

            let createNextNewSimplex (newv: (Vector<float> * float)) = 
                let otherThanLargest = List.tail descendingSimplex
                newv :: otherThanLargest

            let getXYSet (x: Vector<float>) = let xArr = x.ToArray()
                                              in (x, f xArr)

            // Checking Convergence
            let f_L2 =
                if count > 0 then
                    (getVectorFX oldSimplex) - (getVectorFX newSimplex) |> Vector.fold (fun a x -> a + x*x) 0.0 |> System.Math.Sqrt
                else
                    0.0
            if f_L2 < this.Tolerance && count > 0 then (List.head newSimplex) |> (fun x -> ((fst x), (snd x), true)) 
            else if count > this.Iteration && count > 0 then (List.head newSimplex) |> (fun x -> ((fst x), (snd x), false))
            else if System.Double.IsNaN(f_L2) || System.Double.IsInfinity(f_L2) then
                let invalidres = List.init initval.Length (fun i -> System.Double.NaN) |> DenseVector.ofList |> Vector.map (fun x -> x)
                (invalidres, f_L2, false)
            else
                let (best, worst, sndWorst) =
                    (List.head ascendingSimplex, List.head descendingSimplex, (descendingSimplex |> List.tail |> List.head))
                let centroid = List.reduce (+) (List.tail descendingSimplex |> List.unzip |> fst)
                               |> (fun x -> x / (float (descendingSimplex.Length) - 1.0))

                // Reflect
                let r = centroid + this.Rho * (centroid - (fst worst)) |> getXYSet
            
                if (snd r) < (snd best) then       // Expand
                    let e = (fst r) + this.Chi * ((fst r) - centroid) |> getXYSet
                    loopIT ascendingSimplex (createNextNewSimplex e) (count + 1)
                else if (snd sndWorst) <= (snd r) then     // Contract
                    let c =
                        if (snd r) < (snd worst) then
                            centroid + this.Psi * ((fst r) - centroid) |> getXYSet    // Outside Contract
                        else
                            centroid + this.Psi * ((fst worst) - centroid) |> getXYSet    // Inside Contract
                    if (snd c) <= (snd worst) then
                        loopIT newSimplex (createNextNewSimplex c) (count + 1)
                    else                // Shrink
                        let nNewSimplex = List.map (fun (v: (Vector<float> * float)) ->
                                                        let x = (fst best) + this.Sigma * ((fst v) - (fst best))
                                                        let xArr = x.ToArray()
                                                        (x, f xArr)) newSimplex
                        loopIT newSimplex nNewSimplex (count + 1)
                else
                    loopIT newSimplex (createNextNewSimplex r) (count + 1)

        let getResults (rawRes: Vector<float> * float * bool) =
            let (ps, fVal, conv) = rawRes
            in (ps.ToArray(), fVal, conv)
        
        (loopIT [] ss 0) |> getResults

    member this.ResultConvertToType (result:(float[] * float * bool)) =
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


type QuasiNewtonMethodResultFSharp = { Status: int; Parameters: float[]; FunctionValue: System.Nullable<float>; InvertedWeightMatrix: float[,] }

[<CompiledName "TraceOutputFSharp">]
type TraceOutput<'a> =
    TextWriter of 'a
    | StdOut
    | NoTrace 

type internal StepSizeStatus =
    | ValidStep of float
    | InvalidStep
    | NotYet

[<CompiledName "BFGSFSharp">]
type BFGS (f: System.Func<float[], float>, iteration: int, tolerance: float) = 
    let defaultIteration = 100
    let defaultTolerance = 1e-1
    let defaultDerivation = new NumericalJacobian()
    let defaultLineSearch = new LineSearch(f, 1.0, 10.0, 10)

    let mutable m_Iteration = iteration
    let mutable m_Tolerance = tolerance

    let mutable m_DerivationMethod = new System.Func<float[], float[]>(fun x -> defaultDerivation.Evaluate(f, x))
    let mutable m_LineSearch = defaultLineSearch

    let mutable m_FirstTimeStepSizeMuiltiplier = 1.0

    let mutable m_LatestStepSize = NotYet
    let mutable m_LatestXVector = None
    let mutable m_LatestGradientVector = None
    let mutable m_LatestInvertedWeightMatrix = None

    let mutable m_WriteTrace = NoTrace

    let isInvalidFloat (x: float) =
        if System.Double.IsInfinity x || System.Double.IsNaN x then true
        else false

    let insideDerivationMethod (xVec: Vector<float>) =
        let xArr = xVec.ToArray()
        m_DerivationMethod.Invoke(xArr) |> DenseVector.ofArray
         
    let qnsearch = new QuasiNewtonSearchBuilder()

    member this.Iteration with get() = m_Iteration and set v = m_Iteration <- if v <= 0 then defaultIteration else v
    member this.Tolerance with get() = m_Tolerance and set v = m_Tolerance <- if v <= 0.0 then defaultTolerance else v
    member this.DerivationMethod with get() = m_DerivationMethod and set v = m_DerivationMethod <- v
    member this.LineSearch with get() = m_LineSearch and set v = m_LineSearch <- v
    member this.FirstTimeStepSizeMultiplier with get() = m_FirstTimeStepSizeMuiltiplier and set v = m_FirstTimeStepSizeMuiltiplier <- v

    member this.LatestStepSize with get() = match m_LatestStepSize with
                                            | ValidStep(x) -> Some(x)
                                            | _ -> None 
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
        let tRes = x |> insideDerivationMethod
        if Vector.exists isInvalidFloat tRes then
             do m_LatestGradientVector <- None
             GradientInvalid
        else
             do m_LatestGradientVector <- Some(tRes)
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
            do m_LatestInvertedWeightMatrix <- None
            InvertedWeightMatrixInvalid
        else
            do m_LatestInvertedWeightMatrix <- Some(tRes)
            NotConverged(tRes)

    member private this.lineSearch (r: Vector<float>) (g: Vector<float>) =
        let (rArr, gArr) = (r.ToArray(), g.ToArray())
        let tRes = m_LineSearch.Search rArr gArr
        if isInvalidFloat tRes then
            do m_LatestStepSize <- InvalidStep
            LineSearchFailure
        else do m_LatestStepSize <- ValidStep(tRes)
             NotConverged(tRes)

    member this.GetResults(result: QuasiNewtonMethodStatus<Vector<float> * float * Matrix<float>>) =
        match result with
            Converged(x, f, w) -> Converged(x.ToArray(), f, w.ToArray())
            | NotConverged(x, f, w) -> NotConverged(x.ToArray(), f, w.ToArray())
            | FunctionValueInvalid -> FunctionValueInvalid
            | GradientInvalid -> GradientInvalid
            | WeightMatrixInvalid -> WeightMatrixInvalid
            | LineSearchFailure -> LineSearchFailure
            | InvertedWeightMatrixInvalid -> InvertedWeightMatrixInvalid

    member this.FSResultToCSResult(result: QuasiNewtonMethodStatus<float[] * float * float[,]>) =
        match result with
            Converged((x, f, w)) -> { Status = 0; Parameters = x; FunctionValue = new System.Nullable<float>(f); InvertedWeightMatrix = w }
            | NotConverged((x, f, w)) -> { Status = 1; Parameters = x; FunctionValue = new System.Nullable<float>(f); InvertedWeightMatrix = w }
            | FunctionValueInvalid -> { Status = 2; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | GradientInvalid -> { Status = 3; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | WeightMatrixInvalid -> { Status = 4; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | LineSearchFailure -> { Status = 5; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }
            | InvertedWeightMatrixInvalid -> { Status = 6; Parameters = null; FunctionValue = new System.Nullable<float>(); InvertedWeightMatrix = null }

    member private this.Trace (writeline: string -> unit) (sw: System.Diagnostics.Stopwatch) =
        do writeline("---- Tracing Log of BFGS Optimization ----")
        do writeline("Elapsed Time:")
        do writeline(sw.Elapsed.ToString())
        do writeline("Estimated Parameters:")
        match m_LatestXVector with
        | Some(x : Vector<float>) -> writeline(x.ToVectorString())
        | None -> writeline("NaN")
        do writeline("Gradients:")
        match m_LatestGradientVector with
        | Some(x : Vector<float>) -> writeline(x.ToVectorString())
        | None -> writeline("NaN")
        do writeline("Inverted Weight Matrix:")
        match m_LatestInvertedWeightMatrix with
        | Some(x : Matrix<float>) -> writeline(x.ToMatrixString(x.RowCount, x.ColumnCount, null))
        | None -> writeline("NaN")
        do writeline("Step Size:")
        match m_LatestStepSize with
        | ValidStep(x : float) -> writeline(x.ToString())
        | InvalidStep -> writeline("NaN")
        | NotYet -> writeline("Not computed yet.")
        do writeline("")

    member this.Minimize(initVal: float[]) =
        let sw = new System.Diagnostics.Stopwatch()
        do sw.Start()

        let rec search (winv: Matrix<float>) (r: Vector<float>) (g: Vector<float>) count = 
            match m_WriteTrace with
            | StdOut -> do this.Trace (System.Console.WriteLine) sw
            | TextWriter(writer) -> do this.Trace writer.WriteLine sw
                                    do writer.Flush()
            | NoTrace -> ()

            qnsearch {
                let cur_fval = f.Invoke(r.ToArray())

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
            let initX = initVal |> DenseVector.ofArray
            let! initD = this.differentiation initX
            let initW = let mult = if this.FirstTimeStepSizeMultiplier = 1.0 then 1.0
                                   else (1.0 / this.FirstTimeStepSizeMultiplier) * initD.Norm(2.0)
                        mult * DenseMatrix.identity(initX.Count)
            let initWInv = let wDiagInv = initW.Diagonal() |> Vector.map (fun x -> 1.0 / x)
                           in DenseMatrix.initDiag wDiagInv.Count wDiagInv.Count (fun i -> wDiagInv.[i])

            do m_LatestXVector <- Some(initX)
            do m_LatestInvertedWeightMatrix <- Some(initWInv)
            return (search initWInv initX initD 0) |> this.GetResults
        }
