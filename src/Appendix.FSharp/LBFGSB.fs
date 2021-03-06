﻿// Math.NET Numerics Appendix Library License (MIT/X11)
// ===================================================
// 
// Copyright (c) 2016 Fukui Shogo
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

namespace MathNet.Numerics.Appendix.Optimization

open System
//open System.Collections.Generic
open MathNet.Numerics
open MathNet.Numerics.Differentiation
//open MathNet.Numerics.Interpolation
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double


type internal GeneralizedCauchyPointResult<'a, 'b> =
    | Normal of 'a
    | Corner of 'b

module internal GeneralizedCauchyPoint =
    let private BreakPointsOnDirection (values: float array) (gradients: float array) (bounds: (float * float) array) : float array =
        if values.Length <> gradients.Length || values.Length <> bounds.Length
            then failwith "Lengths of \"values\", \"gradients\", and \"bounds\" must agree."
        else if Array.exists2 (fun v bound -> v < (fst bound) || v > (snd bound)) values bounds
            then failwith "One or more values are out of the bound."
        else Array.map3 (fun v g bound ->
                            if g < 0.0 then (v - (snd bound)) / g
                            else if g > 0.0 then (v - (fst bound)) / g
                            else Double.PositiveInfinity) values gradients bounds

    let private BreakPointsSort (breakPoints: float array) : (int * float) array =
        breakPoints |> Array.zip [| 0..(breakPoints.Length-1) |]
                    |> Array.sortBy (fun (i, point) -> point)
        
    let SearchWithBFGSMatrix (values: float array) (gradients: float array) (bounds: (float * float) array) (bfgsMatrix: Matrix<float>) =
        let breaksOnDirecton = BreakPointsOnDirection values gradients bounds 
        
        let initialSortedBreaks = BreakPointsSort breaksOnDirecton |> Array.filter (fun (i, bp) -> bp > 0.0)
        let initialDirections = Array.mapi (fun i bp -> if bp > 0.0 then -gradients.[i] else 0.0) breaksOnDirecton
        let initialStepD =
            let vecDirs = initialDirections |> DenseVector.ofArray
            (vecDirs * vecDirs) / (vecDirs * bfgsMatrix * vecDirs)
        let vecGraidents = gradients |> DenseVector.ofArray

        let rec search xs sortedBreaks dirs stepD oldStep =
            if (Array.isEmpty sortedBreaks) then Corner(xs)
            else if (stepD < (snd sortedBreaks.[0] - oldStep)) then
                let step = oldStep + (max stepD 0.0)
                let cp = Array.map2 (fun x d -> x + step * d) xs dirs
                let freeIndice = Array.unzip sortedBreaks |> fst |> Array.sort
                Normal(cp, freeIndice)
            else
                let (newXs, newDirs) =
                    Array.mapi2 (fun i x dir -> if i <> (fst sortedBreaks.[0]) then (x, dir)
                                                else
                                                    if dirs.[i] > 0.0 then (snd bounds.[i], 0.0)
                                                    else (fst bounds.[i], 0.0) ) xs dirs |> Array.unzip
                let newSortedBreaks = Array.tail sortedBreaks
                let vecZs = Array.map2 (fun x1 x0 -> x1 - x0) newXs xs |> DenseVector.OfArray
                let newStepD =
                    let vecDirs = newDirs |> DenseVector.ofArray
                    let numer = (-vecDirs * vecDirs) + (vecDirs * bfgsMatrix * vecZs)
                    let denom = (vecDirs * bfgsMatrix * vecDirs)
                    if (abs numer) < Double.Epsilon then 0.0 else - numer / denom
                search newXs newSortedBreaks newDirs newStepD (snd sortedBreaks.[0])

        search values initialSortedBreaks initialDirections initialStepD 0.0

    let Search (values: float array) (gradients: float array) (bounds: (float * float) array) (theta: float) (matW: Matrix<float>) (matM: Matrix<float>) =
        let breaksOnDirection = BreakPointsOnDirection values gradients bounds

        let initialSortedBreaks = BreakPointsSort breaksOnDirection |> Array.filter (fun (i, bp) -> bp > 0.0)
        let initialDirections = Array.mapi (fun i bp -> if bp > 0.0 then -gradients.[i] else 0.0) breaksOnDirection
        let vecDirs = initialDirections |> DenseVector.OfArray
        let initialP = matW.Transpose() * vecDirs
        let (initialDerivative1, initialDerivative2) =
            let df = - vecDirs * vecDirs
            let d2f = - theta * df - (initialP * matM * initialP)
            (df, d2f)
        let initialStepD = - initialDerivative1 / initialDerivative2
        let initialC = Array.init matW.ColumnCount (fun _ -> 0.0) |> DenseVector.ofArray
        let vecGraidents = gradients |> DenseVector.ofArray

        let rec search xs sortedBreaks dirs stepD (p: Vector<float>) (c: Vector<float>) derivative1 derivative2 oldStep =
            if Array.isEmpty sortedBreaks then Corner(xs)
            else if stepD < (snd sortedBreaks.[0] - oldStep) then
                let step = oldStep + (max stepD 0.0)
                let gcp = Array.map2 (fun x d -> x + step * d) xs dirs
                let freeIndice = Array.unzip sortedBreaks |> fst |> Array.sort
                Normal(gcp, freeIndice, (c + stepD * p))
            else
                let b = sortedBreaks.[0] |> fst
                let (newXs, newDirs) =
                    Array.mapi2 (fun i x dir -> if i <> b then (x, dir)
                                                else
                                                    if dirs.[i] > 0.0 then (snd bounds.[i], 0.0)
                                                    else (fst bounds.[i], 0.0) ) xs dirs |> Array.unzip
                let newSortedBreaks = Array.tail sortedBreaks
                let newC = c + stepD * p
                let (zb, gb, wb) = ((newXs.[b] - (snd sortedBreaks.[0])), gradients.[0], matW.Row(b))
                let newDF =
                    derivative1 + stepD * derivative1 + gb**2.0 + theta * gb * zb - gb * (wb * matM * newC)
                let newD2F = 
                    derivative2 - theta * (gb**2.0) - 2.0 * gb * (wb * matM * p) - (gb**2.0) * (wb * matM * wb)
                let newP = p + gb * wb
                let newStepD = - newDF / newD2F
                
                search newXs newSortedBreaks newDirs newStepD newP newC newDF newD2F (snd sortedBreaks.[0])
        
        search values initialSortedBreaks initialDirections initialStepD initialP initialC initialDerivative1 initialDerivative2 0.0 


module internal DirectPrimalMethod =
    let private PartialIdentityMatrix (indice: int array) (valueNum: int) =
        Array.map (fun i -> Array.init valueNum (fun j -> if i = j then 1.0 else 0.0)) indice
            |> DenseMatrix.ofColumnArrays

    // Backtracking multiplier calculation        
    let private SimpleBacktracking (values: float array) (valuesGCP: float array) (freeIndice: int array) (bounds: (float * float) array) (directions: Vector<float>) =
        let filter (fi: int) =
            let dir = directions.At(fi)
            if (abs dir) < Double.Epsilon then 0.0
            else 
                let (lbound, ubound) = bounds.[fi]
                let xGCP = valuesGCP.[fi]
                max ((ubound - xGCP) / dir) ((lbound - xGCP) / dir) |> min 1.0
        Array.init freeIndice.Length filter |> Array.max |> (*) directions

    let SearchWithBFGSMatrix (values: float array) (valuesGCP: float array) (freeIndice: int array) (gradients: float array) (bounds: (float * float) array) (bfgsMatrix: Matrix<float>) =
        // Matrix $Z_k$
        let matZ = PartialIdentityMatrix freeIndice values.Length
        let matZT = matZ.Transpose()
        // Matrix $A_k$
        let matA =
            let constrainedIndice = Array.except freeIndice [|0..(values.Length - 1)|]
            if constrainedIndice.Length = 0 then None
            else PartialIdentityMatrix constrainedIndice values.Length |> Some
        
        let vecValues = DenseVector.ofArray values
        let vecValuesGCP = DenseVector.ofArray valuesGCP
        let vecGradients = DenseVector.OfArray gradients

        let reducedGradients = 
            matZT * (vecGradients + bfgsMatrix * (vecValuesGCP - vecValues))
        let reducedBFGSMatrix =
            matZT * bfgsMatrix * matZ

        let dirsCandidate = - (reducedBFGSMatrix.Inverse()) * reducedGradients
        
        let backtrackedDirs = SimpleBacktracking values valuesGCP freeIndice bounds dirsCandidate

        let rec summarize (restFreeIndice: int array) (index: int) (bdirCounter: int) (xs: float array) (ds: float array) =
            if restFreeIndice.Length = 0 then (xs, ds)
            else
                if index = restFreeIndice.[0] then
                    //let (newx, newd) = (xs.[index] + backtrackedDirs.[bdirCounter], backtrackedDirs.[bdirCounter])
                    let (newx, newd) = (xs.[index], backtrackedDirs.[bdirCounter])
                    do xs.[index] <- newx
                    do ds.[index] <- newd
                    summarize (Array.tail restFreeIndice) (index+1) (bdirCounter+1) xs ds
                else
                    summarize (restFreeIndice) (index + 1) bdirCounter xs ds

        summarize freeIndice 0 0 (Array.copy valuesGCP) (Array.init gradients.Length (fun _ -> 0.0)) 

    let Search (values: float array) (valuesGCP: float array) (freeIndice: int array) (gradients: float array) (bounds: (float * float) array)
        (theta: float) (matW: Matrix<float>) (matM: Matrix<float>) (vecC: Vector<float>) =
        let matWT = matW.Transpose()
        // Matrix $Z_k$
        let matZ = PartialIdentityMatrix freeIndice values.Length
        let matZT = matZ.Transpose()
        // Matrix $A_k$
        let matA =
            let constrainedIndice = Array.except freeIndice [|0..(values.Length - 1)|]
            if constrainedIndice.Length = 0 then None
            else PartialIdentityMatrix constrainedIndice values.Length |> Some

        let vecValues = DenseVector.ofArray values
        let vecValuesGCP = DenseVector.ofArray valuesGCP
        let vecGradients = DenseVector.OfArray gradients

        let rc = matZT * (vecGradients + theta * (vecValuesGCP - vecValues) - (matW * matM * vecC))
        let matI = DenseMatrix.identity matM.RowCount
        let matN =
            let temp = (1.0 / theta) * matWT * matZ * matZT * matW
            matI - (matM * temp)
        let nu =
            let temp = matM * (matWT * matZ * rc)
            matN.Inverse() * temp

        let dirsCandidate = (1.0 / theta) * rc + (1.0 / (theta * theta)) * matZT * matW * nu |> (*) (-1.0)

        let backtrackedDirs = SimpleBacktracking values valuesGCP freeIndice bounds dirsCandidate

        let rec summarize (restFreeIndice: int array) (index: int) (bdirCounter: int) (xs: float array) (ds: float array) =
            if restFreeIndice.Length = 0 then (xs, ds)
            else
                if index = restFreeIndice.[0] then
                    //let (newx, newd) = (xs.[index] + backtrackedDirs.[bdirCounter], backtrackedDirs.[bdirCounter])
                    let (newx, newd) = (xs.[index], backtrackedDirs.[bdirCounter])
                    do xs.[index] <- newx
                    do ds.[index] <- newd
                    summarize (Array.tail restFreeIndice) (index+1) (bdirCounter + 1) xs ds
                else
                    summarize (restFreeIndice) (index + 1) bdirCounter xs ds

        summarize freeIndice 0 0 (Array.copy valuesGCP) (Array.init gradients.Length (fun _ -> 0.0)) 


[<CompiledName "SubOptimizationTypeFSharp">]
type SubOptimizationType =
     | DirectPrimal


[<CompiledName "CorrectionFSharp">]
type Correction(size: int) =
    let mutable m_Memory = Array.empty<float array>

    member this.Memory with get() = m_Memory
    member this.MaxSize with get() = size
    member this.CurrentSize with get() = m_Memory.Length
    member this.Item with get(i:int) = this.Memory.[i]

    member this.Add(xs: float array) =
        if this.CurrentSize < this.MaxSize then m_Memory <- Array.append m_Memory [|xs|]
        else m_Memory <- let oldmemory = m_Memory.[1..(this.MaxSize - 1)]
                         Array.append oldmemory [|xs|] 
    member this.GetMatrix() =
        DenseMatrix.ofColumnArrays this.Memory


[<CompiledName "FuncWithBoundsFSharp">]
type FuncWithBounds(f: System.Func<float array, float>, bounds: (float * float) array) =
    member this.VariableNumber with get() = bounds.Length
    member this.Bounds with get() = bounds
    member this.Func with get() = f

    member this.ValidateBounds(xs: float array) =
        Array.forall2 (fun x bound -> x >= (fst bound) && x <= (snd bound)) xs bounds

    member this.EvaluateRaw(xs: float array) = this.Func.Invoke(xs)
    member this.Evaluate(xs: float array) =
        if xs.Length <> bounds.Length then invalidArg "xs" "The lengths of x and bounds are not same."
        else if this.ValidateBounds(xs) = false then invalidArg "xs" "One or more x are outside of the bounds."
        else this.EvaluateRaw(xs)


[<CompiledName "LBFGSBStatusFSharp">]
type LBFGSBStatus<'a> =
    | InProcess of 'a
    | Converged of 'a
    | NotConverged of 'a
    | FunctionValueInvalid
    | GradientInvalid
    | BFGSMatrixInvalid
    | LineSearchFailure
    | InverseBFGSMatrixInvalid
    | NoGeneralizedCauchyPoint
    | BlockOfBFGSMatrixInvalid
    | BlockOfInverseBFGSMatrixInvalid
    | CorrectionHistoryInvalid
    | ConvergedAtCorner of 'a


[<CompiledName "LBFGSBResultFSharp">]
type LBFGSBResult = { Status: int; Values: float[]; FunctionValue: System.Nullable<float>; InverseBFGSMatrix: float[,] }


type internal LBFGSBSearchBuilder() =
    member this.Bind(result: LBFGSBStatus<'a>, rest: 'a -> LBFGSBStatus<'b>) : LBFGSBStatus<'b> =
        match result with
        | InProcess(x) -> rest x
        | Converged(x) -> failwith "\"Converged\" case cannot be passed to bind function."
        | NotConverged(x) -> failwith "\"NotConverged\" case cannot be passed to bind function."
        | FunctionValueInvalid -> FunctionValueInvalid
        | GradientInvalid -> GradientInvalid
        | BFGSMatrixInvalid -> BFGSMatrixInvalid
        | LineSearchFailure -> LineSearchFailure
        | InverseBFGSMatrixInvalid -> InverseBFGSMatrixInvalid
        | NoGeneralizedCauchyPoint -> NoGeneralizedCauchyPoint
        | BlockOfBFGSMatrixInvalid -> BlockOfBFGSMatrixInvalid
        | BlockOfInverseBFGSMatrixInvalid -> BlockOfInverseBFGSMatrixInvalid
        | CorrectionHistoryInvalid -> CorrectionHistoryInvalid
        | ConvergedAtCorner(x)-> failwith "\"ConvergedAtCorner\" case cannot be passed to bind function."

    member this.ReturnFrom(result: LBFGSBStatus<'a>) : LBFGSBStatus<'a> =
        result

    member this.Return(result: 'a) : LBFGSBStatus<'a> =
        NotConverged(result)


[<CompiledName "LBFGSBFSharp">]
type LBFGSB(f: FuncWithBounds, iteration: int, tolerance: float, approxdimension: int) = 
    let (defaultCentralDerivation, defaultForwardDerivation, defaultBackwardDerivation)
        = (new NumericalDerivative(3,1), new NumericalDerivative(3,0), new NumericalDerivative(3,2))
    let mutable m_Bounds = f.Bounds
    let mutable m_DerivationMethod =
        new System.Func<float[], float[]>(fun xs -> Array.init f.VariableNumber
                                                        (fun i -> let (lbound, ubound) = m_Bounds.[i]
                                                                  if xs.[i] - lbound < Double.Epsilon then defaultForwardDerivation.EvaluatePartialDerivative(f.Func, xs, i, 1)
                                                                  else if ubound - xs.[i] < Double.Epsilon then defaultBackwardDerivation.EvaluatePartialDerivative(f.Func, xs, i, 1)
                                                                  else defaultCentralDerivation.EvaluatePartialDerivative(f.Func, xs, i, 1) ) )
    let mutable m_LineSearchMethod = new Appendix.Optimization.LineSearch(f.Func, 1.0, 10.0, 10)
    let mutable m_SubOptimizationMethod = DirectPrimal
    let mutable m_ApproxBFGSMatrixDimension = approxdimension
    // Matrix $S_k$
    let mutable m_ValueCorrections = new Correction(m_ApproxBFGSMatrixDimension)
    // Matrix $Y_k$
    let mutable m_GradientCorrections = new Correction(m_ApproxBFGSMatrixDimension)
    
    let mutable m_WriteTrace = NoTrace

    /// Fields for optimization.
    member this.BoundedFunction = f
    member this.Iteration = iteration
    member this.Tolerance = tolerance
    member this.DerivationMethod with get() = m_DerivationMethod and set v = m_DerivationMethod <- v
    member this.SubOptimizationMethod with get() = m_SubOptimizationMethod and set v = m_SubOptimizationMethod <- v
    // Variable $m$
    member this.ApproxBFGSMatrixDimension with get () = m_ApproxBFGSMatrixDimension and set (v) = m_ApproxBFGSMatrixDimension <- v
    // Bounds of variables
    member this.ValueBounds with get() = m_Bounds
    // Line search method
    member this.LineSearchMethod with get() = m_LineSearchMethod and set (v) = m_LineSearchMethod <- v

    new(f: System.Func<float array, float>, bounds: (float * float) array, iteration: int, tolerance: float, approxdimension: int) =
        let boundedFunc = new FuncWithBounds(f, bounds)
        LBFGSB(boundedFunc, iteration, tolerance, approxdimension)
    new(f: FuncWithBounds, iteration: int, tolerance: float) =
        LBFGSB(f, iteration, tolerance, 3)
    new(f: System.Func<float array, float>, bounds: (float * float) array, iteration: int, tolerance: float) =
        LBFGSB(f, bounds, iteration, tolerance, 3)

    member this.TraceToStdOut = 
        do m_WriteTrace <- StdOut

    member this.TraceToTextWriter (writer: System.IO.TextWriter) =
        do m_WriteTrace <- TextWriter(writer)

    member this.TraceNone =
        do m_WriteTrace <- NoTrace

    member this.FSResultToCSResult(result: LBFGSBStatus<float[] * float * float[,]>) : LBFGSBResult =
        match result with
        | InProcess(_) -> failwith "\"InProcess\" case never become the result of minimization."
        | Converged((x, f, w)) -> { Status = 1; Values = x; FunctionValue = new System.Nullable<float>(f); InverseBFGSMatrix = w }
        | NotConverged((x, f, w)) -> { Status = 2; Values = x; FunctionValue = new System.Nullable<float>(f); InverseBFGSMatrix = w }
        | FunctionValueInvalid -> { Status = 3; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | GradientInvalid -> { Status = 4; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | BFGSMatrixInvalid -> { Status = 5; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | LineSearchFailure -> { Status = 6; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | InverseBFGSMatrixInvalid -> { Status = 7; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | NoGeneralizedCauchyPoint -> { Status = 8; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | BlockOfBFGSMatrixInvalid -> { Status = 9; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | BlockOfInverseBFGSMatrixInvalid -> { Status = 10; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | CorrectionHistoryInvalid -> { Status = 11; Values = null; FunctionValue = new System.Nullable<float>(); InverseBFGSMatrix = null }
        | ConvergedAtCorner((x, f, w)) -> { Status = 12; Values = x; FunctionValue = new System.Nullable<float>(f); InverseBFGSMatrix = w }

    /// Fields and methods for limited memory BFGS matrix computation. 
    // Multiplier $\theta$ calculator.
    member this.Theta() =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0
            then 1.0
        else
            let (newestS, newestY) = (m_ValueCorrections.[m_ValueCorrections.CurrentSize - 1], m_GradientCorrections.[m_GradientCorrections.CurrentSize - 1])
            (AppendixFunctions.InnerProductOfTwoArrays newestY newestY) / (AppendixFunctions.InnerProductOfTwoArrays newestY newestS)        

    // Matrix $L_k$ calculator.
    member this.LMatrix() =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0
            then failwith "There exist no correnctions about value/gradient."
        else
            DenseMatrix.init m_ValueCorrections.CurrentSize m_GradientCorrections.CurrentSize
                (fun i j -> if i > j then AppendixFunctions.InnerProductOfTwoArrays m_ValueCorrections.[i] m_GradientCorrections.[j] else 0.0)

    // Matrix $R_k$ calculator.
    member this.RMatrix() =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0
            then failwith "There exist no correnctions about value/gradient."
        else
            DenseMatrix.init m_ValueCorrections.CurrentSize m_GradientCorrections.CurrentSize
                (fun i j -> if i <= j then AppendixFunctions.InnerProductOfTwoArrays m_ValueCorrections.[i] m_GradientCorrections.[j] else 0.0)        

    // Matrix $D_k$ calculator.
    member this.DMatrix() =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0
            then failwith "There exist no correnctions about value/gradient."
        else
            DenseMatrix.initDiag m_ValueCorrections.CurrentSize m_GradientCorrections.CurrentSize
                (fun i -> AppendixFunctions.InnerProductOfTwoArrays m_ValueCorrections.[i] m_GradientCorrections.[i])

    // Matrix $W_k$ calculator.
    member this.WMatrix(theta: float) =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0
            then failwith "There exist no correnctions about value/gradient."
        else
            let (matY, thetaMatS) = (m_GradientCorrections.GetMatrix(), theta * m_ValueCorrections.GetMatrix())
            matY.Append(thetaMatS)

    // Matrix $\bar{W}_k$ calculator.
    member this.BarWMatrix(theta: float) =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0
            then failwith "There exist no correnctions about value/gradient."
        else
            let (rThetaMatY, MatS) = ((1.0 / theta) * m_GradientCorrections.GetMatrix(), m_ValueCorrections.GetMatrix())
            rThetaMatY.Append(MatS)

    // An matrix($M_k$) in BFGS matrix approximation.
    member this.MMatrix(negMatD: Matrix<float>, matL: Matrix<float>, matS: Matrix<float>, theta: float) : LBFGSBStatus<Matrix<float>> =
        if negMatD.RowCount <> matL.ColumnCount || matL.RowCount <> matS.ColumnCount || negMatD.ColumnCount <> matL.ColumnCount || matL.RowCount <> matS.ColumnCount
            then failwith "One or more sets of block have different rows/columns each other."
        else
            let matLT = matL.Transpose()
            let thetaSS = matS.Transpose() * matS * theta

            let inverseNegMatD = DenseMatrix.initDiag negMatD.RowCount negMatD.ColumnCount (fun i -> 1.0 / negMatD.At(i, i))
            let temp1 = (thetaSS - (matL * inverseNegMatD * matLT)).Inverse()
            if Matrix.exists AppendixFunctions.IsInvalidFloat temp1 then BlockOfBFGSMatrixInvalid
            else
                AppendixFunctions.CreateMatrixFrom2x2Blocks
                    (inverseNegMatD + inverseNegMatD * matLT * temp1 * matL * inverseNegMatD) (-inverseNegMatD * matLT * temp1)
                    (-temp1 * matL * inverseNegMatD)                                          (temp1)
                    |> InProcess
    
    member this.BarMMatrix(matD: Matrix<float>, matR: Matrix<float>, matY: Matrix<float>, theta: float) : LBFGSBStatus<Matrix<float>> =
        if matD.RowCount <> matR.RowCount || matY.ColumnCount <> matD.ColumnCount
            then failwith "One or more sets of block have different rows/columns each other."
        else
            let invMatR = matR.Inverse()
            if Matrix.exists AppendixFunctions.IsInvalidFloat invMatR then BlockOfInverseBFGSMatrixInvalid
            else
                let negInvMatR = (-1.0) * invMatR
                let temp1 = matD + (1.0 / theta) * (matY.Transpose() * matY)
                let block22 = invMatR.Transpose() * temp1 * invMatR
                AppendixFunctions.CreateMatrixFrom2x2Blocks
                    (DenseMatrix.init negInvMatR.RowCount negInvMatR.RowCount (fun _ _ -> 0.0))  (negInvMatR)
                    (negInvMatR.Transpose())                                                     (block22)
                    |> InProcess

    // BFGS matrix($B_k$) computation
    member this.BFGSMatrix() : LBFGSBStatus<Matrix<float>> =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0 then DenseMatrix.identity f.VariableNumber |> InProcess
        else
            // Variable $\theta$
            let theta = this.Theta()
            // Matrix $I$
            let matI = DenseMatrix.identity f.VariableNumber
            // Matrix $L_k$ and $D_k$
            let (matL, matD) = (this.LMatrix(), this.DMatrix())
            // Matrix $W_k$
            let matW = this.WMatrix(theta)
            // Matrix $M$ calculated from other matrices.
            let optMatM = this.MMatrix(-matD, matL, m_ValueCorrections.GetMatrix(), theta)
            
            match optMatM with
            | InProcess(matM) -> theta * matI - matW * matM * matW.Transpose() |> InProcess
            | BlockOfBFGSMatrixInvalid -> BFGSMatrixInvalid
            | _ -> failwith "Unexpected error in matrix $M$ calculation."
    
    member this.InverseBFGSMatrix() : LBFGSBStatus<Matrix<float>> =
        if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0 then DenseMatrix.identity f.VariableNumber |> InProcess
        else
            // Variable $\theta$
            let theta = this.Theta()
            // Matrix $I$
            let matI = DenseMatrix.identity f.VariableNumber
            // Matrix $R_k$ and $D_k$
            let (matR, matD) = (this.RMatrix(), this.DMatrix())
            // Matrix $\bar{W}_k$
            let matBarW = this.BarWMatrix(theta)
            // Matrix $\bar{M}$ calculated from other matrices.
            let optMatBarM = this.BarMMatrix(matD, matR, m_GradientCorrections.GetMatrix(), theta)

            match optMatBarM with
            | InProcess(matBarM) -> (1.0 / theta) * matI + matBarW * matBarM * matBarW.Transpose() |> InProcess
            | BlockOfInverseBFGSMatrixInvalid -> InverseBFGSMatrixInvalid
            | _ -> failwith "Unexpected error in matrix $M$ calculation."

    member private this.Trace (writeline: string -> unit) (values: float array) (gradients: float array) (sw: System.Diagnostics.Stopwatch) =
        let (vecValues, vecGradients, inverseBFGSMatrix) =
            (DenseVector.ofArray values), (DenseVector.ofArray gradients), (this.InverseBFGSMatrix())
        do writeline("---- Tracing Log of BFGS Optimization ----")
        do writeline("Elapsed Time:")
        do writeline(sw.Elapsed.ToString())
        do writeline("Estimated Parameters:")
        writeline(vecValues.ToVectorString())
        do writeline("Gradients:")
        writeline(vecGradients.ToVectorString())
        do writeline("Inverted BFGS Matrix:")
        match inverseBFGSMatrix with
        | InProcess(mat) -> writeline(mat.ToMatrixString(mat.RowCount, mat.ColumnCount, null))
        | _ -> writeline("Inverted BFGS matrix could not be calculated.")
        do writeline("")

    member private this.LineSearch(subOptimizedValues: float array, subOptimizedDirections: float array, count: int) : LBFGSBStatus<(float array) * (float array)> =  
        try
            let step =
                if (Array.forall (fun d -> (abs d) < Double.Epsilon) subOptimizedDirections ) || count = 0 then 0.0
                else
                    this.LineSearchMethod.Search subOptimizedValues subOptimizedDirections
            let newValues = Array.map3 (fun x d bounds -> let newValue = x + step * d
                                                          let (lbound, ubound) = bounds
                                                          if newValue < lbound then lbound
                                                          else if newValue > ubound then ubound
                                                          else newValue) subOptimizedValues subOptimizedDirections this.ValueBounds
            let newGradients = this.DerivationMethod.Invoke(newValues)
            InProcess(newValues, newGradients)     
        with
            | ex -> LineSearchFailure

    member private this.CorrectionUpdate(currentValues: float array, currentGradients: float array, newValues: float array, newGradients: float array) : unit =
        let (newValuedDiff, newGradientsDiff) = ((AppendixFunctions.PairwiseDifferenceOfTwoArrays newValues currentValues),
                                                 (AppendixFunctions.PairwiseDifferenceOfTwoArrays newGradients currentGradients))
        if (AppendixFunctions.InnerProductOfTwoArrays newValuedDiff newGradientsDiff > Double.Epsilon * (AppendixFunctions.L2NormOfArray newGradientsDiff)) then
            do m_ValueCorrections.Add(newValuedDiff)
            do m_GradientCorrections.Add(newGradientsDiff)

    // Minimization
    member this.Minimize(initVal: float array) =
        let lbfgsbSearch = new LBFGSBSearchBuilder()
        let sw = new System.Diagnostics.Stopwatch()
        do sw.Start()

        let initGradients = this.DerivationMethod.Invoke(initVal)

        let rec search (currentValues: float array) (currentGradients: float array) (currentFreeVariableIndice: int array) (count: int) =
            let freeVarDirections = Array.map (fun fi -> currentGradients.[fi]) currentFreeVariableIndice

            match m_WriteTrace with
            | StdOut -> do this.Trace (System.Console.WriteLine) currentValues currentGradients sw
            | TextWriter(writer) -> do this.Trace writer.WriteLine currentValues currentGradients sw
                                    do writer.Flush()
            | NoTrace -> ()

            lbfgsbSearch {
                if freeVarDirections.Length > 0 && this.Tolerance > AppendixFunctions.InfinityNormOfArray freeVarDirections then
                    do sw.Stop()
                    let! inverseBFGSMatrix = this.InverseBFGSMatrix()
                    return! Converged(currentValues, this.BoundedFunction.EvaluateRaw(currentValues), inverseBFGSMatrix.ToArray())
                else if count >= this.Iteration then
                    do sw.Stop()
                    let! inverseBFGSMatrix = this.InverseBFGSMatrix()
                    return! NotConverged(currentValues, this.BoundedFunction.EvaluateRaw(currentValues), inverseBFGSMatrix.ToArray())
                else
                    match this.SubOptimizationMethod with
                    | DirectPrimal ->
                        if m_ValueCorrections.CurrentSize = 0 && m_GradientCorrections.CurrentSize = 0 then
                            let! bfgsMatrix = this.BFGSMatrix()
                            match GeneralizedCauchyPoint.SearchWithBFGSMatrix currentValues currentGradients this.ValueBounds bfgsMatrix with
                            | Corner(xs) -> if Array.forall2 (fun oldX newX -> abs (oldX - newX) < Double.Epsilon) currentValues xs then
                                                let! inverseBFGSMatrix = this.InverseBFGSMatrix()
                                                return! ConvergedAtCorner(xs, this.BoundedFunction.EvaluateRaw(xs), inverseBFGSMatrix.ToArray())
                                            else
                                                let newGradients = this.DerivationMethod.Invoke(xs)
                                                return! search xs newGradients Array.empty (count + 1)
                            | Normal(gcp, freeVarIndice) ->
                                let (subXs, subDirs) = DirectPrimalMethod.SearchWithBFGSMatrix currentValues gcp freeVarIndice currentGradients this.ValueBounds bfgsMatrix
                                let! (newValues, newGradients) = this.LineSearch(subXs, subDirs, count)
                                do this.CorrectionUpdate(currentValues, currentGradients, newValues, newGradients)
                                return! search newValues newGradients freeVarIndice (count + 1)
                        else if m_ValueCorrections.CurrentSize = 0 || m_GradientCorrections.CurrentSize = 0 then
                            return! CorrectionHistoryInvalid
                        else
                            let (theta, matD, matL, matS) = (this.Theta(), this.DMatrix(), this.LMatrix(), m_ValueCorrections.GetMatrix())
                            let matW = this.WMatrix(theta)
                            let! matM = this.MMatrix(-matD, matL, matS, theta)
                            match GeneralizedCauchyPoint.Search currentValues currentGradients this.ValueBounds theta matW matM with
                            | Corner(xs) -> if Array.forall2 (fun oldX newX -> abs (oldX - newX) < Double.Epsilon) currentValues xs then
                                                let! inverseBFGSMatrix = this.InverseBFGSMatrix()
                                                return! ConvergedAtCorner(xs, this.BoundedFunction.EvaluateRaw(xs), inverseBFGSMatrix.ToArray())
                                            else
                                                let newGradients = this.DerivationMethod.Invoke(xs)
                                                return! search xs newGradients Array.empty (count + 1)
                            | Normal(gcp, freeVarIndice, vecC) ->
                                let (subXs, subDirs) = DirectPrimalMethod.Search currentValues gcp freeVarIndice currentGradients this.ValueBounds theta matW matM vecC
                                let! (newValues, newGradients) = this.LineSearch(subXs, subDirs, count)
                                do this.CorrectionUpdate(currentValues, currentGradients, newValues, newGradients)
                                return! search newValues newGradients freeVarIndice (count + 1)
            }

        search initVal initGradients (Array.init initVal.Length (fun i -> i)) 0

