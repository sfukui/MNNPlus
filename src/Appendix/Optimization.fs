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

type LineSearch (f: (Vector<float> -> float), xinit, xmax) =
    let c1, c2 = 10.0**(-4.0), 0.9
    let dmin = 1e-16
    let rescuestepsize = 1.0
    let searchmaxstepmult = 0.5
    let searchmaxstepmax = 100

    let cubicInterpolation((a_old, f_old, df_old), (a_cur, f_cur, df_cur)) =
        if abs( (a_old - a_cur) / a_cur) < dmin then Some(a_cur)
        else
            let d1 = df_old + df_cur - 3.0 * ((f_old - f_cur) / (a_old - a_cur))
            let d2temp = (d1**2.0 - (df_old * df_cur)) |> sqrt
            let d2sign = (a_cur - a_old) |> sign |> float
            let d2 = d2sign * d2temp
            let tres = a_cur - (a_cur - a_old) * ((df_cur + d2 - d1) / (df_cur - df_old + 2.0 * d2))
            let (amin, amax) = ((min a_old a_cur), (max a_old a_cur))

            if System.Double.IsNaN(tres) then None
            else if tres < amin then Some(amin)
            else if tres > amax then Some(amax)
            else Some(tres)

    let initialMaxStep (v: Vector<float>) (d: Vector<float>) =
        if Vector.exists System.Double.IsNaN v || Vector.exists System.Double.IsInfinity v then None
        else if Vector.exists System.Double.IsNaN d || Vector.exists System.Double.IsInfinity d then None
        else
            let rec searchMaxStep (cmax: float, searchnum: int) =
                if searchnum >= searchmaxstepmax then None
                else
                    let y = f (v + cmax * d)
                    if System.Double.IsInfinity(y) || System.Double.IsNaN(y) then let nextxmax = cmax * searchmaxstepmult
                                                                                  searchMaxStep (nextxmax, searchnum + 1)
                    else Some(cmax)

            searchMaxStep(xmax, 0)
            
    member x.Search (v: Vector<float>) (d: Vector<float>) =
        let actualMaxStep = initialMaxStep v d
        match actualMaxStep with
        | None -> System.Double.NaN
        | Some(maxstep) -> 
            let actualInitStep = if maxstep < xinit then 0.5 * maxstep else xinit
            let phi = (fun (a: Vector<float>) -> f (v + a.[0] * d))
            let dphi = (fun a -> let res = Differentiation.Gradient(phi, a)
                                 res.[0])
            let phi0, dphi0 = DenseVector(1, 0.0) |> (fun x -> (phi x, dphi x))
            let phii, dphii = DenseVector(1, actualInitStep) |> (fun x -> (phi x, dphi x))
            let phimax, dphimax = DenseVector(1, maxstep) |> (fun x -> (phi x, dphi x))

            let rec zoom ((a_l, phi_l, dphi_l) as low) ((a_h, phi_h, dphi_h) as high) =
                if a_h < dmin then rescuestepsize
                else
                    let int_a = (low, high) |> cubicInterpolation

                    match int_a with
                    | None -> rescuestepsize
                    | Some a_j
                        -> let zoomnotforsamepoint newlow newhigh =
                               if (newlow = low) && (newhigh = high) then a_j
                               else zoom newlow newhigh

                           let phi_j, dphi_j = [a_j] |> DenseVector.ofList |> (fun x -> (phi x, dphi x))
            
                           match phi_j with
                           | _ when phi_j > (phi0 + c1 * a_j * dphi0) || (phi_j >= phi_l) -> zoomnotforsamepoint (a_l, phi_l, dphi_l) (a_j, phi_j, dphi_j)
                           | _ -> match dphi_j with
                                  | _ when (abs dphi_j) <= -c2 * dphi0 -> a_j
                                  | _ when dphi_j * (a_h - a_l) >= 0.0 -> zoomnotforsamepoint (a_j, phi_j, dphi_j) (a_l, phi_l, dphi_l)
                                  | _ -> zoomnotforsamepoint (a_j, phi_j, dphi_j) (a_h, phi_h, dphi_h)

            let rec search a_cur (a_old, phi_old, dphi_old as old) =
                let a_curvec = [a_cur] |> DenseVector.ofList
                let phi_cur, dphi_cur = (phi a_curvec, dphi a_curvec)

                match phi_cur with
                | _ when (phi_cur > phi0 + c1 * a_cur * dphi0) || (phi_cur >= phi_old && a_old <> actualInitStep)
                    -> zoom (a_old, phi_old, dphi_old) (a_cur, phi_cur, dphi_cur)
                | _ -> match dphi_cur with
                       | _ when (abs dphi_cur) <= -c2 * dphi0 -> a_cur
                       | _ when dphi_cur >= 0.0 -> zoom (a_cur, phi_cur, dphi_cur) (a_old, phi_old, dphi_old)
                       | _ -> let a_new = cubicInterpolation ((a_cur, phi_cur, dphi_cur), (maxstep, phimax, dphimax))
                              match a_new with
                              | None -> rescuestepsize
                              | Some a_new ->
                                  if a_new >= maxstep then maxstep
                                  else if a_new <= a_cur then a_cur
                                  else (search a_new (a_cur, phi_cur, dphi_cur))

            match cubicInterpolation((0.0, phi0, dphi0), (maxstep, phimax, dphimax)) with
            | None -> rescuestepsize
            | Some x -> search actualInitStep (0.0, phi0, dphi0)


type NelderMead (f:(Vector<float> -> float), iteration: int, tolerance: float) =
    let defaultIteration = 100
    let defaultTolerance = 1e-3

    let mutable iteration = iteration
    let mutable tolerance = tolerance
    let mutable zdelta = 0.00025
    let mutable delta = 0.05
    let mutable rho = 1.0
    let mutable chi = 2.0
    let mutable psi = 0.5
    let mutable sigma = 0.5

    member x.Iteration with get() = iteration and set v = iteration <- if v <= 0 then defaultIteration else v
    member x.Tolerance with get() = tolerance and set v = tolerance <- if v <= 0.0 then defaultTolerance else v
    member x.ZDelta with get() = zdelta and set v = zdelta <- v
    member x.Delta with get() = delta and set v = delta <- v
    member x.Rho with get() = rho and set v = rho <- v
    member x.Chi with get() = chi and set v = chi <- v
    member x.Psi with get() = psi and set v = psi <- v
    member x.Sigma with get() = sigma and set v = sigma <- v
                
    member x.Minimize(initval: Vector<float>) =
        // Creating Simplex
        let rec loop_cs (vec: Vector<float>) acc simplex =
            if acc = vec.Count then simplex
            else vec |> Vector.mapi (fun i x -> if i = acc then match x with
                                                                | 0.0 -> zdelta
                                                                | x -> (1.0 + delta) * x
                                                else x) 
                     |> (fun s -> loop_cs vec (acc+1) (s :: simplex))
        let ss = initval :: (loop_cs initval 0 []) |> List.map (fun x -> (x, (f x)))
        
        let rec loop_it (oldsimplex: (Vector<float> * float) list) (newsimplex: (Vector<float> * float) list) (count: int) =
            let ascendingsimplex = List.sortBy (fun elem -> (snd elem)) newsimplex
            let descendingsimplex = List.rev ascendingsimplex

            let getVectorFX (s: (Vector<float> * float) list) =
                s |> List.unzip |> snd |> DenseVector.ofList

            let createNextNewSimplex (newv: (Vector<float> * float)) = 
                let otherthanlargest = List.tail descendingsimplex
                newv :: otherthanlargest

            // Checking Convergence
            let f_L2 =
                if count > 0 then
                    (getVectorFX oldsimplex) - (getVectorFX newsimplex) |> Vector.fold (fun a x -> a + x*x) 0.0 |> System.Math.Sqrt
                else
                    0.0
            if f_L2 < x.Tolerance && count > 0 then (List.head newsimplex) |> (fun x -> ((fst x), (snd x), true)) 
            else if count > x.Iteration && count > 0 then (List.head newsimplex) |> (fun x -> ((fst x), (snd x), false))          
            else
                let (best, worst, sndworst) =
                    (List.head ascendingsimplex, List.head descendingsimplex, (descendingsimplex |> List.tail |> List.head))
                let centroid = List.reduce (+) (List.tail descendingsimplex |> List.unzip |> fst)
                               |> (fun x -> x / (float (descendingsimplex.Length) - 1.0))

                // Reflect
                let r = centroid + rho * (centroid - (fst worst)) |> (fun x -> (x, f x))
            
                if (snd r) < (snd best) then       // Expand
                    let e = (fst r) + chi * ((fst r) - centroid) |> (fun x -> (x, f x))
                    loop_it ascendingsimplex (createNextNewSimplex e) (count + 1)
                else if (snd sndworst) <= (snd r) then     // Contract
                    let c =
                        if (snd r) < (snd worst) then
                            centroid + psi * ((fst r) - centroid) |> (fun x -> (x, f x))    // Outside Contract
                        else
                            centroid + psi * ((fst worst) - centroid) |> (fun x -> (x, f x))    // Inside Contract
                    if (snd c) <= (snd worst) then
                        loop_it newsimplex (createNextNewSimplex c) (count + 1)
                    else                // Shrink
                        let nnewsimplex = List.map (fun (v: (Vector<float> * float)) ->
                                                        let x = (fst best) + sigma * ((fst v) - (fst best))
                                                        (x, f x)) newsimplex
                        loop_it newsimplex nnewsimplex (count + 1)
                else
                    loop_it newsimplex (createNextNewSimplex r) (count + 1)

        loop_it [] ss 0

type QuasiNewtonMethodStatus<'a> =
    Ongoing of 'a
    | Converged of 'a
    | NotConverged of 'a
    | FunctionValueInvalid
    | GradientInvalid
    | WeightMatrixInvalid
    | LineSearchFailure

type QuasiNewtonSearchBuilder() =
    member this.Bind(result: QuasiNewtonMethodStatus<'a>, rest: 'a -> QuasiNewtonMethodStatus<'b>) =
        match result with
        | Ongoing(x)
        | Converged(x)
        | NotConverged(x) -> rest x
        | FunctionValueInvalid -> FunctionValueInvalid
        | GradientInvalid -> GradientInvalid
        | WeightMatrixInvalid -> WeightMatrixInvalid
        | LineSearchFailure -> LineSearchFailure

    member this.Return(result) = result

type BFGS (f:(Vector<float> -> float), iteration: int, tolerance: float) =  
    let defaultIteration = 100
    let defaultTolerance = 1e-3

    let mutable iteration = iteration
    let mutable tolerance = tolerance
    let mutable initialstepsize = 1.0
    let mutable maxstepsize = 5.0
    let mutable derivationmethod = (fun x -> Differentiation.Gradient(f, x))
    let mutable firsttimestepsizemuiltiplier = 1.0

    let isInvalidFloat (x: float) =
        if System.Double.IsInfinity x || System.Double.IsNaN x then true
        else false
         
    let qnsearch = new QuasiNewtonSearchBuilder()
    let derivative x =
        let tres = derivationmethod x
        if Vector.exists isInvalidFloat tres then GradientInvalid
        else Ongoing(tres)

    let BFGSWeightMatrix (w: Matrix<float>) (xds: Vector<float>) (dds: Vector<float>) =
        let tres = w + ((1.0 / (dds * xds)) * Vector.OuterProduct(dds, dds)) - ((1.0 / (xds * w * xds)) * Vector.OuterProduct((w * xds), (w * xds)))
        if Matrix.exists isInvalidFloat tres then WeightMatrixInvalid
        else Ongoing(tres)

    let linesearch r g =
        let ls = LineSearch(f, initialstepsize, maxstepsize)
        let tres = ls.Search r g
        if isInvalidFloat tres then LineSearchFailure
        else Ongoing(tres)

    member x.Iteration with get() = iteration and set v = iteration <- if v <= 0 then defaultIteration else v
    member x.Tolerance with get() = tolerance and set v = tolerance <- if v <= 0.0 then defaultTolerance else v
    member x.DerivationMethod with get() = derivationmethod and set v = derivationmethod <- v
    member x.InitialStepSize with get() = initialstepsize and set v = initialstepsize <- v
    member x.MaxStepSize with get() = maxstepsize and set v = maxstepsize <- v
    member x.FirstTimeStepSizeMultiplier with get() = firsttimestepsizemuiltiplier and set v = firsttimestepsizemuiltiplier <- v

    // To Do: Rewrite with computation Expression
    member x.Minimize(initval: Vector<float>) =
        let rec search (w: Matrix<float>) (r: Vector<float>) (g: Vector<float>) count =
            qnsearch {
                let wi = w.Inverse()

                if g.Norm(2.0) < x.Tolerance && count > 0 then return Converged(r, (f r), wi)
                else if (count >= x.Iteration) then return NotConverged(r, (f r), wi)
                else let! step = linesearch r ((-1.0) * (wi * g))
                     let newr = r - step * (wi * g)
                     let! newg = derivative newr
                     let! neww = BFGSWeightMatrix w (newr - r) (newg - g)
                     return search neww newr newg (count + 1)
            }

        qnsearch {
            let initx = initval.Clone()
            let! initd = derivative initx
            let w = let mult = if firsttimestepsizemuiltiplier = 1.0 then 1.0
                               else (1.0 / firsttimestepsizemuiltiplier) * initd.Norm(2.0)
                    mult * DenseMatrix.Identity(initx.Count)

            return search w initx initd 0
        }
