open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Optimization
open MathNet.Numerics.Differentiation
open MathNet.Numerics.Integration
open MathNet.Numerics.LinearAlgebra.Generic
open MathNet.Numerics.LinearAlgebra.Double
open MathNet.Numerics.Statistics.Mcmc
open MathNet.Numerics.Integration.Algorithms

// Adaptive rejection metropolis sampling from normal distribution
let normalPdf (theta: Vector<float>) (x: float)
    = (1.0 / (sqrt (2.0 * theta.[1] * System.Math.PI))) * exp (-0.5 * ((x - theta.[0])**2.0) / (theta.[1]**2.0))
let normalPdfLn (theta: Vector<float>) (x: float)
    = -log (sqrt (2.0 * theta.[1] * System.Math.PI)) + (-0.5 * ((x - theta.[0])**2.0) / (theta.[1]**2.0))
let normTheta = [1.0; 2.0] |> DenseVector.ofList
let normARMSampler = AdaptiveRejectionMetropolisSampler((normalPdfLn normTheta), -7.0, 9.0, -2.0, 4.0)
let normSample = normARMSampler.Sample(5.0, 10000)
let normMean = List.sum normSample / (List.length normSample |> float)
let normSampleV = DenseVector.ofList normSample
let normMeanV = normSampleV.Sum() / ( normSample.Length |> float)
let normVar = normSample |> List.map (fun x -> (x - normMean)**2.0) |> List.sum |> (*) (1.0 / (List.length normSample - 1 |> float))
let normSd = sqrt normVar

// Adaptive rejection metropolis sampling from generalized beta2 distribution
let gbeta2PdfLn (theta: Vector<float>) x = 
    let lnumer = (log theta.[0]) + (theta.[0] * theta.[2] - 1.0) * (log x)
    let ldenom = ( theta.[0] * theta.[2] ) * (log theta.[1]) + SpecialFunctions.BetaLn(theta.[2], theta.[3]) + 
                 ( theta.[2] + theta.[3] ) * (log ( 1.0 + (x / theta.[1])**(theta.[0]) ) )
    lnumer - ldenom
let gbeta2Theta = [1.58840897986439; 1315.45595003162; 2.12302790284458; 4.62820357269144] |> DenseVector.ofList
let gbeta2ARMSampler = new AdaptiveRejectionMetropolisSampler((gbeta2PdfLn gbeta2Theta), 0.0, 40000.0, 200.0, 1500.0)
let gbeta2Sample = gbeta2ARMSampler.Sample(800.0, 200000)
let gbeta2SimMean = List.average gbeta2Sample
let gbeta2SimVar = gbeta2Sample |> List.map (fun x -> (x - gbeta2SimMean)**2.0) |> List.sum
                |> (*) (1.0 / ((List.length gbeta2Sample |> float) - 1.0))
let gbeta2SimSd = sqrt gbeta2SimVar

let deQuadrature = new DoubleExponentialTransformation()
let gbeta2Mean = deQuadrature.Integrate((fun x -> x * exp(gbeta2PdfLn gbeta2Theta x)), 0.0, 30000.0, 0.001)
let gbeta2Var = deQuadrature.Integrate((fun x -> (x - gbeta2Mean)**2.0 * exp(gbeta2PdfLn gbeta2Theta x)), 0.0, 30000.0, 0.001)
let gbeta2Sd = sqrt gbeta2Var




// Optimization: Nelder-Mead and BFGS
let targetsf (x: Generic.Vector<float>) =
    if System.Math.Abs(x.[0]) >= 1.0 then System.Double.NaN else System.Math.Sqrt(x.[0])

let x = [0.5] |> Double.DenseVector.ofList

let dif = Differentiation.Gradient(targetsf, x) 

let targetfunction (x: Generic.Vector<float>) = (x.[0] - 1.0) * (x.[0] - 1.0) + (x.[1] - 1.0) * (x.[1] - 1.0)
// Minimized at (x,y) = (1.0,1.0)
// df / dx = 4.0 and df / dy = 4.0 at (3.0,3.0)
let dataset = [(20.5, 12.0); (31.5, 16.0); (47.7, 18.0); (26.2, 16.0); (44.0, 12.0)]
let loglikelihood (p: Generic.Vector<float>) = dataset
                                               |> List.fold (fun acc (x,y) -> acc + (-0.5) * (System.Math.Log(2.0 * System.Math.PI) + (System.Math.Log(p.[2]))) - ( (y - p.[0] - p.[1] * x)**2.0 / ( 2.0 * p.[2] ) ) ) 0.0 |> (*) (-1.0)

let nm = NelderMead(targetfunction, 100, 0.001)
let bfgs = BFGS(targetfunction, 100, 0.1)
let init = [0.1;10.0] |> Double.DenseVector.ofList
let resnm = nm.Minimize(init)
let resbfgs = bfgs.Minimize(init)

let nm2 = NelderMead(loglikelihood, 10, 0.1)
let bfgs2 = BFGS(loglikelihood, 1000, 0.1)
let init2 = [0.1;0.1;10.0] |> Double.DenseVector.ofList
let (p2, y2, _) = nm2.Minimize(init2)
do bfgs2.MaxStepSize <- 5.0
let resbfgs2 = bfgs2.Minimize(p2)

let df1 = Differentiation.Gradient(targetfunction, init, ([0.01;0.01] |> Double.DenseVector.ofList))
let df2 = Differentiation.Gradient(targetfunction, init)

let funcA (x: Generic.Vector<float>) = ((x.[0] + x.[0]) + 3.0) * x.[0] + 100000.0 * System.Math.Exp(x.[1]) + 1.0
// Hessian = (4.0, 0), (0, 100000 * e^y)
let funcB (x: Generic.Vector<float>) = (1.0 / System.Math.Sqrt(2.0 * System.Math.PI) * 2.0 * x.[0])
                                        * System.Math.Exp((-1.0) * (System.Math.Log(x.[0] - 10.0))**2.0 / (2.0 * 2.0 ** 2.0))

let df3 = Differentiation.Gradient(funcA, [1000000.0;1.0] |> Double.DenseVector.ofList)
let df4 = Differentiation.Gradient(funcA, [0.1;10.0] |> Double.DenseVector.ofList)
//let h1 = Differentiation.Hessian(funcA, [0.1;10.0] |> Double.DenseVector.ofList)
//let h2 = Differentiation.Hessian(funcA, [1000000.0;10.0] |> Double.DenseVector.ofList)

let funcA2 (x: float) = [x;1.0] |> Double.DenseVector.ofList |> funcA
let int1 = MathNet.Numerics.Integration.Integrate.OnClosedInterval((fun x -> funcA2 x), 1.0, 10.0)


printfn "End of the program."

