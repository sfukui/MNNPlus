namespace MathNet.Numerics.Appendix.Optimization

open MathNet.Numerics
open MathNet.Numerics.Differentiation
open MathNet.Numerics.Interpolation
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double

module AppendixFunctions =
    // $L^2$ norm of an array
    let L2NormOfArray (values: float array) =
        values |> Array.map (fun v -> v * v) |> Array.sum |> sqrt

    // $L^\infty$ norm of an array
    let InfinityNormOfArray (values: float array) =
        values |> Array.map (fun v -> abs v) |> Array.max

    // $L^\infty$ norm of a vector
    let InfinityNormOfVector (values: Vector<float>) = Vector.maxAbs values

    // Pairwise difference between two arrays.
    let PairwiseDifferenceOfTwoArrays (values1: float array) (values2: float array) : float array =
        Array.map2 (fun x y -> x - y) values1 values2

    // Inner product between two arrays.
    let InnerProductOfTwoArrays (values1: float array) (values2: float array) =
                Array.map2 (fun a1 a2 -> a1 * a2) values1 values2 |> Array.sum

    // Create matrix from 2x2 block matrices.
    let CreateMatrixFrom2x2Blocks (block11: Matrix<float>) (block12: Matrix<float>) (block21: Matrix<float>) (block22: Matrix<float>) =
        let blocks = [| [| block11 ; block12 |];
                        [| block21 ; block22 |] |]
        DenseMatrix.init (blocks.[0].[0].RowCount + blocks.[1].[0].RowCount) (blocks.[0].[1].ColumnCount + blocks.[1].[1].ColumnCount)
            (fun i j -> let (matI,subI) = if i < blocks.[0].[0].RowCount then (0,i) else (1, i - blocks.[0].[0].RowCount)
                        let (matJ,subJ) = if j < blocks.[0].[0].ColumnCount then (0,j) else (1, j - blocks.[0].[0].ColumnCount)
                        let curMat = blocks.[matI].[matJ]
                        curMat.At(subI,subJ))

    // Check whether a float number is invalid.  
    let IsInvalidFloat (x: float) =
        if System.Double.IsInfinity x || System.Double.IsNaN x then true
        else false
