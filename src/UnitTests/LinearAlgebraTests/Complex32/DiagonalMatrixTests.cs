// <copyright file="DiagonalMatrixTests.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
//
// Copyright (c) 2009-2013 Math.NET
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
// </copyright>

using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex32;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.Complex32
{
    using Numerics;

    /// <summary>
    /// Diagonal matrix tests.
    /// </summary>
    public class DiagonalMatrixTests : MatrixTests
    {
        /// <summary>
        /// Setup test matrices.
        /// </summary>
        [SetUp]
        public override void SetupMatrices()
        {
            TestData2D = new Dictionary<string, Complex32[,]>
                {
                    {"Singular3x3", new[,] {{new Complex32(1.0f, 1), Complex32.Zero, Complex32.Zero}, {Complex32.Zero, Complex32.Zero, Complex32.Zero}, {Complex32.Zero, Complex32.Zero, new Complex32(3.0f, 1)}}},
                    {"Square3x3", new[,] {{new Complex32(-1.1f, 1), Complex32.Zero, Complex32.Zero}, {Complex32.Zero, new Complex32(1.1f, 1), Complex32.Zero}, {Complex32.Zero, Complex32.Zero, new Complex32(6.6f, 1)}}},
                    {"Square4x4", new[,] {{new Complex32(-1.1f, 1), Complex32.Zero, Complex32.Zero, Complex32.Zero}, {Complex32.Zero, new Complex32(1.1f, 1), Complex32.Zero, Complex32.Zero}, {Complex32.Zero, Complex32.Zero, new Complex32(6.2f, 1), Complex32.Zero}, {Complex32.Zero, Complex32.Zero, Complex32.Zero, new Complex32(-7.7f, 1)}}},
                    {"Singular4x4", new[,] {{new Complex32(-1.1f, 1), Complex32.Zero, Complex32.Zero, Complex32.Zero}, {Complex32.Zero, new Complex32(-2.2f, 1), Complex32.Zero, Complex32.Zero}, {Complex32.Zero, Complex32.Zero, Complex32.Zero, Complex32.Zero}, {Complex32.Zero, Complex32.Zero, Complex32.Zero, new Complex32(-4.4f, 1)}}},
                    {"Tall3x2", new[,] {{new Complex32(-1.1f, 1), Complex32.Zero}, {Complex32.Zero, new Complex32(1.1f, 1)}, {Complex32.Zero, Complex32.Zero}}},
                    {"Wide2x3", new[,] {{new Complex32(-1.1f, 1), Complex32.Zero, Complex32.Zero}, {Complex32.Zero, new Complex32(1.1f, 1), Complex32.Zero}}}
                };

            TestMatrices = new Dictionary<string, Matrix<Complex32>>();
            foreach (var name in TestData2D.Keys)
            {
                TestMatrices.Add(name, CreateMatrix(TestData2D[name]));
            }
        }

        /// <summary>
        /// Creates a matrix for the given number of rows and columns.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="columns">The number of columns.</param>
        /// <returns>A matrix with the given dimensions.</returns>
        protected override Matrix<Complex32> CreateMatrix(int rows, int columns)
        {
            return new DiagonalMatrix(rows, columns);
        }

        /// <summary>
        /// Creates a matrix from a 2D array.
        /// </summary>
        /// <param name="data">The 2D array to create this matrix from.</param>
        /// <returns>A matrix with the given values.</returns>
        protected override Matrix<Complex32> CreateMatrix(Complex32[,] data)
        {
            return DiagonalMatrix.OfArray(data);
        }

        /// <summary>
        /// Creates a vector of the given size.
        /// </summary>
        /// <param name="size">The size of the vector to create.
        /// </param>
        /// <returns>The new vector. </returns>
        protected override Vector<Complex32> CreateVector(int size)
        {
            return new SparseVector(size);
        }

        /// <summary>
        /// Creates a vector from an array.
        /// </summary>
        /// <param name="data">The array to create this vector from.</param>
        /// <returns>The new vector. </returns>
        protected override Vector<Complex32> CreateVector(Complex32[] data)
        {
            return SparseVector.OfEnumerable(data);
        }

        /// <summary>
        /// Can create a matrix from a diagonal array.
        /// </summary>
        [Test]
        public void CanCreateMatrixFromDiagonalArray()
        {
            var testData = new Dictionary<string, Matrix<Complex32>>
                {
                    {"Singular3x3", new DiagonalMatrix(3, 3, new[] {new Complex32(1.0f, 1), Complex32.Zero, new Complex32(3.0f, 1)})},
                    {"Square3x3", new DiagonalMatrix(3, 3, new[] {new Complex32(-1.1f, 1), new Complex32(1.1f, 1), new Complex32(6.6f, 1)})},
                    {"Square4x4", new DiagonalMatrix(4, 4, new[] {new Complex32(-1.1f, 1), new Complex32(1.1f, 1), new Complex32(6.2f, 1), new Complex32(-7.7f, 1)})},
                    {"Tall3x2", new DiagonalMatrix(3, 2, new[] {new Complex32(-1.1f, 1), new Complex32(1.1f, 1)})},
                    {"Wide2x3", new DiagonalMatrix(2, 3, new[] {new Complex32(-1.1f, 1), new Complex32(1.1f, 1)})},
                };

            foreach (var name in testData.Keys)
            {
                Assert.That(testData[name], Is.EqualTo(TestMatrices[name]));
            }
        }

        /// <summary>
        /// Matrix from array is a reference.
        /// </summary>
        [Test]
        public void MatrixFrom1DArrayIsReference()
        {
            var data = new[] {new Complex32(1.0f, 1), new Complex32(2.0f, 1), new Complex32(3.0f, 1), new Complex32(4.0f, 1), new Complex32(5.0f, 1)};
            var matrix = new DiagonalMatrix(5, 5, data);
            matrix[0, 0] = new Complex32(10.0f, 1);
            Assert.AreEqual(new Complex32(10.0f, 1), data[0]);
        }

        /// <summary>
        /// Can create a matrix from two-dimensional array.
        /// </summary>
        /// <param name="name">Matrix name.</param>
        [TestCase("Singular3x3")]
        [TestCase("Singular4x4")]
        [TestCase("Square3x3")]
        [TestCase("Square4x4")]
        [TestCase("Tall3x2")]
        [TestCase("Wide2x3")]
        public void CanCreateMatrixFrom2DArray(string name)
        {
            var matrix = DiagonalMatrix.OfArray(TestData2D[name]);
            for (var i = 0; i < TestData2D[name].GetLength(0); i++)
            {
                for (var j = 0; j < TestData2D[name].GetLength(1); j++)
                {
                    Assert.AreEqual(TestData2D[name][i, j], matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Can create a matrix with uniform values.
        /// </summary>
        [Test]
        public void CanCreateMatrixWithUniformValues()
        {
            var matrix = new DiagonalMatrix(10, 10, new Complex32(10.0f, 1));
            for (var i = 0; i < matrix.RowCount; i++)
            {
                Assert.AreEqual(matrix[i, i], new Complex32(10.0f, 1));
            }
        }

        /// <summary>
        /// Can create an identity matrix.
        /// </summary>
        [Test]
        public void CanCreateIdentity()
        {
            var matrix = DiagonalMatrix.CreateIdentity(5);
            for (var i = 0; i < matrix.RowCount; i++)
            {
                for (var j = 0; j < matrix.ColumnCount; j++)
                {
                    Assert.AreEqual(i == j ? Complex32.One : Complex32.Zero, matrix[i, j]);
                }
            }
        }

        /// <summary>
        /// Identity with wrong order throws <c>ArgumentOutOfRangeException</c>.
        /// </summary>
        /// <param name="order">The size of the square matrix</param>
        [TestCase(0)]
        [TestCase(-1)]
        public void IdentityWithWrongOrderThrowsArgumentOutOfRangeException(int order)
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => DiagonalMatrix.CreateIdentity(order));
        }

        /// <summary>
        /// Can multiply a matrix with matrix.
        /// </summary>
        /// <param name="nameA">Matrix A name.</param>
        /// <param name="nameB">Matrix B name.</param>
        public override void CanMultiplyMatrixWithMatrixIntoResult(string nameA, string nameB)
        {
            var matrixA = TestMatrices[nameA];
            var matrixB = TestMatrices[nameB];
            var matrixC = new SparseMatrix(matrixA.RowCount, matrixB.ColumnCount);
            matrixA.Multiply(matrixB, matrixC);

            Assert.AreEqual(matrixC.RowCount, matrixA.RowCount);
            Assert.AreEqual(matrixC.ColumnCount, matrixB.ColumnCount);

            for (var i = 0; i < matrixC.RowCount; i++)
            {
                for (var j = 0; j < matrixC.ColumnCount; j++)
                {
                    AssertHelpers.AlmostEqualRelative(matrixA.Row(i)*matrixB.Column(j), matrixC[i, j], 15);
                }
            }
        }

        /// <summary>
        /// Permute matrix rows throws <c>InvalidOperationException</c>.
        /// </summary>
        [Test]
        public void PermuteMatrixRowsThrowsInvalidOperationException()
        {
            var matrixp = CreateMatrix(TestData2D["Singular3x3"]);
            var permutation = new Permutation(new[] {2, 0, 1});
            Assert.Throws<InvalidOperationException>(() => matrixp.PermuteRows(permutation));
        }

        /// <summary>
        /// Permute matrix columns throws <c>InvalidOperationException</c>.
        /// </summary>
        [Test]
        public void PermuteMatrixColumnsThrowsInvalidOperationException()
        {
            var matrixp = CreateMatrix(TestData2D["Singular3x3"]);
            var permutation = new Permutation(new[] {2, 0, 1});
            Assert.Throws<InvalidOperationException>(() => matrixp.PermuteColumns(permutation));
        }

        /// <summary>
        /// Can pointwise divide matrices into a result matrix.
        /// </summary>
        public override void CanPointwiseDivideIntoResult()
        {
            foreach (var data in TestMatrices.Values)
            {
                var other = data.Clone();
                var result = data.Clone();
                data.PointwiseDivide(other, result);
                var min = Math.Min(data.RowCount, data.ColumnCount);
                for (var i = 0; i < min; i++)
                {
                    Assert.AreEqual(data[i, i]/other[i, i], result[i, i]);
                }

                result = data.PointwiseDivide(other);
                for (var i = 0; i < min; i++)
                {
                    Assert.AreEqual(data[i, i]/other[i, i], result[i, i]);
                }
            }
        }

        /// <summary>
        /// Can compute Frobenius norm.
        /// </summary>
        public override void CanComputeFrobeniusNorm()
        {
            var matrix = TestMatrices["Square3x3"];
            var denseMatrix = DenseMatrix.OfArray(TestData2D["Square3x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.FrobeniusNorm(), matrix.FrobeniusNorm(), 14);

            matrix = TestMatrices["Wide2x3"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Wide2x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.FrobeniusNorm(), matrix.FrobeniusNorm(), 14);

            matrix = TestMatrices["Tall3x2"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Tall3x2"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.FrobeniusNorm(), matrix.FrobeniusNorm(), 14);
        }

        /// <summary>
        /// Can compute Infinity norm.
        /// </summary>
        public override void CanComputeInfinityNorm()
        {
            var matrix = TestMatrices["Square3x3"];
            var denseMatrix = DenseMatrix.OfArray(TestData2D["Square3x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.InfinityNorm(), matrix.InfinityNorm(), 14);

            matrix = TestMatrices["Wide2x3"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Wide2x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.InfinityNorm(), matrix.InfinityNorm(), 14);

            matrix = TestMatrices["Tall3x2"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Tall3x2"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.InfinityNorm(), matrix.InfinityNorm(), 14);
        }

        /// <summary>
        /// Can compute L1 norm.
        /// </summary>
        public override void CanComputeL1Norm()
        {
            var matrix = TestMatrices["Square3x3"];
            var denseMatrix = DenseMatrix.OfArray(TestData2D["Square3x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.L1Norm(), matrix.L1Norm(), 14);

            matrix = TestMatrices["Wide2x3"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Wide2x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.L1Norm(), matrix.L1Norm(), 14);

            matrix = TestMatrices["Tall3x2"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Tall3x2"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.L1Norm(), matrix.L1Norm(), 14);
        }

        /// <summary>
        /// Can compute L2 norm.
        /// </summary>
        public override void CanComputeL2Norm()
        {
            var matrix = TestMatrices["Square3x3"];
            var denseMatrix = DenseMatrix.OfArray(TestData2D["Square3x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.L2Norm(), matrix.L2Norm(), 14);

            matrix = TestMatrices["Wide2x3"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Wide2x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.L2Norm(), matrix.L2Norm(), 14);

            matrix = TestMatrices["Tall3x2"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Tall3x2"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.L2Norm(), matrix.L2Norm(), 14);
        }

        /// <summary>
        /// Can compute determinant.
        /// </summary>
        [Test]
        public void CanComputeDeterminant()
        {
            var matrix = TestMatrices["Square3x3"];
            var denseMatrix = DenseMatrix.OfArray(TestData2D["Square3x3"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.Determinant(), matrix.Determinant(), 14);

            matrix = TestMatrices["Square4x4"];
            denseMatrix = DenseMatrix.OfArray(TestData2D["Square4x4"]);
            AssertHelpers.AlmostEqualRelative(denseMatrix.Determinant(), matrix.Determinant(), 14);
        }

        /// <summary>
        /// Determinant of  non-square matrix throws <c>ArgumentException</c>.
        /// </summary>
        [Test]
        public void DeterminantNotSquareMatrixThrowsArgumentException()
        {
            var matrix = TestMatrices["Tall3x2"];
            Assert.Throws<ArgumentException>(() => matrix.Determinant());
        }

        /// <summary>
        /// Can check if a matrix is symmetric.
        /// </summary>
        [Test]
        public override void CanCheckIfMatrixIsSymmetric()
        {
            var matrix = TestMatrices["Square3x3"];
            Assert.IsTrue(matrix.IsSymmetric);
        }
    }
}
