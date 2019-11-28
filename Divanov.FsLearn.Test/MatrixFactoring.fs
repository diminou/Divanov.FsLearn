module Divanov.FsLearn.Test.MatrixFactoring

open Expecto
open MathNet.Numerics.LinearAlgebra
open Divanov.FsLearn.MatrixFactoring
open Divanov.FsLearn.Common

[<Tests>]
let tests = 
  testList "Matrix factoring" [
    testCase "SVD float" <| fun _ ->
      let mat1: Matrix<double> = DenseMatrix.randomStandard 10 15 |> Matrix.map abs
      let ps = { SVD.Params.Empty with MaxDimension = Absolute (NonNegativeInt 3) }
      let svd = SVD.factor ps mat1
      let approx = svd.Approximate
      let l2 = approx - mat1 |> (fun x -> x.L2Norm())
      Expect.isLessThan l2 5.0 "The approximate is of low quality for SVD"

    testCase "ALS float" <| fun _ ->
      let mat1: Matrix<double> = DenseMatrix.randomStandard 10 15 |> Matrix.map abs
      let ps: ALS.Params = { MaxDimension = NonNegativeInt 3 }
      let initialGuess = ALS.Output<double>.Init(ps, mat1)
      let finalGuess = ALS.iterateN 22 0.002 initialGuess mat1
      let approx = finalGuess.Approximate
      let l2 = approx - mat1 |> (fun x -> x.L2Norm())
      Expect.isLessThan l2 5.0 "The approximate is of low quality for ALS"
  ]
