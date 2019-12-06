module Divanov.FsLearn.Test.Clustering

open Expecto
open MathNet.Numerics.LinearAlgebra
open Divanov.FsLearn.Clustering
open Divanov.FsLearn.Common

[<Tests>]
let tests =
  testList "Clustering" [
    testCase "KMeans" <| fun _ ->
      let sr = Divanov.Prelude.Misc.SyncRandom(System.Random())
      let data =
        [| 1 .. 100 |]
        |> Array.map (fun x -> [| 1 .. 5 |] |> Array.map (fun _ -> float x / 10.0 |> round |> (*) 10.0 ))
        |> Array.map (Array.map (fun x -> x + sr.NextDouble()))
        |> DenseMatrix.ofRowArrays
      let parameters = { KMeans.Params.Empty with MaxIter = 222 }
      let fittedCentroids = KMeans.fit parameters 10 data
      let classifications = Result.bind (KMeans.classify data) fittedCentroids
      Expect.isOk classifications "Should compute k-means"
      let counts = Result.map (fun x -> Array.countBy id x |> Array.map snd) classifications
      Result.map
        (fun x -> Expect.all x ((<) 0) "All clusters contain 1+ elements")
        counts
      |> ignore
      Result.map
        (fun x -> Expect.all x ((>) 33) "No clusters of 33+ elements")
        counts
      |> ignore
  ]
