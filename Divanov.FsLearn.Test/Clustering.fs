module Divanov.FsLearn.Test.Clustering

open Expecto
open MathNet.Numerics.LinearAlgebra
open Divanov.FsLearn.Clustering
open Divanov.FsLearn.Common


let sr = Divanov.Prelude.Misc.SyncRandom(System.Random())
let data =
        [| 1 .. 100 |]
        |> Array.map (fun x -> [| 1 .. 5 |] |> Array.map (fun _ -> float x / 10.0 |> round |> (*) 10.0 ))
        |> Array.map (Array.map (fun x -> x + sr.NextDouble()))
        |> DenseMatrix.ofRowArrays


[<Tests>]
let tests =
  testList "Clustering" [
    testCase "KMeans" <| fun _ ->
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
        (fun x -> Expect.all x ((>=) 50) "No clusters of 50+ elements")
        counts
      |> ignore

    testCase "Spectral clustering" <| fun _ ->
      let arrData = data.ToRowArrays()
      let dist (v1: double []) (v2: double []) =
        Array.zip v1 v2 |> Array.map (fun (x, y) -> (x - y) * (x - y)) |> Array.sum
      let dists =
        arrData
        |> Array.map (fun p -> arrData |> Array.map (fun q -> dist p q))
        |> DenseMatrix.ofRowArrays
      let parameters = { Spectral.Params.Empty with MaxIter = 222 }
      let classifications = Divanov.FsLearn.Clustering.Spectral.fitPredict parameters 10 dists
      let counts = Result.map (fun x -> Array.countBy id x |> Array.map snd) classifications
      Expect.isOk classifications "Should compute spectral clustering"
      Result.map
        (fun x -> Expect.all x ((<) 0) "All clusters contain 1+ elements")
        counts |> ignore
      Result.map
        (fun x -> Expect.all x ((>=) 40) "No clusters of 40+ elements")
        counts |> ignore
  ]
