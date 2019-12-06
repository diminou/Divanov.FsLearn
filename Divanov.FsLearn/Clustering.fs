module Divanov.FsLearn.Clustering

open MathNet.Numerics.LinearAlgebra
open Divanov.FsLearn.Common
open Divanov.Prelude.Misc


module KMeans =
  type Params =
    {
      MaxIter: int; MinStep: double; Rand: SyncRandom
    } with
    static member Empty =
      {
        MaxIter = System.Int32.MaxValue; MinStep = 0.0
        Rand = SyncRandom(System.Random())
      }

  let voronoi (means: Matrix<double>) (vector: Vector<double>): Result<int, System.Exception> =
    try
      let replicae =
        [1 .. means.RowCount] |> List.map (fun _ -> vector)
        |> DenseMatrix.ofRows
      let diffs = means - replicae
      let scalars = diffs.RowNorms(2.0)
      let distances = scalars.ToArray()
      Divanov.Prelude.Array.argMin distances |> Ok
    with
    | exc -> Error exc

  let private accResult (elts: Result<'a, 'b> list): Result<'a list, 'b> =
    let init: Result<'a list, 'b> = Ok []
    let update accum elt =
      accum |> Result.bind (fun x -> elt |> Result.map (fun e -> e :: x))
      |> Result.map List.rev
    List.fold update init elts
  
  let adjust (rows: Matrix<double>) (means: Matrix<double>): Result<Matrix<double>, System.Exception> =
    let avg (vecs: Vector<double> []): Vector<double> =
      let len = Array.length vecs |> double
      let summed = Array.reduce (+) vecs
      summed.Multiply(1./len)

    let assignments: Result<(int * Vector<double>) [] , System.Exception> =
      [| 0 .. rows.RowCount - 1 |]
      |> Array.Parallel.map (fun x -> voronoi means (rows.Row x) |> Result.map (fun y -> (y, rows.Row x)))
      |> Array.toList
      |> accResult
      |> Result.map Array.ofList

    let newMeans: Result<Matrix<double>, System.Exception> =
      assignments
      |> Result.map (Array.groupBy fst)
      |> Result.map (Array.map (fun (x, y) -> (x, y |> Array.map snd |> avg)))
      |> Result.map (Array.sortBy fst >> Array.map snd >> List.ofArray)
      |> Result.map (DenseMatrix.ofRows)
    
    newMeans

  let init (parameters: Params) (k: int) (data: Matrix<double>): Matrix<double> =
    let indices =
      [| 0 .. data.RowCount |]
      |> Array.sortBy (fun _ -> parameters.Rand.Next())
      |> Array.truncate k
    indices |> Array.map data.Row |> List.ofArray |> DenseMatrix.ofRows

  let fit (parameters: Params) (k: int) (data: Matrix<double>): Result<Matrix<double>, System.Exception> =
    let initial = Ok (init parameters k data)
    let stagnation (mat1: Matrix<float>) (mat2: Matrix<float>): bool =
      let diff = mat1 - mat2
      diff.Exists(System.Func<float, bool>((<) parameters.MinStep))
    let newElement (i, mat) =
      if i >= parameters.MaxIter then None
      else
        match mat with
        | Error e -> None
        | Ok o ->
          let newMat = adjust data o
          let s = Result.map (stagnation o) newMat
          match s with
          | Error e2 -> Some (newMat, (i + 1, newMat))
          | Ok b ->
            if b then None
            else Some (newMat, (i + 1, newMat))
    Seq.unfold newElement (1, initial) |> Seq.fold (fun _ elt -> elt) initial

  let classify (rows: Matrix<double>) (centroids: Matrix<double>): Result<int [], System.Exception> =
    [| 0 .. rows.RowCount - 1 |]
    |> Array.Parallel.map (fun x -> voronoi centroids (rows.Row x))
    |> Array.toList
    |> accResult
    |> Result.map Array.ofList

    
