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

module Spectral =
  type Params = KMeans.Params
  
  let isNonNegative (mat: Matrix<double>): bool =
    mat.ForAll(System.Func<double, bool>((<=) 0.0))

  let private unsafeDMatrix(mat: Matrix<double>): Matrix<double> =
    [| 0 .. mat.ColumnCount - 1 |]
    |> Array.map (fun i -> (mat.Column i).Sum())
    |> SparseMatrix.ofDiagArray
    

  let dMatrix (mat: Matrix<double>): Result<Matrix<double>, System.Exception> =
    if not (isNonNegative mat)
    then Error <| System.Exception("distance matrix with negative values")
    elif mat.ColumnCount = mat.RowCount
    then Error <| System.Exception "asymmetric distance matrix"
    elif not <| mat.IsSymmetric()
    then Error <| System.Exception "asymmetric distance matrix"
    else Ok <| unsafeDMatrix mat

  let lMatrix (mat: Matrix<double>): Result<Matrix<double>, System.Exception> =
    let d = dMatrix mat
    Result.map (fun dmat -> dmat - mat) d

  let lRw (mat: Matrix<double>): Result<Matrix<double>, System.Exception> =
    let identity = SparseMatrix.identity mat.ColumnCount
    let invD = dMatrix mat |> Result.map Matrix.inverse
    Result.map (fun dinv -> identity - (dinv * mat)) invD

  let private clusterSubv (parameters: Params) (k: int) (lrw: Matrix<double>) =
    if k <= 0 then Error <| System.Exception "K lower than 0"
    else
    let eigen = Matrix.eigen lrw
    let vecs = eigen.EigenVectors
    let subVecs = vecs.SubMatrix(0, vecs.RowCount, 0, k)
    let clusters = KMeans.fit parameters k subVecs
    let classifications = Result.bind (fun centroids -> KMeans.classify subVecs centroids) clusters
    classifications

  let fit (parameters: Params) (k: int) (distances: Matrix<double>) =
    let lrw = lRw distances
    let clusters = Result.bind (clusterSubv parameters k) lrw
    clusters
