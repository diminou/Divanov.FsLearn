module Divanov.FsLearn.MatrixFactoring

open MathNet.Numerics.LinearAlgebra
open Divanov.Prelude
open Common



module SVD =
  type Params =
    {
      MaxDimension: Common.AbsOrProp
      SR: Misc.SyncRandom
    } with
    static member Empty =
      {
        MaxDimension = Proportion <| DoubleProportion 1.0
        SR = Misc.SyncRandom(System.Random())
      }

  type Output<'a when
                'a: struct and
                'a: (new: unit -> 'a) and
                'a :> System.IEquatable<'a> and
                'a :> System.IFormattable and
                'a :> System.ValueType > =
    {
      U: Matrix<'a>; S: Matrix<'a>; Vt: Matrix<'a>
    } with
      member o.Approximate: Matrix<'a> = o.U * o.S * o.Vt

  module Output =
    let truncate (maxDim: int) (o: Output<'a>) =
      let S' = o.S.SubMatrix(0, min maxDim o.S.RowCount, 0, min maxDim o.S.ColumnCount)
      let U' = o.U.SubMatrix(0, o.U.RowCount, 0, min maxDim o.U.ColumnCount)
      let Vt' = o.Vt.SubMatrix(0, min maxDim o.Vt.RowCount, 0, o.Vt.ColumnCount)
      { U = U'; S = S'; Vt = Vt' }

  let factor (parameters: Params) (mat: Matrix<'a>): Output<'a> =
    let svd = mat.Svd(true)
    match parameters.MaxDimension with
    | Absolute a -> Output.truncate a.Value { U = svd.U; S = svd.W; Vt = svd.VT }
    | Proportion p ->
      let initialDim = min svd.W.ColumnCount svd.W.RowCount
      let targetDim = float initialDim |> (*) p.Value |> round |> int
      Output.truncate targetDim { U = svd.U; S = svd.W; Vt = svd.VT }

  let recommendUs (maxItems: NonNegativeInt ) (output: Output<'a>) (v: Vector<'a>): int [] =
    let lower = 0
    let higher = output.U.RowCount
    let maxItems' = max lower (min maxItems.Value higher)
    let similarities = output.U * v
    similarities.ToArray()
    |> Array.argMaxN maxItems'

  let recommendVs (maxItems: NonNegativeInt) (output: Output<'a>) (u: Vector<'a>): int [] =
    let higher = output.Vt.RowCount
    let maxItems' = max 0 (min maxItems.Value higher)
    let similarities = output.Vt * u
    similarities.ToArray()
    |> Array.argMaxN maxItems'

  let fetchUs (maxItems: NonNegativeInt) (output: Output<'a>) (v: NonNegativeInt): Result<int [], string> =
    let higher = output.Vt.RowCount
    if v.Value > higher
    then Error "index v is higher than there are elements in Vt"
    else Ok (recommendUs maxItems output (output.Vt.Row v.Value))

  let fetchVs (maxItems: NonNegativeInt) (output: Output<'a>) (u: NonNegativeInt): Result<int [], string> =
    let higher = output.U.RowCount
    if u.Value > higher
    then Error "index u is higher than there are elements in U"
    else Ok (recommendUs maxItems output (output.U.Row u.Value))

module ALS =
  type Params =
    {
      MaxDimension: NonNegativeInt
    }
  type Output<'a when
                'a: struct and
                'a: (new: unit -> 'a) and
                'a :> System.IEquatable<'a> and
                'a :> System.IFormattable and
                'a :> System.ValueType > =
    {
        U: Matrix<'a>; V: Matrix<'a>
    } with
      static member Init(parameters: Params, target: Matrix<double>): Output<double> =
          let matU: Matrix<double> =
            DenseMatrix.randomStandard parameters.MaxDimension.Value target.RowCount
            |> Matrix.map abs
          let matV: Matrix<double> =
            DenseMatrix.randomStandard parameters.MaxDimension.Value target.ColumnCount
            |> Matrix.map abs
          { U = matU; V = matV }

      member o.Approximate = o.U.Transpose() * o.V

  let private reg1 (lambda: 'a) (data: Vector<'a>) (target: 'a) =
    let xtx = (data.OuterProduct data) + (SparseMatrix.diag data.Count lambda)
    (xtx.Inverse() * data).Multiply target

  let private regress (lambda: 'a) (data: Matrix<'a>) (target: Matrix<'a>) =
    let xtx = data.Transpose() * data + (SparseMatrix.diag data.ColumnCount lambda)
    xtx.Inverse() * data.Transpose() * target

  let private uStep (lambda: 'a) (coeffMat: Matrix<'a>) (dataMat: Matrix<'a>) (target: Matrix<'a>) =
    let localCoeff: Matrix<'a> = coeffMat
    let dataMatCols = [ 0 .. dataMat.ColumnCount - 1 ]
    let coeffMatCols = [| 0 .. coeffMat.ColumnCount - 1 |]
    let unitaryXtx index =
      let col = dataMat.Column index
      col.OuterProduct col

    let xtxs = dataMatCols |> List.map (fun i -> Lazy.Create(fun () -> unitaryXtx i))
    let fstXtx = (List.head xtxs).Value
    let inv: Matrix<'a> =
      List.tail xtxs
      |> List.fold (fun accum elt -> accum + elt.Value) fstXtx
      |> (+) (SparseMatrix.diag fstXtx.ColumnCount lambda)
      |> fun x -> x.Inverse()
    let xr coeffIndex dataIndex =
      let col: Vector<'a> = dataMat.Column dataIndex
      let r: 'a = target.At(coeffIndex, dataIndex)
      col.Multiply r
    let xrs coeffIndex =
      let fstXr = xr coeffIndex dataMatCols.Head
      let rest = List.tail dataMatCols
      rest |> List.fold (fun accum elt -> accum + (xr coeffIndex elt)) fstXr
    coeffMatCols |> Array.Parallel.iter (fun i -> localCoeff.[*,i] <- (inv * (xrs i)))
    localCoeff
    

  let step (lambda: 'a) (o: Output<'a>) (target: Matrix<'a>) =
    let mutable u = o.U
    let mutable v = o.V
    let targetT = target.Transpose()
    u <- uStep lambda u v target
    v <- uStep lambda v u targetT
    { U = u ; V = v }

  let iterateN (n: int) (lambda: 'a) (o: Output<'a>) (target: Matrix<'a>) =
    [1 .. n]
    |> List.fold (fun accum _ -> step lambda accum target) o
   
module RandomProjections =
  
  type Params =
    {
      MaxDimension: NonNegativeInt
      Rg: System.Random
    } with

    static member Default(data: Matrix<float>, expectedDim: NonNegativeInt): Params =
      let REL_EPSILON = 0.05
      let absEpsilon = data.FrobeniusNorm() * REL_EPSILON
      let recommendedT = (float expectedDim.Value) / (absEpsilon * absEpsilon) |> round |> int
      { MaxDimension = NonNegativeInt recommendedT ; Rg = System.Random() }
  
  let signMatrix (parameters: Params) (dataDim: NonNegativeInt): Matrix<float> =
      
    let coeff (b: bool): float =
      let t = float parameters.MaxDimension.Value
      if b then 1.0/(sqrt t) else -1.0/(sqrt t)
    
    let signBits =
      DenseMatrix.init
        dataDim.Value
        parameters.MaxDimension.Value
        (fun _ _ -> parameters.Rg.Next(2) |> (<) 0 |> coeff)
    
    signBits

  let project (parameters: Params) (data: Matrix<float>): Matrix<float> =
    let projectionMatrix = signMatrix parameters (NonNegativeInt data.ColumnCount)
    data * projectionMatrix
