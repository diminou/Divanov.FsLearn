module Divanov.FsLearn.Common


type NonNegativeInt(i: int) =
  member __.Value = abs i

type DoubleProportion(d: double) =
  member __.Value =
    max 0.0 d |> min 1.0

type AbsOrProp =
  | Absolute of NonNegativeInt
  | Proportion of DoubleProportion
