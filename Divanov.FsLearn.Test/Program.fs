module Divanov.FsLearn.Test.Program
// Learn more about F# at http://fsharp.org

open Expecto

[<EntryPoint>]
let main argv =
    let writeResults = TestResults.writeNUnitSummary ("TestResults.xml", "Expecto.Tests")
    let config = defaultConfig.appendSummaryHandler writeResults
    runTestsInAssembly config argv
