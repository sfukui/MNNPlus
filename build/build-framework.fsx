//  __  __       _   _       _   _ ______ _______
// |  \/  |     | | | |     | \ | |  ____|__   __|
// | \  / | __ _| |_| |__   |  \| | |__     | |
// | |\/| |/ _` | __| '_ \  | . ` |  __|    | |
// | |  | | (_| | |_| | | |_| |\  | |____   | |
// |_|  |_|\__,_|\__|_| |_(_)_| \_|______|  |_|
//
// Math.NET Numerics - https://numerics.mathdotnet.com
// Copyright (c) Math.NET - Open Source MIT/X11 License
//
// Build Framework using FAKE (http://fsharp.github.io/FAKE)
//

module BuildFramework

#I "../packages/build/FAKE/tools"
#r "../packages/build/FAKE/tools/FakeLib.dll"

open FSharp.Core
open Fake
open Fake.ReleaseNotesHelper
open Fake.StringHelper
open Fake.Testing.NUnit3
open System
open System.IO

let rootDir = Path.GetFullPath (Path.Combine (__SOURCE_DIRECTORY__ + "/../"))
Environment.CurrentDirectory <- rootDir
trace rootDir

let header = ReadFile(__SOURCE_DIRECTORY__ </> __SOURCE_FILE__) |> Seq.take 10 |> Seq.map (fun s -> s.Substring(2)) |> toLines

type Release =
    { Title: string
      AssemblyVersion: string
      PackageVersion: string
      ReleaseNotes: string
      ReleaseNotesFile: string }

type Package =
    { Id: string
      Release: Release
      Title: string
      Summary: string
      Description: string
      Tags: string
      FsLoader: bool
      Authors: string list
      Dependencies: NugetFrameworkDependencies list
      Files: (string * string option * string option) list }

type Bundle =
    { Id: string
      Release: Release
      Title: string
      Packages: Package list }

let release title releaseNotesFile : Release =
    let info = LoadReleaseNotes releaseNotesFile
    let buildPart = "0"
    let assemblyVersion = info.AssemblyVersion + "." + buildPart
    let packageVersion = info.NugetVersion
    let notes = info.Notes |> List.map (fun l -> l.Replace("*","").Replace("`","")) |> toLines
    { Title = title
      AssemblyVersion = assemblyVersion
      PackageVersion = packageVersion
      ReleaseNotes = notes
      ReleaseNotesFile = releaseNotesFile }

let traceHeader (releases:Release list) =
    trace header
    let titleLength = releases |> List.map (fun r -> r.Title.Length) |> List.max
    for release in releases do
        trace ([ " "; release.Title.PadRight titleLength; "  v"; release.PackageVersion ] |> String.concat "")
    trace ""


// --------------------------------------------------------------------------------------
// TARGET FRAMEWORKS
// --------------------------------------------------------------------------------------

let libnet35 = "lib/net35"
let libnet40 = "lib/net40"
let libnet45 = "lib/net45"
let netstandard13 = "lib/netstandard1.3"
let netstandard16 = "lib/netstandard1.6"
let netstandard20 = "lib/netstandard2.0"
let libpcl7 = "lib/portable-net45+netcore45+MonoAndroid1+MonoTouch1"
let libpcl47 = "lib/portable-net45+sl5+netcore45+MonoAndroid1+MonoTouch1"
let libpcl78 = "lib/portable-net45+netcore45+wp8+MonoAndroid1+MonoTouch1"
let libpcl259 = "lib/portable-net45+netcore45+wpa81+wp8+MonoAndroid1+MonoTouch1"
let libpcl328 = "lib/portable-net4+sl5+netcore45+wpa81+wp8+MonoAndroid1+MonoTouch1"


// --------------------------------------------------------------------------------------
// .Net SDK
// --------------------------------------------------------------------------------------

let msbuild targets configuration project =
    MSBuildHelper.build (fun p ->
        { p with
            NoLogo = true
            NodeReuse = true
            Targets = targets
            Properties = [ "Configuration", configuration ]
            RestorePackagesFlag = false
            Verbosity = Some MSBuildVerbosity.Minimal
        }) project

let dotnet workingDir command =
    DotNetCli.RunCommand
        (fun c -> { c with WorkingDir = workingDir})
        command


// --------------------------------------------------------------------------------------
// PREPARE
// --------------------------------------------------------------------------------------

let patchVersionInAssemblyInfo path (release:Release) =
    BulkReplaceAssemblyInfoVersions path (fun f ->
        { f with
            AssemblyVersion = release.AssemblyVersion
            AssemblyFileVersion = release.AssemblyVersion
            AssemblyInformationalVersion = release.PackageVersion })

let private regexes_sl = new System.Collections.Generic.Dictionary<string, System.Text.RegularExpressions.Regex>()
let private getRegexSingleLine pattern =
    match regexes_sl.TryGetValue pattern with
    | true, regex -> regex
    | _ -> (new System.Text.RegularExpressions.Regex(pattern, System.Text.RegularExpressions.RegexOptions.Singleline))
let regex_replace_singleline pattern (replacement : string) text = (getRegexSingleLine pattern).Replace(text, replacement)

let patchVersionInResource path (release:Release) =
    ReplaceInFile
        (regex_replace @"\d+\.\d+\.\d+\.\d+" release.AssemblyVersion
         >> regex_replace @"\d+,\d+,\d+,\d+" (replace "." "," release.AssemblyVersion))
        path

let patchVersionInProjectFile path (release:Release) =
    let semverSplit = release.PackageVersion.IndexOf('-')
    let prefix = if semverSplit <= 0 then release.PackageVersion else release.PackageVersion.Substring(0, semverSplit)
    let suffix = if semverSplit <= 0 then "" else release.PackageVersion.Substring(semverSplit+1)
    ReplaceInFile
        (regex_replace """\<PackageVersion\>.*\</PackageVersion\>""" (sprintf """<PackageVersion>%s</PackageVersion>""" release.PackageVersion)
        >> regex_replace """\<Version\>.*\</Version\>""" (sprintf """<Version>%s</Version>""" release.PackageVersion)
        >> regex_replace """\<AssemblyVersion\>.*\</AssemblyVersion\>""" (sprintf """<AssemblyVersion>%s</AssemblyVersion>""" release.AssemblyVersion)
        >> regex_replace """\<FileVersion\>.*\</FileVersion\>""" (sprintf """<FileVersion>%s</FileVersion>""" release.AssemblyVersion)
        >> regex_replace """\<VersionPrefix\>.*\</VersionPrefix\>""" (sprintf """<VersionPrefix>%s</VersionPrefix>""" prefix)
        >> regex_replace """\<VersionSuffix\>.*\</VersionSuffix\>""" (sprintf """<VersionSuffix>%s</VersionSuffix>""" suffix)
        >> regex_replace_singleline """\<PackageReleaseNotes\>.*\</PackageReleaseNotes\>""" (sprintf """<PackageReleaseNotes>%s</PackageReleaseNotes>""" release.ReleaseNotes))
        path


// --------------------------------------------------------------------------------------
// BUILD
// --------------------------------------------------------------------------------------

let clean project = msbuild [ "Clean" ] "Release" project
let restore project = msbuild [ "Restore" ] "Release" project
let build project = msbuild [ (if hasBuildParam "incremental" then "Build" else "Rebuild") ] "Release" project
let pack project = dotnet rootDir (sprintf "pack %s --configuration Release --no-restore --no-build" project)

//let buildConfig config subject = MSBuild "" (if hasBuildParam "incremental" then "Build" else "Rebuild") [ "Configuration", config ] subject |> ignore
//let build subject = buildConfig "Release" subject
//let buildSigned subject = buildConfig "Release-Signed" subject
let buildConfig32 config subject = MSBuild "" (if hasBuildParam "incremental" then "Build" else "Rebuild") [("Configuration", config); ("Platform","Win32")] subject |> ignore
let buildConfig64 config subject = MSBuild "" (if hasBuildParam "incremental" then "Build" else "Rebuild") [("Configuration", config); ("Platform","x64")] subject |> ignore


// --------------------------------------------------------------------------------------
// TEST
// --------------------------------------------------------------------------------------

let test testsDir testsProj framework =
    dotnet testsDir (sprintf "run -p %s --configuration Release --framework %s --no-restore --no-build" testsProj framework)


// --------------------------------------------------------------------------------------
// PACKAGES
// --------------------------------------------------------------------------------------

let provideLicense path =
    ReadFileAsString "LICENSE.md"
    |> ConvertTextToWindowsLineBreaks
    |> ReplaceFile (path </> "license.txt")

let provideReadme title (release:Release) path =
    String.concat Environment.NewLine [header; " " + title; ""; ReadFileAsString release.ReleaseNotesFile]
    |> ConvertTextToWindowsLineBreaks
    |> ReplaceFile (path </> "readme.txt")

let provideFsLoader includes path =
    // inspired by FsLab/tpetricek
    let fullScript = ReadFile "src/FSharp/MathNet.Numerics.fsx" |> Array.ofSeq
    let startIndex = fullScript |> Seq.findIndex (fun s -> s.Contains "***MathNet.Numerics.fsx***")
    let extraScript = fullScript .[startIndex + 1 ..] |> List.ofSeq
    let assemblies = [ "MathNet.Numerics.dll"; "MathNet.Numerics.FSharp.dll" ]
    let nowarn = ["#nowarn \"211\""]
    let references = [ for assembly in assemblies -> sprintf "#r \"%s\"" assembly ]
    ReplaceFile (path </> "MathNet.Numerics.fsx") (nowarn @ includes @ references @ extraScript |> toLines)

let provideFsIfSharpLoader path =
    let fullScript = ReadFile "src/FSharp/MathNet.Numerics.IfSharp.fsx" |> Array.ofSeq
    let startIndex = fullScript |> Seq.findIndex (fun s -> s.Contains "***MathNet.Numerics.IfSharp.fsx***")
    ReplaceFile (path </> "MathNet.Numerics.IfSharp.fsx") (fullScript .[startIndex + 1 ..] |> toLines)

let provideZipExtraFiles path (bundle:Bundle) =
    provideLicense path
    provideReadme (sprintf "%s v%s" bundle.Title bundle.Release.PackageVersion) bundle.Release path
    if bundle.Packages |> List.exists (fun p -> p.FsLoader) then
        let includes = [ for root in [ ""; "../"; "../../" ] -> sprintf "#I \"%sNet40\"" root ]
        provideFsLoader includes path
        provideFsIfSharpLoader path

let provideNuGetExtraFiles path (bundle:Bundle) (pack:Package) =
    provideLicense path
    provideReadme (sprintf "%s v%s" pack.Title pack.Release.PackageVersion) bundle.Release path
    if pack.FsLoader then
        let includes = [ for root in [ ""; "../"; "../../"; "../../../" ] do
                         for package in bundle.Packages do
                         yield sprintf "#I \"%spackages/%s/lib/net40/\"" root package.Id
                         yield sprintf "#I \"%spackages/%s.%s/lib/net40/\"" root package.Id package.Release.PackageVersion ]
        provideFsLoader includes path
        provideFsIfSharpLoader path

// SIGN

let sign fingerprint timeserver files =
    let fileArgs = files |> Seq.map (sprintf "\"%s\"") |> String.concat " "
    let optionsArgs = sprintf """/v /fd sha256 /sha1 "%s" /tr "%s" /td sha256""" fingerprint timeserver
    let arguments = sprintf """sign %s %s""" optionsArgs fileArgs
    let result =
        ExecProcess (fun info ->
            info.FileName <- findToolInSubPath "signtool.exe" """C:\Program Files (x86)\Windows Kits\10\bin\x64"""
            info.Arguments <- arguments) TimeSpan.MaxValue
    if result <> 0 then
        failwithf "Error during SignTool call "


// ZIP

let zip zipDir filesDir filesFilter (bundle:Bundle) =
    CleanDir "obj/Zip"
    let workPath = "obj/Zip/" + bundle.Id
    CopyDir workPath filesDir filesFilter
    provideZipExtraFiles workPath bundle
    Zip "obj/Zip/" (zipDir </> sprintf "%s-%s.zip" bundle.Id bundle.Release.PackageVersion) !! (workPath + "/**/*.*")
    CleanDir "obj/Zip"

// NUGET

let updateNuspec (pack:Package) outPath symbols updateFiles spec =
    { spec with ToolPath = "packages/build/NuGet.CommandLine/tools/NuGet.exe"
                OutputPath = outPath
                WorkingDir = "obj/NuGet"
                Version = pack.Release.PackageVersion
                ReleaseNotes = pack.Release.ReleaseNotes
                Project = pack.Id
                Title = pack.Title
                Summary = pack.Summary
                Description = pack.Description
                Tags = pack.Tags
                Authors = pack.Authors
                DependenciesByFramework = pack.Dependencies
                SymbolPackage = symbols
                Files = updateFiles pack.Files
                Publish = false }

let nugetPack (bundle:Bundle) outPath =
    CleanDir "obj/NuGet"
    for pack in bundle.Packages do
        provideNuGetExtraFiles "obj/NuGet" bundle pack
        let withLicenseReadme f = [ "license.txt", None, None; "readme.txt", None, None; ] @ f
        let withoutSymbolsSources f =
            List.choose (function | (_, Some (target:string), _) when target.StartsWith("src") -> None
                                  | (s, t, None) -> Some (s, t, Some ("**/*.pdb"))
                                  | (s, t, Some e) -> Some (s, t, Some (e + ";**/*.pdb"))) f
        // first pass - generates symbol + normal package. NuGet does drop the symbols from the normal package, but unfortunately not the sources.
        // NuGet (updateNuspec pack outPath NugetSymbolPackage.Nuspec withLicenseReadme) "build/MathNet.Numerics.nuspec"
        // second pass - generate only normal package, again, but this time explicitly drop the sources (and the debug symbols)
        NuGet (updateNuspec pack outPath NugetSymbolPackage.None (withLicenseReadme >> withoutSymbolsSources)) "build/MathNet.Numerics.nuspec"
        CleanDir "obj/NuGet"

let nugetPackExtension (bundle:Bundle) outPath =
    CleanDir "obj/NuGet"
    for pack in bundle.Packages do
        provideNuGetExtraFiles "obj/NuGet" bundle pack
        let withLicenseReadme f = [ "license.txt", None, None; "readme.txt", None, None; ] @ f
        NuGet (updateNuspec pack outPath NugetSymbolPackage.None withLicenseReadme) "build/MathNet.Numerics.Extension.nuspec"
        CleanDir "obj/NuGet"


// --------------------------------------------------------------------------------------
// Documentation
// --------------------------------------------------------------------------------------

let provideDocExtraFiles extraDocs (releases:Release list) =
    for (fileName, docName) in extraDocs do CopyFile ("docs/content" </> docName) fileName
    let menu = releases |> List.map (fun r -> sprintf "[%s](%s)" r.Title (r.ReleaseNotesFile |> replace "RELEASENOTES" "ReleaseNotes" |> replace ".md" ".html")) |> String.concat " | "
    for release in releases do
        String.concat Environment.NewLine
          [ "# " + release.Title + " Release Notes"
            menu
            ""
            ReadFileAsString release.ReleaseNotesFile ]
        |> ReplaceFile ("docs/content" </> (release.ReleaseNotesFile |> replace "RELEASENOTES" "ReleaseNotes"))

let buildDocumentationTarget fsiargs target =
    trace (sprintf "Building documentation (%s), this could take some time, please wait..." target)
    let fakePath = "packages" </> "build" </> "FAKE" </> "tools" </> "FAKE.exe"
    let fakeStartInfo script workingDirectory args fsiargs environmentVars =
        (fun (info: System.Diagnostics.ProcessStartInfo) ->
            info.FileName <- System.IO.Path.GetFullPath fakePath
            info.Arguments <- sprintf "%s --fsiargs -d:FAKE %s \"%s\"" args fsiargs script
            info.WorkingDirectory <- workingDirectory
            let setVar k v =
                info.EnvironmentVariables.[k] <- v
            for (k, v) in environmentVars do
                setVar k v
            setVar "MSBuild" msBuildExe
            setVar "GIT" Git.CommandHelper.gitPath
            setVar "FSI" fsiPath)
    let executeFAKEWithOutput workingDirectory script fsiargs envArgs =
        let exitCode =
            ExecProcessWithLambdas
                (fakeStartInfo script workingDirectory "" fsiargs envArgs)
                TimeSpan.MaxValue false ignore ignore
        System.Threading.Thread.Sleep 1000
        exitCode
    let exit = executeFAKEWithOutput "docs/tools" "build-docs.fsx" fsiargs ["target", target]
    if exit <> 0 then
        failwith "Generating documentation failed"
    ()

let generateDocs fail local =
    let args = if local then "" else "--define:RELEASE"
    try
        buildDocumentationTarget args "Default"
        traceImportant "Documentation generated"
    with
    | e when not fail ->
        failwith "Generating documentation failed"


// --------------------------------------------------------------------------------------
// Publishing
// Requires permissions; intended only for maintainers
// --------------------------------------------------------------------------------------

let publishReleaseTag title prefix (release:Release) =
    // inspired by Deedle/tpetricek
    let tagName = prefix + "v" + release.PackageVersion
    let tagMessage = String.concat Environment.NewLine [title + " v" + release.PackageVersion; ""; release.ReleaseNotes ]
    let cmd = sprintf """tag -a %s -m "%s" """ tagName tagMessage
    Git.CommandHelper.runSimpleGitCommand "." cmd |> printfn "%s"
    let _, remotes, _ = Git.CommandHelper.runGitCommand "." "remote -v"
    let main = remotes |> Seq.find (fun s -> s.Contains("(push)") && s.Contains("mathnet/mathnet-numerics"))
    let remoteName = main.Split('\t').[0]
    Git.Branches.pushTag "." remoteName tagName

let publishNuGet packageFiles =
    // TODO: Migrate to NuGet helper once it supports direct (non-integrated) operations
    let rec impl trials file =
        trace ("NuGet Push: " + System.IO.Path.GetFileName(file) + ".")
        try
            let args = sprintf """push "%s" -Source https://api.nuget.org/v3/index.json -T 900""" (FullName file)
            let result =
                ExecProcess (fun info ->
                    info.FileName <- "packages/build/NuGet.CommandLine/tools/NuGet.exe"
                    info.WorkingDirectory <- FullName "obj/NuGet"
                    info.Arguments <- args) (TimeSpan.FromMinutes 10.)
            if result <> 0 then failwith "Error during NuGet push."
        with exn ->
            if trials > 0 then impl (trials-1) file
            else ()
    Seq.iter (impl 3) packageFiles

let publishMirrors () =
    let repo = "../mirror-numerics"
    Git.CommandHelper.runSimpleGitCommand repo "remote update" |> printfn "%s"
    Git.CommandHelper.runSimpleGitCommand repo "push mirrors" |> printfn "%s"

let publishDocs (release:Release) =
    let repo = "../mathnet-websites"
    Git.Branches.pull repo "origin" "master"
    CopyRecursive "out/docs" "../mathnet-websites/numerics" true |> printfn "%A"
    Git.Staging.StageAll repo
    Git.Commit.Commit repo (sprintf "Numerics: %s docs update" release.PackageVersion)
    Git.Branches.pushBranch repo "origin" "master"

let publishApi (release:Release) =
    let repo = "../mathnet-websites"
    Git.Branches.pull repo "origin" "master"
    CleanDir "../mathnet-websites/numerics/api"
    CopyRecursive "out/api" "../mathnet-websites/numerics/api" true |> printfn "%A"
    Git.Staging.StageAll repo
    Git.Commit.Commit repo (sprintf "Numerics: %s api update" release.PackageVersion)
    Git.Branches.pushBranch repo "origin" "master"

let publishNuGetToArchive archivePath nupkgFile (pack:Package) =
    let tempDir = Path.GetTempPath() </> Path.GetRandomFileName()
    let archiveDir = archivePath </> pack.Id </> pack.Release.PackageVersion
    CleanDirs [tempDir; archiveDir]
    nupkgFile |> CopyFile archiveDir
    use sha512 = System.Security.Cryptography.SHA512.Create()
    let hash = File.ReadAllBytes nupkgFile |> sha512.ComputeHash |> Convert.ToBase64String
    File.WriteAllText ((archiveDir </> (Path.GetFileName(nupkgFile) + ".sha512")), hash)
    ZipHelper.Unzip tempDir nupkgFile
    !! (tempDir </> "*.nuspec") |> Copy archiveDir
    DeleteDir tempDir

let publishArchive zipOutPath nugetOutPath (bundles:Bundle list) =
    let archivePath = (environVarOrFail "MathNetReleaseArchive") </> "Math.NET Numerics"
    if directoryExists archivePath |> not then failwith "Release archive directory does not exists. Safety Check failed."
    for bundle in bundles do
        let zipFile = zipOutPath </> sprintf "%s-%s.zip" bundle.Id bundle.Release.PackageVersion
        if FileSystemHelper.fileExists zipFile then
            zipFile |> CopyFile (archivePath </> "Zip")
    for package in bundles |> List.collect (fun b -> b.Packages) do
        let nupkgFile = nugetOutPath </> sprintf "%s.%s.nupkg" package.Id package.Release.PackageVersion
        if FileSystemHelper.fileExists nupkgFile then
            trace nupkgFile
            publishNuGetToArchive (archivePath </> "NuGet") nupkgFile package
        let symbolsFile = nugetOutPath </> sprintf "%s.%s.symbols.nupkg" package.Id package.Release.PackageVersion
        if FileSystemHelper.fileExists symbolsFile then
            symbolsFile |> CopyFile (archivePath </> "Symbols")
