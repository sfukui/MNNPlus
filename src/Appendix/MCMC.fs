// Math.NET Numerics Appendix Library License (MIT/X11)
// ===================================================
// 
// Copyright (c) 2013 Fukui Shogo
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

namespace MathNet.Numerics.Statistics.Mcmc

open MathNet.Numerics
open MathNet.Numerics.Differentiation
open MathNet.Numerics.Interpolation.Algorithms
open MathNet.Numerics.LinearAlgebra.Generic
open MathNet.Numerics.LinearAlgebra.Double
open MathNet.Numerics.Random

type ZoneInfo = {
    XRange : (float * float);
    HCoeff : (float * float);
    HArea : float;
    Prob : float
}

type AdaptiveRejectionMetropolisSampler (pdfLn:(float -> float), xMin:float, xMax:float, x1:float, xn:float) =
    let mutable abscissas = List.empty
    let mutable proposalInfos = List.empty
    let mutable burnIn = 0

    do abscissas <- [x1; (x1 + xn) * 0.5; xn]

    let mtSampler = new MersenneTwister();

    let lineab (a:float, b:float) = (fun x -> a + b * x)

    let h x =
        let proposalInfosArray = proposalInfos |> List.toArray

        let rec binSearch onex lowind upind =
            let cind = (lowind + upind) / 2.0
            let oneProposalInfo = proposalInfosArray.[(round cind) |> int]
            if (onex >= (fst oneProposalInfo.XRange) && onex < (snd oneProposalInfo.XRange)) ||
                onex = (snd oneProposalInfo.XRange) && (snd oneProposalInfo.XRange) = xMax
                then onex |> (lineab oneProposalInfo.HCoeff)
            elif onex < (fst oneProposalInfo.XRange) then binSearch onex lowind cind
            else // onex >= (snd zoneInfo.XRange)
                binSearch onex cind upind

        binSearch x 0.0 ((proposalInfos.Length |> float) - 1.0)

    let calcProposal absc =
        let xranges = List.zip (xMin :: absc) (absc @ [xMax])

        let secantab (xl, xu) =
            let (yl, yu) = (pdfLn xl, pdfLn xu)
            (yl - xl * (yu - yl) / (xu - xl), (yu - yl) / (xu - xl))

        let rec getEnvelope acc ranges ls revnewranges revhs =
            match acc with
            | x when x = xranges.Length -> (List.rev revnewranges), (List.rev revhs)
            | x when x = xranges.Length - 1 ->
                let rightmostab =
                    if (lineab (List.head ls)) (List.head ranges |> snd) < pdfLn (List.head ranges |> snd) then ls.[1]
                    else List.head ls
                getEnvelope (x+1) (List.tail ranges) (List.tail ls) ((List.head ranges) :: revnewranges) (rightmostab :: revhs)
            | x when x = xranges.Length - 2 ->
                getEnvelope (x+1) (List.tail ranges) (List.tail ls) ((List.head ranges) :: revnewranges) ((List.head ls) :: revhs)
            | 0 -> 
                let leftmostab =
                    if (lineab ls.[1]) (List.head ranges |> fst) < pdfLn (List.head ranges |> fst) then ls.[0]
                    else ls.[1]
                getEnvelope 1 (List.tail ranges) ls ((List.head ranges) :: revnewranges) (leftmostab :: revhs)
            | 1 -> getEnvelope 2 (List.tail ranges) (List.tail ls) ((List.head ranges) :: revnewranges) (ls.[2] :: revhs)
            | _ -> let ((la, lb), curab, (ua, ub), (xl, xu)) = ((List.head ls), ls.[1], ls.[2], (List.head ranges))
                   let xint = (ua - la) / (lb - ub)
                   if xint > xl && xint < xu && (lineab curab) xint < (lineab (la, lb)) xint
                       then getEnvelope (acc+1) (List.tail ranges) (List.tail ls)
                                ([(xint, xu); (xl, xint)] @ revnewranges) ([(ua, ub); (la, lb)] @ revhs)
                   else getEnvelope (acc+1) (List.tail ranges) (List.tail ls) ((xl,xu) :: revnewranges) (curab :: revhs)
        
        let secants = List.mapi (fun i (l, u) ->
                                if l <> (-infinity) && u <> infinity then secantab (l, u)
                                else (nan, nan) ) xranges

        let (newxranges, hs) = getEnvelope 0 xranges secants [] []

        let areas = List.map2 (fun (l,u) (a,b) -> (System.Math.Exp(a + b * u) - System.Math.Exp(a + b * l)) / b) newxranges hs
        let probs = List.map (fun x -> x / (List.sum areas)) areas
        
        List.map3 (fun range hi (area, prob) -> { XRange = range; HCoeff = hi; HArea = area; Prob = prob }) 
            newxranges hs (List.zip areas probs)
        
    let sampleProposal () =
        let ubetween = mtSampler.NextDouble()

        let rec searchrange uval (uprobs: float array) i =
            if uprobs.Length <= 1 then i
            else let upcenter = (uprobs.Length / 2)
                 let loweruprobs = Array.sub uprobs 0 upcenter
                 let upperuprobs = Array.sub uprobs upcenter (uprobs.Length - loweruprobs.Length)
                 let lowerupsum = loweruprobs |> Array.sum
                 if uval < lowerupsum then searchrange uval loweruprobs (i + 0)
                 else searchrange (uval-lowerupsum) upperuprobs (i + upcenter)
            
        let probs = List.map (fun x -> x.Prob) proposalInfos |> List.toArray
        let targetPropInfo = proposalInfos.[(searchrange ubetween probs 0)]

        let (range, coef, area) = (targetPropInfo.XRange, targetPropInfo.HCoeff, targetPropInfo.HArea) 
        let uwithin = mtSampler.NextDouble()
        let res = System.Math.Log( (uwithin * (snd coef) * area / System.Math.Exp(fst coef)) + System.Math.Exp((snd coef) * (fst range)) ) / (snd coef)
        res

    let rec onedraw xcur =         
        // step 1
        let xcand = sampleProposal ()
        let hxcand = h xcand
        
        // step 2
        let uar = mtSampler.NextDouble()

        // step 3
        if (log uar) > (pdfLn xcand) - hxcand then
            abscissas <- (xcand :: abscissas) |> List.sort
            proposalInfos <- calcProposal abscissas
            onedraw xcur
        else
            // step 4
            let umh = mtSampler.NextDouble()

            // step 5 
            let hxcur = h xcur 

            if (log umh) > min 0.0 ( ((pdfLn xcand) + (min (pdfLn xcur) hxcur)) - ((pdfLn xcur) + (min (pdfLn xcand) hxcand)) ) then 
                xcur        
            else xcand

    member x.BurnIn
           with get() = burnIn
           and set(value) = burnIn <- value

    member x.Sample(x0:float, iteration:int) =
        // step 0
        proposalInfos <- calcProposal abscissas

        let rec draw xlst acc iter =
            if acc >= iter then xlst
            else 
                let onex = onedraw (List.head xlst)
                draw (onex :: xlst) (acc+1) iter

        draw [x0] 0 burnIn |> ignore
        let res = draw [x0] 0 iteration
        res

