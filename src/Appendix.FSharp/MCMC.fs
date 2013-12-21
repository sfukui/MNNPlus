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
open MathNet.Numerics.Interpolation
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double
open MathNet.Numerics.Random
open MathNet.Numerics.Integration

[<CompiledName "ZoneInfoFSharp">]
type ZoneInfo = {
    XRange : (float * float);
    HCoeff : (float * float);
    HArea : float;
    Prob : float
}

[<CompiledName "AdaptiveRejectionMetropolisSamplerFSharp">]
type AdaptiveRejectionMetropolisSampler =
    val m_LnPdf : float -> float
    val m_XMin : float
    val m_XMax : float
    val m_X1 : float
    val m_Xn : float
    val mutable Abscissas : List<float>
    val mutable ProposalInfos : List<ZoneInfo>
    val m_Sampler : AbstractRandomNumberGenerator

    static member private calcMoment (pdfLn: (float -> float)) (xMin: float) (xMax: float) (mFunc: (float -> float)) =
        let denom = Integrate.OnClosedInterval((fun x -> pdfLn x |> exp), xMin, xMax)
        Integrate.OnClosedInterval((fun x -> (mFunc x) * (pdfLn x |> exp)), xMin, xMax) / denom

    member private this.lineAB (a:float, b:float) = (fun x -> a + b * x)

    member private this.h x =
        let proposalInfosArray = this.ProposalInfos |> List.toArray

        let rec binSearch oneX lowInd upInd =
            let cInd = (lowInd + upInd) / 2.0
            let oneProposalInfo = proposalInfosArray.[(round cInd) |> int]
            if (oneX >= (fst oneProposalInfo.XRange) && oneX < (snd oneProposalInfo.XRange)) ||
                oneX = (snd oneProposalInfo.XRange) && (snd oneProposalInfo.XRange) = this.m_XMax
                then oneX |> (this.lineAB oneProposalInfo.HCoeff)
            elif oneX < (fst oneProposalInfo.XRange) then binSearch oneX lowInd cInd
            else // onex >= (snd zoneInfo.XRange)
                binSearch oneX cInd upInd

        binSearch x 0.0 ((this.ProposalInfos.Length |> float) - 1.0)

    member private this.calcProposal absc =
        let xRanges = List.zip (this.m_XMin :: absc) (absc @ [this.m_XMax])

        let secantAB (xl, xu) =
            let (yl, yu) = (this.m_LnPdf xl, this.m_LnPdf xu)
            (yl - xl * (yu - yl) / (xu - xl), (yu - yl) / (xu - xl))

        let rec getEnvelope acc ranges ls revNewRanges revHs =
            match acc with
            | x when x = xRanges.Length -> (List.rev revNewRanges), (List.rev revHs)
            | x when x = xRanges.Length - 1 ->
                let rightmostab =
                    if (this.lineAB (List.head ls)) (List.head ranges |> snd) < this.m_LnPdf (List.head ranges |> snd) then ls.[1]
                    else List.head ls
                getEnvelope (x+1) (List.tail ranges) (List.tail ls) ((List.head ranges) :: revNewRanges) (rightmostab :: revHs)
            | x when x = xRanges.Length - 2 ->
                getEnvelope (x+1) (List.tail ranges) (List.tail ls) ((List.head ranges) :: revNewRanges) ((List.head ls) :: revHs)
            | 0 -> 
                let leftmostab =
                    if (this.lineAB ls.[1]) (List.head ranges |> fst) < this.m_LnPdf (List.head ranges |> fst) then ls.[0]
                    else ls.[1]
                getEnvelope 1 (List.tail ranges) ls ((List.head ranges) :: revNewRanges) (leftmostab :: revHs)
            | 1 -> getEnvelope 2 (List.tail ranges) (List.tail ls) ((List.head ranges) :: revNewRanges) (ls.[2] :: revHs)
            | _ -> let ((la, lb), curab, (ua, ub), (xl, xu)) = ((List.head ls), ls.[1], ls.[2], (List.head ranges))
                   let xint = (ua - la) / (lb - ub)
                   if xint > xl && xint < xu && (this.lineAB curab) xint < (this.lineAB (la, lb)) xint
                       then getEnvelope (acc+1) (List.tail ranges) (List.tail ls)
                                ([(xint, xu); (xl, xint)] @ revNewRanges) ([(ua, ub); (la, lb)] @ revHs)
                   else getEnvelope (acc+1) (List.tail ranges) (List.tail ls) ((xl,xu) :: revNewRanges) (curab :: revHs)
        
        let secants = List.mapi (fun i (l, u) ->
                                if l <> (-infinity) && u <> infinity then secantAB (l, u)
                                else (nan, nan) ) xRanges

        let (newXRanges, hs) = getEnvelope 0 xRanges secants [] []

        let areas = List.map2 (fun (l,u) (a,b) -> (System.Math.Exp(a + b * u) - System.Math.Exp(a + b * l)) / b) newXRanges hs
        let probs = List.map (fun x -> x / (List.sum areas)) areas
        
        List.map3 (fun range hi (area, prob) -> { XRange = range; HCoeff = hi; HArea = area; Prob = prob }) 
            newXRanges hs (List.zip areas probs)
        
    member private this.sampleProposal () =
        let uBetween = this.m_Sampler.NextDouble()

        let rec searchRange uVal (uProbs: float array) i =
            if uProbs.Length <= 1 then i
            else let uPCenter = (uProbs.Length / 2)
                 let lowerUProbs = Array.sub uProbs 0 uPCenter
                 let upperUProbs = Array.sub uProbs uPCenter (uProbs.Length - lowerUProbs.Length)
                 let lowerUPSum = lowerUProbs |> Array.sum
                 if uVal < lowerUPSum then searchRange uVal lowerUProbs (i + 0)
                 else searchRange (uVal-lowerUPSum) upperUProbs (i + uPCenter)
            
        let probs = List.map (fun x -> x.Prob) this.ProposalInfos |> List.toArray
        let targetPropInfo = this.ProposalInfos.[(searchRange uBetween probs 0)]

        let (range, coef, area) = (targetPropInfo.XRange, targetPropInfo.HCoeff, targetPropInfo.HArea) 
        let uWithin = this.m_Sampler.NextDouble()
        let res = System.Math.Log( (uWithin * (snd coef) * area / System.Math.Exp(fst coef)) + System.Math.Exp((snd coef) * (fst range)) ) / (snd coef)
        res

    member private this.oneDraw xcur = 
        let rec draw cur =         
            // step 1
            let xCand = this.sampleProposal ()
            let hXCand = this.h xCand
        
            // step 2
            let uAR = this.m_Sampler.NextDouble()

            // step 3
            if (log uAR) > (this.m_LnPdf xCand) - hXCand then
                this.Abscissas <- (xCand :: this.Abscissas) |> List.sort
                this.ProposalInfos <- this.calcProposal this.Abscissas
                draw cur
            else
                // step 4
                let uMH = this.m_Sampler.NextDouble()

                // step 5 
                let hXCur = this.h cur 

                if (log uMH) > min 0.0 ( ((this.m_LnPdf xCand) + (min (this.m_LnPdf cur) hXCur)) - ((this.m_LnPdf cur) + (min (this.m_LnPdf xCand) hXCand)) ) then 
                    cur        
                else xCand

        draw xcur

    member this.Sample(x0:float, iteration:int) =
        // step 0
        this.ProposalInfos <- this.calcProposal this.Abscissas

        let rec draw xlst acc iter =
            if acc >= iter then
                xlst
            else 
                let oneX = this.oneDraw (List.head xlst)
                draw (oneX :: xlst) (acc+1) iter

        let res = draw [x0] 0 iteration
        res |> List.rev |> List.tail
        
    member this.Sample(iteration: int) =
        let mean = AdaptiveRejectionMetropolisSampler.calcMoment this.m_LnPdf this.m_XMin this.m_XMax (fun x -> x)
        this.Sample(mean, iteration)

    new(lnPdf:(float -> float), xMin:float, xMax:float, x1:float, xn:float) as this =
        { m_LnPdf = lnPdf; m_XMin = xMin; m_XMax = xMax; m_X1 = x1; m_Xn = xn; Abscissas = List.empty; ProposalInfos = List.empty; m_Sampler = new MersenneTwister()}
        then
            this.Abscissas <- [x1; (x1 + xn) * 0.5; xn]
        
    new(lnPdf:(float -> float), xMin:float, xMax:float, x1:float, xn:float, seed:int) as this =
        { m_LnPdf = lnPdf; m_XMin = xMin; m_XMax = xMax; m_X1 = x1; m_Xn = xn; Abscissas = List.empty; ProposalInfos = List.empty; m_Sampler = new MersenneTwister(seed)}
        then
            this.Abscissas <- [x1; (x1 + xn) * 0.5; xn]

    new(lnPdf:(float -> float), xMin:float, xMax:float, x1:float, xn:float, sampler:AbstractRandomNumberGenerator) as this =
        { m_LnPdf = lnPdf; m_XMin = xMin; m_XMax = xMax; m_X1 = x1; m_Xn = xn; Abscissas = List.empty; ProposalInfos = List.empty; m_Sampler = sampler}
        then
            this.Abscissas <- [x1; (x1 + xn) * 0.5; xn]
                
    new(lnPdf:(float -> float), xMin: float, xMax: float) as this =
        let mean = AdaptiveRejectionMetropolisSampler.calcMoment lnPdf xMin xMax (fun y -> y)
        let sd = AdaptiveRejectionMetropolisSampler.calcMoment lnPdf xMin xMax (fun x -> (x - mean)**2.0) |> sqrt
        let x1t = (max (mean - 2.0*sd) (0.5*(xMin + mean)))
        let xnt = (min (mean + 2.0*sd) (0.5*(mean + xMax)))
        { m_LnPdf = lnPdf; m_XMin = xMin; m_XMax = xMax; m_X1 = x1t; m_Xn = xnt; Abscissas = List.empty; ProposalInfos = List.empty; m_Sampler = new MersenneTwister()}  
        then
            this.Abscissas <- [this.m_X1; (this.m_X1 + this.m_Xn) * 0.5; this.m_Xn]

    new(lnPdf:(float -> float), xMin: float, xMax: float, seed: int) as this =
        let mean = AdaptiveRejectionMetropolisSampler.calcMoment lnPdf xMin xMax (fun y -> y)
        let sd = AdaptiveRejectionMetropolisSampler.calcMoment lnPdf xMin xMax (fun x -> (x - mean)**2.0) |> sqrt
        let x1t = (max (mean - 2.0*sd) (0.5*(xMin + mean)))
        let xnt = (min (mean + 2.0*sd) (0.5*(mean + xMax)))
        { m_LnPdf = lnPdf; m_XMin = xMin; m_XMax = xMax; m_X1 = x1t; m_Xn = xnt; Abscissas = List.empty; ProposalInfos = List.empty; m_Sampler = new MersenneTwister(seed)}  
        then
            this.Abscissas <- [this.m_X1; (this.m_X1 + this.m_Xn) * 0.5; this.m_Xn]
