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
open MathNet.Numerics.Integration

type ZoneInfo = {
    XRange : (float * float);
    HCoeff : (float * float);
    HArea : float;
    Prob : float
}

type AdaptiveRejectionMetropolisSampler =
    val m_pdfLn : float -> float
    val m_xMin : float
    val m_xMax : float
    val m_x1 : float
    val m_xn : float
    val mutable abscissas : List<float>
    val mutable proposalInfos : List<ZoneInfo>
    val m_Sampler : AbstractRandomNumberGenerator

    member private this.lineab (a:float, b:float) = (fun x -> a + b * x)

    static member private calcMoment (pdfLn: (float -> float)) (xMin: float) (xMax: float) (mfunc: (float -> float)) =
        let denom = Integrate.OnClosedInterval((fun x -> pdfLn x |> exp), xMin, xMax)
        Integrate.OnClosedInterval((fun x -> (mfunc x) * (pdfLn x |> exp)), xMin, xMax) / denom

    member private this.h x =
        let proposalInfosArray = this.proposalInfos |> List.toArray

        let rec binSearch onex lowind upind =
            let cind = (lowind + upind) / 2.0
            let oneProposalInfo = proposalInfosArray.[(round cind) |> int]
            if (onex >= (fst oneProposalInfo.XRange) && onex < (snd oneProposalInfo.XRange)) ||
                onex = (snd oneProposalInfo.XRange) && (snd oneProposalInfo.XRange) = this.m_xMax
                then onex |> (this.lineab oneProposalInfo.HCoeff)
            elif onex < (fst oneProposalInfo.XRange) then binSearch onex lowind cind
            else // onex >= (snd zoneInfo.XRange)
                binSearch onex cind upind

        binSearch x 0.0 ((this.proposalInfos.Length |> float) - 1.0)

    member private this.calcProposal absc =
        let xranges = List.zip (this.m_xMin :: absc) (absc @ [this.m_xMax])

        let secantab (xl, xu) =
            let (yl, yu) = (this.m_pdfLn xl, this.m_pdfLn xu)
            (yl - xl * (yu - yl) / (xu - xl), (yu - yl) / (xu - xl))

        let rec getEnvelope acc ranges ls revnewranges revhs =
            match acc with
            | x when x = xranges.Length -> (List.rev revnewranges), (List.rev revhs)
            | x when x = xranges.Length - 1 ->
                let rightmostab =
                    if (this.lineab (List.head ls)) (List.head ranges |> snd) < this.m_pdfLn (List.head ranges |> snd) then ls.[1]
                    else List.head ls
                getEnvelope (x+1) (List.tail ranges) (List.tail ls) ((List.head ranges) :: revnewranges) (rightmostab :: revhs)
            | x when x = xranges.Length - 2 ->
                getEnvelope (x+1) (List.tail ranges) (List.tail ls) ((List.head ranges) :: revnewranges) ((List.head ls) :: revhs)
            | 0 -> 
                let leftmostab =
                    if (this.lineab ls.[1]) (List.head ranges |> fst) < this.m_pdfLn (List.head ranges |> fst) then ls.[0]
                    else ls.[1]
                getEnvelope 1 (List.tail ranges) ls ((List.head ranges) :: revnewranges) (leftmostab :: revhs)
            | 1 -> getEnvelope 2 (List.tail ranges) (List.tail ls) ((List.head ranges) :: revnewranges) (ls.[2] :: revhs)
            | _ -> let ((la, lb), curab, (ua, ub), (xl, xu)) = ((List.head ls), ls.[1], ls.[2], (List.head ranges))
                   let xint = (ua - la) / (lb - ub)
                   if xint > xl && xint < xu && (this.lineab curab) xint < (this.lineab (la, lb)) xint
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
        
    member private this.sampleProposal () =
        let ubetween = this.m_Sampler.NextDouble()

        let rec searchrange uval (uprobs: float array) i =
            if uprobs.Length <= 1 then i
            else let upcenter = (uprobs.Length / 2)
                 let loweruprobs = Array.sub uprobs 0 upcenter
                 let upperuprobs = Array.sub uprobs upcenter (uprobs.Length - loweruprobs.Length)
                 let lowerupsum = loweruprobs |> Array.sum
                 if uval < lowerupsum then searchrange uval loweruprobs (i + 0)
                 else searchrange (uval-lowerupsum) upperuprobs (i + upcenter)
            
        let probs = List.map (fun x -> x.Prob) this.proposalInfos |> List.toArray
        let targetPropInfo = this.proposalInfos.[(searchrange ubetween probs 0)]

        let (range, coef, area) = (targetPropInfo.XRange, targetPropInfo.HCoeff, targetPropInfo.HArea) 
        let uwithin = this.m_Sampler.NextDouble()
        let res = System.Math.Log( (uwithin * (snd coef) * area / System.Math.Exp(fst coef)) + System.Math.Exp((snd coef) * (fst range)) ) / (snd coef)
        res

    member private this.onedraw xcur = 
        let rec draw cur =         
            // step 1
            let xcand = this.sampleProposal ()
            let hxcand = this.h xcand
        
            // step 2
            let uar = this.m_Sampler.NextDouble()

            // step 3
            if (log uar) > (this.m_pdfLn xcand) - hxcand then
                this.abscissas <- (xcand :: this.abscissas) |> List.sort
                this.proposalInfos <- this.calcProposal this.abscissas
                draw cur
            else
                // step 4
                let umh = this.m_Sampler.NextDouble()

                // step 5 
                let hxcur = this.h cur 

                if (log umh) > min 0.0 ( ((this.m_pdfLn xcand) + (min (this.m_pdfLn cur) hxcur)) - ((this.m_pdfLn cur) + (min (this.m_pdfLn xcand) hxcand)) ) then 
                    cur        
                else xcand

        draw xcur

    member this.Sample(x0:float, iteration:int) =
        // step 0
        this.proposalInfos <- this.calcProposal this.abscissas

        let rec draw xlst acc iter =
            if acc >= iter then
                xlst
            else 
                let onex = this.onedraw (List.head xlst)
                draw (onex :: xlst) (acc+1) iter

        let res = draw [x0] 0 iteration
        //new System.Collections.Generic.List<float>(res)
        res |> List.rev |> List.tail
        
    member this.Sample(iteration: int) =
        let mean = AdaptiveRejectionMetropolisSampler.calcMoment this.m_pdfLn this.m_xMin this.m_xMax (fun x -> x)
        this.Sample(mean, iteration)

    new(pdfLn:(float -> float), xMin:float, xMax:float, x1:float, xn:float) as this =
        //let fspdfln = (fun x -> pdfLn.Invoke(x))
        { m_pdfLn = pdfLn; m_xMin = xMin; m_xMax = xMax; m_x1 = x1; m_xn = xn; abscissas = List.empty; proposalInfos = List.empty; m_Sampler = new MersenneTwister()}
        then
            this.abscissas <- [x1; (x1 + xn) * 0.5; xn]
    
    new(pdfLn:(float -> float), xMin:float, xMax:float, x1:float, xn:float, burnIn:int) as this =
        //let fspdfln = (fun x -> pdfLn.Invoke(x))
        { m_pdfLn = pdfLn; m_xMin = xMin; m_xMax = xMax; m_x1 = x1; m_xn = xn; abscissas = List.empty; proposalInfos = List.empty; m_Sampler = new MersenneTwister()}
        then
            this.abscissas <- [x1; (x1 + xn) * 0.5; xn]
            do this.Sample(burnIn) |> ignore
                
    new(pdfLn:(float -> float), xMin: float, xMax: float) as this =
        //let fspdfln = (fun x -> pdfLn.Invoke(x))
        let mean = AdaptiveRejectionMetropolisSampler.calcMoment pdfLn xMin xMax (fun y -> y)
        let sd = AdaptiveRejectionMetropolisSampler.calcMoment pdfLn xMin xMax (fun x -> (x - mean)**2.0) |> sqrt
        let x1t = (max (mean - 2.0*sd) (0.5*(xMin + mean)))
        let xnt = (min (mean + 2.0*sd) (0.5*(mean + xMax)))
        { m_pdfLn = pdfLn; m_xMin = xMin; m_xMax = xMax; m_x1 = x1t; m_xn = xnt; abscissas = List.empty; proposalInfos = List.empty; m_Sampler = new MersenneTwister()}  
        then
            this.abscissas <- [this.m_x1; (this.m_x1 + this.m_xn) * 0.5; this.m_xn]

    new(pdfLn:(float -> float), xMin: float, xMax: float, burnIn: int) as this =
        //let fspdfln = (fun x -> pdfLn.Invoke(x))
        let mean = AdaptiveRejectionMetropolisSampler.calcMoment pdfLn xMin xMax (fun x -> x)
        let sd = AdaptiveRejectionMetropolisSampler.calcMoment pdfLn xMin xMax (fun x -> (x - mean)**2.0) |> sqrt
        let x1t = (max (mean - 2.0*sd) (0.5*(xMin + mean)))
        let xnt = (min (mean + 2.0*sd) (0.5*(mean + xMax)))
        { m_pdfLn = pdfLn; m_xMin = xMin; m_xMax = xMax; m_x1 = x1t; m_xn = xnt; abscissas = List.empty; proposalInfos = List.empty; m_Sampler = new MersenneTwister()}  
        then
            this.abscissas <- [this.m_x1; (this.m_x1 + this.m_xn) * 0.5; this.m_xn]
            do this.Sample(burnIn) |> ignore
