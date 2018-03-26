#!/usr/bin/env/python

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)

from optparse import OptionParser
import sys
import glob
from math import sqrt

class Sample :
    def __init__(self, input_file = "", name = "", displayname = "") :
        self.file = input_file
        self.name = name
        self.displayname = displayname
        self.dsid = dsid_from_file(input_file)
        self.nominal = False
        if "nom" in name :
            self.nominal = True

        self.tree = self.get_chain()

        self.counts_dict = {} # { region_name : { weight_number : [ yield, stat_err ] } }
        self.tf_dict = {} # { weight_number : [ tf, tf_err ] }

    def get_chain(self) :
        chain = r.TChain("truth")
        chain.AddFile(self.file)
        return chain

    def entries(self) :
        return self.tree.GetEntries()

class Region :
    def __init__(self, name = "", displayname = "") :
        self.name = name
        self.displayname = displayname
        self.tcut = ""

def dsid_from_file(filename) :

    dsid = filename[filename.find("_truth_") + len("_truth_") : filename.find("_truth_") + len("_truth_") + 6]
    return dsid

def get_name(dsid = "") :

    names = {}
    names["410009"] = ["ttbar_nom", "t#bar{t} nominal"]
    names["410001"] = ["ttbar_radHi", "t#bar{t} radHi"]
    names["410002"] = ["ttbar_radLo", "t#bar{t} radLo"]
    names["410003"] = ["ttbar_aMCatNLOHpp", "t#bar{t} aMCatNLO+Hpp"]
    names["410004"] = ["ttbar_PowhegHpp", "t#bar{t} Powheg+Hpp"]
    names["410189"] = ["ttbar_sherpa_dilep", "t#bar{t} Sherpa (dilep)"]

    names["410015"] = ["wt_nom", "Wt nominal"]
    names["410145"] = ["wt_PowhegHpp", "Wt Powheg+Hpp"]
    names["410164"] = ["wt_aMCatNLOHpp", "Wt aMCatNLO+Hpp"]
    names["410099"] = ["wt_radLo", "Wt radLo"]
    names["410100"] = ["wt_radHi", "Wt radHi"]

    if dsid not in names :
        print "get_name    dsid %s not found" % dsid
        sys.exit()

    return names[dsid]

def get_region(region_name = "") :

    regions = {}

    # CRttbar 
    r = Region("cr_crtt", "CR-tt")
    r.tcut = "nBJets==2 && (mbb>100 && mbb<140) && (mt2_llbb>100 && mt2_llbb<140) && (dRll>1.5 && dRll<3.0) && (ht2ratio>0.4 && ht2ratio<0.6)"
    regions["cr_crtt"] = r

    # CRWt
    r = Region("cr_wt", "CR-wt")
    r.tcut = "nBJets==2 && (mbb>140) && (mt2_bb>150) && (ht2ratio>0.6 && ht2ratio<0.8)"
    regions["cr_crwt"] = r

    # hhNonRes
    r = Region("sr_hhNonRes", "hhNonRes")
    r.tcut = "nBJets==2 && (mbb>100 && mbb<140) && (mt2_llbb>100 && mt2_llbb<140) && (ht2ratio>0.8) && dRll<0.9" # && (mt2_bb>100)"#(dRll<0.9) && (ht2ratio>0.9)"# && (mt2_bb>150)"
    #r.tcut = "nBJets==2 && (mbb>100 && mbb<140) && (mt2_llbb>100 && mt2_llbb<140) && (dRll<0.9) && (ht2ratio>0.9)"# && (mt2_bb>150)"
    regions["sr_hhNonRes"] = r

    if region_name not in regions :
        print "get_region    region %s not found" % region_name
        return None

    return regions[region_name]

def get_cross_section(dsid = "") :

    xsec = {}

    xsec["410009"] = 696.11*1.1949*0.543
    xsec["410001"] = 783.73*1.0613*0.543 
    xsec["410002"] = 783.73*1.0613*0.543 
    xsec["410003"] = 694.59*1.1975*0.543 
    xsec["410004"] = 696.32*1.1926*0.543 
    xsec["410189"] = 76.333*1.1484*1.0

    xsec["410015"] = 3.5835*1.054*1.*1.
    xsec["410164"] = 7.8714*0.9437*1.*1. 
    xsec["410099"] = 34.917*1.027*1.*1. 
    xsec["410100"] = 33.407*1.073*1.*1.
    xsec["410145"] = 4.0000*0.9443*1.*1. 

    if dsid not in xsec :
        print "get_cross_section    dsid %s not found" % dsid
        sys.exit()

    return xsec[dsid]

def get_sumw(sample, weight_index = "") :

    rfile = r.TFile.Open(sample.file)
    return rfile.Get("CutflowWeighted").GetBinContent( int(weight_index) + 1 )

def get_yields(samples, region_name, weights) :

    # get the region object
    region = get_region(region_name)
    if not region :
        print "WARNING skipping %s" % region_name
        return

    # initialize the counts dict for each sample
    for s in samples :
        s.counts_dict[region.name] = {}
        for w in weights :
            w_idx = str(w)
            s.counts_dict[region.name][w_idx] = []
            s.tf_dict[w_idx] = []


    # yields
    for s in samples :
        for w in weights :
            w_idx = str(w)

            sumw = get_sumw(s, w_idx)
            xsec = get_cross_section(s.dsid)
            lumi = 35.0 * 1000.

            h = r.TH1F("h_%s_%s_%s" % (s.name, region_name, w_idx), "", 10, 0, 10)
            cut_string = "(mcEventWeights[%s] * %f * %f / %f) * (%s)" % (w_idx, xsec, lumi, sumw, region.tcut)
            s.tree.Draw("nLeptons>>%s" % h.GetName(), cut_string, "goff")

            err = r.Double(0.0)
            sample_yield = h.IntegralAndError(0,-1, err)

            print "process %s   dsid %s  region %s   weight %s  :  %.2f +/- %.2f" % (s.name, s.dsid, region_name, w_idx, sample_yield, err)
            s.counts_dict[region.name][w_idx] = [sample_yield, err]

def calculate_transfer_factors(samples, weights) :
    print 70 * "-"
    print " calculate_transfer_factors"

    for s in samples :
        s.transfer_factors = {} # { 

    cr_name = ""
    sr_name = ""

    for w_idx in weights :
        print " > weight %s" % w_idx

        for s in samples :
            cr_counts = []
            sr_counts = []
            for region_name in s.counts_dict :
                region_yield = s.counts_dict[region_name][w_idx][0]
                region_yield_err = s.counts_dict[region_name][w_idx][1]
                if "cr_" in region_name :
                    cr_counts += [region_yield, region_yield_err]
                    if cr_name == "" :
                        cr_name = region_name
                elif "sr_" in region_name :
                    sr_counts += [region_yield, region_yield_err]
                    if sr_name == "" :
                        sr_name = region_name

            cr_yield = cr_counts[0]  
            cr_err = cr_counts[1]
            sr_yield = sr_counts[0]
            sr_err = sr_counts[1]

            tf = float(sr_yield) / float(cr_yield)
            tf_err = tf * sqrt( ( sr_err / sr_yield ) ** 2 + ( cr_err / cr_yield ) ** 2 ) 

            s.tf_dict[w_idx] = [tf, tf_err]

    print 75 * "-"
    print "transfer factors for (SR/CR) = (%s/%s)" % (sr_name, cr_name)
    for w_idx in weights :
        print 35 * '- '
        print " > weight %s" % w_idx
        nom_sample = None
        for s in samples :
            if s.nominal :
                nom_sample = s
        for s in samples :
            numerator = abs(s.tf_dict[w_idx][0] - nom_sample.tf_dict[w_idx][0])
            denominator = nom_sample.tf_dict[w_idx][0]
            tf_uncertainty = float(numerator) / float(denominator)
            
            print "   TF %s : %.4f +/- %.4f --> TF uncertainty : %.4f" % ( s.name, s.tf_dict[w_idx][0], s.tf_dict[w_idx][1], tf_uncertainty )
        
        

def main() :

    parser = OptionParser()
    parser.add_option("-i", "--indir",  default = "",  help = "directory where input files are located")
    parser.add_option("-d", "--dsid",   default = "",  help = "input dsid(s) (comma separated list)")
    parser.add_option("-r", "--region", default = "",  help = "region to get yields for (can be comma separated list)")
    parser.add_option("-w", "--weights", default = "0", help = "provide weights to use (comma separated) [default: nominal]")
    (options, args) = parser.parse_args()

    if options.indir == "" :
        print "did not provide input directory"
        sys.exit()
    if options.dsid == "" :
        print "did not provide dsid(s)"
        sys.exit()
    if options.region == "" :
        print "did not provide region(s)"
        sys.exit()

    files = glob.glob("%s/*truth*.root" % options.indir) 
    regions = options.region.split(",")
    dsids = options.dsid.split(",")
    weights = options.weights.split(",")

    if len(regions)>2 :
        print "cannot provide more than 2 regions"
        sys.exit()

    if len(regions)==2 :
        sr_found = False
        cr_found = False
        for reg in regions :
            if "cr_" in reg :
                cr_found = True
            if "sr_" in reg :
                sr_found = True
        sr_cr_found = (sr_found and cr_found)
        if not sr_cr_found :
            print "you provided 2 regions but not in the combination of a CR+SR"

    samples = []
    for dsid in dsids :
        input_file = ""
        for f in files :
            if str(dsid) in f :
                input_file = f
                break
        if input_file == "" :
            print "did not find file for dsid %s in directory %s" % (dsid,options.indir)
            sys.exit()
        s = Sample(input_file, get_name(dsid)[0], get_name(dsid)[1])
        samples.append(s)

    found_nom = False
    for s in samples :
        if "nom" in s.name :
            found_nom = True

    # if we do not compare to nominal, just assume that the first sample is the one to compare to
    # (this will be the one in the denominator of the TF calculation, so it can affect the TF
    # uncertainty!!)
    if not found_nom :
        samples[0].nominal = True
        for isample, sample in enumerate(samples) :
            if isample == 0 : continue
            sample.nominal = False

    print "Loaded %d samples" % len(samples)

    for region_name in regions :
        print 50 * '-'
        print "getting yields for %s" % region_name

        get_yields(samples, region_name, weights)

    if len(regions)==2 :
        calculate_transfer_factors(samples, weights)

    

if __name__ == '__main__' :
    main()
