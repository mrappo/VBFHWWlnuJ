import os,commands
import sys
from optparse import OptionParser
import subprocess





import glob
import math
import array
import ROOT
import ntpath




from array import array
from datetime import datetime

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite




parser = OptionParser()

parser.add_option('--channel', action="store", type="string", dest="channel", default="mu")
#parser.add_option('--sample', action="store", type="string", dest="sample", default="BulkGraviton")
parser.add_option('--ntuple', action="store", type="string", dest="ntuple", default="WWTree_22sep_jecV7_lowmass")
parser.add_option('--lumi', action="store", type="float", dest="lumi", default="2300")
parser.add_option('--nodata', action='store_true', dest='nodata', default=False)
(options, args) = parser.parse_args()
currentDir = os.getcwd();

###########################################################################################
######## GLOBAL VARIABLE DEFINITION
###########################################################################################

number_Events_type=11;
masses_BulkGraviton=[600,800,1000];
masses_Higgs=[650,1000];
       
BulkGraviton_xsec=[0.177400,0.0331548,0.008993];
VBF_BulkGraviton_xsec=[0.01089,0.00217,0.000655];
Higgs_xsec=[0.33639,0.06765];
VBF_Higgs_xsec=[0.03354,0.02375];
		    
NumEntriesBefore_BulkGraviton=[49600,50000,50000];
NumEntriesBefore_VBF_BulkGraviton=[50000,50000,50000];
NumEntriesBefore_Higgs=[399600,400000];
NumEntriesBefore_VBF_Higgs=[398400,400000];
		    
ScaleFactor_BulkGraviton=[900,2000,6000];
ScaleFactor_Higgs=[25,120]


#cuts_itemize=["issignal==1 && abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001","v_pt>200","pfMET>40","l_pt>40","ungroomed_jet_pt>200","nBTagJet_medium <1","jet_mass_pr > 65 && jet_mass_pr < 105","jet_tau2tau1 < 0.45"];
cuts_itemize=["deltaR_lak8jet>(TMath::Pi()/2.0) && abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001","TMath::Abs(deltaphi_METak8jet)>2.0","TMath::Abs(deltaphi_Vak8jet)>2.0","v_pt>200","pfMET>40","l_pt>40","ungroomed_jet_pt>200","nBTagJet_medium <1","jet_tau2tau1 < 0.45","jet_mass_pr > 65 && jet_mass_pr < 105 "];
#cuts_itemize=["issignal==1 && v_pt>200","pfMET>40 && l_pt>40 && ungroomed_jet_pt>200 && nBTagJet_medium <1 && jet_mass_pr > 65 && jet_mass_pr < 105 && abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001","jet_tau2tau1 < 0.45"];


###########################################################################################
######## FUNCTION DEFINITION
###########################################################################################

def make_Input_file(ntuple,channel,lumi,cuts_file,SampleList_filename,table,k,ScaleFactor_table,cut_number,cuttype):
    
    sample=table[0][k][0];
    mass=table[0][k][1];        
    mass_str=str("%.0f"%mass)
    reducedggName=table[0][k][3];
    reducedVBFName=table[1][k][3];
    
    ScaleFactor=ScaleFactor_table[cuttype-1][k][cut_number-1];
    
    
    if channel=="mu":
       leptonT="muon";
    elif channel=="el":
       leptonT="electron";   
    else:
       leptonT="emu";
    
        
    if options.nodata:
       withData="true";
    else:
       withData="false";
    
    cts=str("%.0f"%cuttype);
    cns=str("%.0f"%cut_number);
    data_save="output/run2/MCDATAComparisonPlot_mu_%s_%s%s_%s_%s/"%(ntuple,sample,mass,cts,cns);
       
    filename="cfg/DataMCComparison_InputCfgFile/MATTEO_DataMCComparison_InputCfgFile_%s%s_%.0f_%.0f.cfg"%(sample,mass_str,cuttype,cut_number)
    in_file=open(filename,'w+');
    in_file.write("[Input]\n\n");
    in_file.write(("InputDirectory = /afs/cern.ch/user/l/lbrianza/work/public/%s/WWTree_%s/\n")%(ntuple,channel));
    in_file.write("TreeName = otree\n");
    in_file.write(("LeptonType = %s\n")%leptonT);
    in_file.write(("InputSampleList = %s\n")%(SampleList_filename));
    in_file.write("InputVariableList = cfg/DataMCComparison_InputCfgFile/Run2_Variables_76x_data2015.txt\n");
    in_file.write(("InputCutList = %s\n")%cuts_file);
    in_file.write(("SignalqqHName = %s\n")%reducedVBFName);
    in_file.write(("SignalggHName = %s\n")%reducedggName);
    in_file.write(("WithoutData = %s\n\n\n")%withData);
    in_file.write("[Option]\n\n");
    in_file.write("BackgroundWeight = genWeight\n");
    in_file.write("BackgroundWeightMCatNLO = 1\n");
    in_file.write("SignalggHWeight = 1\n");
    in_file.write("SignalqqHWeight = 1\n");
    in_file.write("SignalGravitonWeight = genWeight\n");
    in_file.write(("Lumi = %f\n")%options.lumi);
    in_file.write("ttbarControlplots = false\n");
    in_file.write("SignalScaleFactor = %f\n"%(ScaleFactor));
    in_file.write("NormalizeSignalToData = false\n");
    in_file.write("NormalizeBackgroundToData = false\n\n\n");
    in_file.write("[Output]\n\n");
    in_file.write(("OutputRootDirectory = %s\n"%(data_save)));
    in_file.write("OutputRootFile = Run2_MCDataComparisonRSGraviton2000_mu.root\n");
    in_file.write("\n");
    in_file.close();
    
    info=[filename,data_save];
    return info;


    
    
    





def run_log(input_file_info,in_cc_dir,table,k,lumi,output_summaryFile,Events_table_log,ctype,cnumber,latex_file,cut_total_number_log):
    input_file=input_file_info[0];
    cuttype=ctype-1;
    sample=table[0][k][0];
    mass=table[0][k][1];        
    mass_str=str("%.0f"%mass)
    process_name=sample+mass_str
    lumi_str=str("%.0f"%lumi)
    Events_table_log[cuttype][cnumber][k][0]=sample;
    Events_table_log[cuttype][cnumber][k][1]=mass;
    
    latex_file.write("%s\n"%sample);
    latex_file.write("%s\n"%mass_str);
    latex_file.write("%s\n"%cuttype);
    latex_file.write("%.0f\n"%cnumber);
            
            
    log_dir=in_cc_dir+"/log";
    if not os.path.isdir(log_dir):
           pd1 = subprocess.Popen(['mkdir',log_dir]);
           pd1.wait();
            
    data_dir=in_cc_dir+"/data";
    if not os.path.isdir(data_dir):
           pd2 = subprocess.Popen(['mkdir',data_dir]);
           pd2.wait();
    
    cn_print=cnumber+1; 
    if cuttype:
       print_string=("\tSingle CUT: %.0f of %.0f")%(cn_print,cut_total_number_log);     
    else:
       print_string=("\tRecursive CUTS: %.0f of %.0f")%(cn_print,cut_total_number_log);
    
    log_file1=log_dir+"/Log_ControlPlots_%s_%s%s.txt"%(options.channel,sample,mass_str)
    log_file2=log_dir+"/Error_Log_ControlPlots_%s_%s%s.txt"%(options.channel,sample,mass_str)
     
          
    output_log1=open(log_file1,'w+')
    output_log1.write("\n\n--------------------------------\n\n")
    output_log1.write("STARTING\t\t")
    output_log1.write(process_name)
    output_log1.write(print_string)
    output_log1.write("\n\n--------------------------------\n\n")
            
    output_log2=open(log_file2,'w+')
    output_log2.write("\n\n--------------------------------\n\n")
    output_log2.write("STARTING\t\t")
    output_log2.write(process_name)
    output_log2.write(print_string)
    output_log2.write("\n\n--------------------------------\n\n")
          
    sys.stdout.write("\n\n------------------------------------------------------------------\n\n")
    sys.stdout.write("STARTING\t\t")
    sys.stdout.write(process_name)
    sys.stdout.write(print_string)
    sys.stdout.write("\n\n------------------------------------------------------------------\n\n")
    
    #output_summaryFile.write("\n\n--------------------------------\n\n")
    output_summaryFile.write("STARTING\t\t")
    output_summaryFile.write(process_name)
    output_summaryFile.write(print_string)
    output_summaryFile.write("\n\n--------------------------------\n\n")                
    
    Wjets_Pythia_Events_d="Jets Pythia  Events:           ";
    Wjets_Herwig_Events_d="Jets Herwig  Events:           ";
    TTbar_Powegh_Events_d="Tbar Powegh Events:           ";
    TTbar_MC_Events_d="Tbar mc@nlo Events:           ";
    VV_QCD_Events_d="V Events QCD:              ";
    WW_EWK_Events_d="W Events EWK:              ";
    STop_Events_d="Top Events:            ";
    All_bkg_Pythia_d="l Backgrounds Events Pythia: ";
    All_bkg_Herwig_d="l Backgrounds Events Herwig: ";
    Signal_gg_d="ignal ggH:             ";
    Signal_VBF_d="ignal qqH:             ";
    
     

    
    Events_type=[Wjets_Pythia_Events_d,Wjets_Herwig_Events_d,TTbar_Powegh_Events_d,TTbar_MC_Events_d,VV_QCD_Events_d,WW_EWK_Events_d,STop_Events_d,All_bkg_Pythia_d,All_bkg_Herwig_d,Signal_gg_d,Signal_VBF_d];
    
    
    
    p1 = subprocess.Popen(['./bin/DataMCComparisonPlot.exe',input_file],stdout=subprocess.PIPE,stderr=output_log2)
    p1.wait()
    start=0;
    for line in p1.stdout:
        #sys.stdout.write(line)
        output_log1.write(line)
        if line.find('Event Scaled To Lumi') !=-1:
           start=1;
        if start:   
           output_summaryFile.write(line);
        cn=0;
        for ev in Events_type:
            if line.find(ev) !=-1:
               cut_string = line.split(ev);
               new_string = cut_string[1]
               #print(new_string)
               print line
               val=float(new_string);
               Events_table_log[cuttype][cnumber][k][cn+2]=val;
               #print "cuttype: %.0f\t CutNumber:%.0f\t SampleNumber:%.0f\t Line:%.0f\t VALUE:%f"%(cuttype,cnumber,k,cn+2,val)
               latex_file.write("%f\n"%val);
               
            cn=cn+1;
    p1.wait();    

    sig_Pytia=Events_table_log[cuttype][cnumber][k][number_Events_type+1]/(1+TMath.Sqrt(Events_table_log[cuttype][cnumber][k][number_Events_type-2]));
    sig_Herwig=Events_table_log[cuttype][cnumber][k][number_Events_type+1]/(1+TMath.Sqrt(Events_table_log[cuttype][cnumber][k][number_Events_type-1]));
    Events_table_log[cuttype][cnumber][k][number_Events_type+2]=sig_Pytia;
    Events_table_log[cuttype][cnumber][k][number_Events_type+3]=sig_Herwig;
    print "\n\n------------Calculate Significance----------------\n\n"
    print "\nqq signal: %f\n"% Events_table_log[cuttype][cnumber][k][number_Events_type+1]
    print "\nPythia Bkg: %f\n"% Events_table_log[cuttype][cnumber][k][number_Events_type-2]
    print "\nHewig Bkg: %f\n"% Events_table_log[cuttype][cnumber][k][number_Events_type-1]
    print "\nSig Pythia: %f\n"% sig_Pytia
    print "\nSig Herwig: %f\n"% sig_Herwig    
    #print "\n\n--------------------------------------------------\n\n"
    latex_file.write("%f\n"%sig_Pytia);
    latex_file.write("%f\n"%sig_Herwig);
    output_summaryFile.write("Significance Pythia:    %f\n"%sig_Pytia);
    output_summaryFile.write("Significance Herwig:    %f\n"%sig_Herwig);
    
    
    latex_file.write("\n\n");
    
    
        

    path_dir_in_tmp=input_file_info[1];
    path_dir_in=path_dir_in_tmp+"Run2_MCDataComparisonRSGraviton2000_mu_plot"; 
    #path_dir_in="output/run2/MCDATAComparisonPlot_mu_22sep_%s/Run2_MCDataComparisonRSGraviton2000_mu_plot"%process_name;
    path_dir_out=in_cc_dir;
    p2 = subprocess.Popen(['cp','-r',path_dir_in,path_dir_out])
    p2.wait()



    data_in=in_cc_dir+"/Run2_MCDataComparisonRSGraviton2000_mu_plot";
    data_out=in_cc_dir+"/data";
    p3 = subprocess.Popen(['mv',data_in,data_out])
    p3.wait()

    root_in=path_dir_in_tmp+"Run2_MCDataComparisonRSGraviton2000_mu.root"
    root_out=data_out+"/Root_ControlPlots_out_%s.root"%process_name;
    p4 = subprocess.Popen(['cp',root_in,root_out])
    p4.wait()

    if os.path.isdir(data_in):
       p5=subprocess.Popen(['rm','-r',data_in])
       p5.wait()
                     
    
            
    output_log1.write("\n\n--------------------------------\n\n")
    output_log1.write("ENDED\t\t")
    output_log1.write(process_name)
    output_log1.write(print_string)
    output_log1.write("\n\n--------------------------------\n\n")
    output_log1.close()
            
    output_log2.write("\n\n--------------------------------\n\n")
    output_log2.write("ENDED\t\t")
    output_log2.write(process_name)
    output_log2.write(print_string)
    output_log2.write("\n\n--------------------------------\n\n")
    output_log2.close()
            
                        
    sys.stdout.write("\n\n------------------------------------------------------------------\n\n")
    sys.stdout.write("ENDED\t\t")
    sys.stdout.write(process_name)
    sys.stdout.write(print_string)
    sys.stdout.write("\n\n------------------------------------------------------------------\n\n")
    
    output_summaryFile.write("\n\n--------------------------------\n\n")
    output_summaryFile.write("ENDED\t\t")
    output_summaryFile.write(process_name)
    output_summaryFile.write(print_string)
    output_summaryFile.write("\n\n--------------------------------\n\n")
    
    
    return Events_table_log;
    
     
def make_sample_table(cuts_number):
    
    
             
   
            counter_BulkGraviton=0;
            for line in masses_BulkGraviton:
                counter_BulkGraviton=counter_BulkGraviton+1;
            
            counter_Higgs=0;
            for line in masses_Higgs:
                counter_Higgs=counter_Higgs+1;
            
            
            
            reducedName_BulkGraviton=[0 for z in range(counter_BulkGraviton)];
            reducedName_VBF_BulkGraviton=[0 for z in range(counter_BulkGraviton)];
            counter=0;
            for counter in range(counter_BulkGraviton):
                reducedName_BulkGraviton[counter]="BG"+str(masses_BulkGraviton[counter]);
                reducedName_VBF_BulkGraviton[counter]="VBF_BG"+str(masses_BulkGraviton[counter]);
                
                                
            reducedName_Higgs=[0 for k in range(counter_Higgs)];
            reducedName_VBF_Higgs=[0 for k in range(counter_Higgs)];
            counter=0;
            for counter in range(counter_Higgs):
                reducedName_Higgs[counter]="HIGGS"+str(masses_Higgs[counter]);
                reducedName_VBF_Higgs[counter]="VBF_HIGGS"+str(masses_Higgs[counter]);
            
            signalName_BulkGraviton=[0 for z in range(counter_BulkGraviton)];
            signalName_VBF_BulkGraviton=[0 for z in range(counter_BulkGraviton)];
            counter=0;
            for counter in range(counter_BulkGraviton):
                signalName_BulkGraviton[counter]="WWTree_BulkGraviton"+str(masses_BulkGraviton[counter]);
                signalName_VBF_BulkGraviton[counter]="WWTree_VBFBulkGraviton"+str(masses_BulkGraviton[counter]);
                
            signalName_Higgs=[0 for k in range(counter_Higgs)];
            signalName_VBF_Higgs=[0 for k in range(counter_Higgs)];
            counter=0;
            for counter in range(counter_Higgs):
                signalName_Higgs[counter]="WWTree_Higgs"+str(masses_Higgs[counter]);
                signalName_VBF_Higgs[counter]="WWTree_VBFHiggs"+str(masses_Higgs[counter]);
            
            
            attributes=["Signal Name","Reduced Name","Color","XSec(pb)","NumEntriesBefore"];
            number_attributes=0;
            for c in attributes:
                number_attributes=number_attributes+1;
                
            VBF_color=1;
            gg_color=13;
            
            
            number_attributes=number_attributes+2;
            total_number_samples=counter_Higgs+counter_BulkGraviton;
            

            Events_table=[[[[0 for k in range(number_Events_type+4)] for j in range(total_number_samples)] for a in range(cuts_number) ]for i in range(2)];
            ScaleFactor_table=[[[0 for k in range(cuts_number)] for j in range(total_number_samples)] for i in range(2)];
            i=j=k=0;
            for i in range(2):
                for j in range(total_number_samples):
                    for k in range(cuts_number):
                        if j<3:
                           ScaleFactor_table[i][j][k]=ScaleFactor_BulkGraviton[j];
                        else:
                           ScaleFactor_table[i][j][k]=ScaleFactor_Higgs[j-3];
                           
                           
                           
                           
                           
            xsec_table= [[[0 for j in range(number_attributes)] for i in range(total_number_samples) ]for k in range(2)  ];
            
            
            counter=0;
            for counter in range(counter_BulkGraviton):
                xsec_table[0][counter][0]="BulkGraviton";
                xsec_table[0][counter][1]=masses_BulkGraviton[counter];
                xsec_table[0][counter][2]=signalName_VBF_BulkGraviton[counter];
                xsec_table[0][counter][3]=reducedName_VBF_BulkGraviton[counter];
                xsec_table[0][counter][4]=VBF_color;
                xsec_table[0][counter][5]=VBF_BulkGraviton_xsec[counter];
                xsec_table[0][counter][6]=NumEntriesBefore_VBF_BulkGraviton[counter];
             
            counter=0;
            for counter in range(counter_Higgs):
                xsec_table[0][counter+counter_BulkGraviton][0]="Higgs";
                xsec_table[0][counter+counter_BulkGraviton][1]=masses_Higgs[counter];
                xsec_table[0][counter+counter_BulkGraviton][2]=signalName_VBF_Higgs[counter];
                xsec_table[0][counter+counter_BulkGraviton][3]=reducedName_VBF_Higgs[counter];
                xsec_table[0][counter+counter_BulkGraviton][4]=VBF_color;
                xsec_table[0][counter+counter_BulkGraviton][5]=VBF_Higgs_xsec[counter];
                xsec_table[0][counter+counter_BulkGraviton][6]=NumEntriesBefore_VBF_Higgs[counter];
                
            counter=0;    
            for counter in range(counter_BulkGraviton):
                xsec_table[1][counter][0]="BulkGraviton";
                xsec_table[1][counter][1]=masses_BulkGraviton[counter];
                xsec_table[1][counter][2]=signalName_BulkGraviton[counter];
                xsec_table[1][counter][3]=reducedName_BulkGraviton[counter];
                xsec_table[1][counter][4]=gg_color;
                xsec_table[1][counter][5]=BulkGraviton_xsec[counter];
                xsec_table[1][counter][6]=NumEntriesBefore_BulkGraviton[counter];
            
            counter=0;
            for counter in range(counter_Higgs):
                xsec_table[1][counter+counter_BulkGraviton][0]="Higgs";
                xsec_table[1][counter+counter_BulkGraviton][1]=masses_Higgs[counter];
                xsec_table[1][counter+counter_BulkGraviton][2]=signalName_Higgs[counter];
                xsec_table[1][counter+counter_BulkGraviton][3]=reducedName_Higgs[counter];
                xsec_table[1][counter+counter_BulkGraviton][4]=gg_color;
                xsec_table[1][counter+counter_BulkGraviton][5]=Higgs_xsec[counter];
                xsec_table[1][counter+counter_BulkGraviton][6]=NumEntriesBefore_Higgs[counter];
            
            value=[xsec_table,counter_BulkGraviton,counter_Higgs,number_attributes,ScaleFactor_table,Events_table]
            return value;
            



def make_SampleList_file(ntuple,sample_number,sample_value):

    sample=sample_value[0][sample_number][0];
    mass=str(sample_value[0][sample_number][1]);
        
      
    filename="cfg/DataMCComparison_InputCfgFile/MATTEO_SampleList_76x_%s_%s%s.txt"%(ntuple,sample,mass);
    
    
    in_file=open(filename,'w+');
    in_file.write("###########################################################################################################################################\n");
    in_file.write("# Sample Name												Reduced Name	Color				XSec (pb)		NumEntriesBefore\n");
    in_file.write("###########################################################################################################################################\n");
    in_file.write("WWTree_data_golden_2p1										DATA			1					1				1\n");
    in_file.write("#WWTree_WJets												W+Jets			2					61526.7			24156124\n");
    in_file.write("WWTree_WJets100												W+Jets          2                  	1347            10205377\n");
    in_file.write("WWTree_WJets200												W+Jets          2                  	360.           	4949568\n");
    in_file.write("WWTree_WJets400												W+Jets          2                  	48.9            1943664\n");
    in_file.write("#WWTree_WJets600												W+Jets          2                  	18.77           1041358\n");
    in_file.write("WWTree_WJets600bis											W+Jets          2                  	12.8            3767766\n");
    in_file.write("WWTree_WJets800												W+Jets          2                  	5.26            1568277\n");
    in_file.write("WWTree_WJets1200												W+Jets          2                  	1.33            246239\n");
    in_file.write("WWTree_WJets2500												W+Jets          2                  	0.03089         251982\n");
    in_file.write("#WWTree_WW													WW              4                  	118.7           988418\n");
    in_file.write("#WWTree_WZ													WZ              4                  	47.13           1000000\n");
    in_file.write("#WWTree_ZZ													ZZ              4                   16.523          985600\n");
    in_file.write("WWTree_WW_excl												WW              4                  	49.997          1924400\n");
    in_file.write("WWTree_WZ_excl												WZ              4                  	10.71           25704656\n");
    in_file.write("WWTree_ZZ_excl												ZZ              4                 	3.22            15301695\n");
    in_file.write("WWTree_sch													STop            7                  	3.65792         998400\n");
    in_file.write("WWTree_tch_bar												STop            7                 	26.0659         1630900\n");
    in_file.write("WWTree_tch													STop            7                  	43.79844        3299200\n");
    in_file.write("WWTree_tWch													STop            7                  	35.6            1000000\n");
    in_file.write("WWTree_tWch_bar												STop            7                  	35.6            999400\n");
    in_file.write("#WWTree_TTbar_amcatnlo										tt_bar          210                	831.76          38475776\n");
    in_file.write("#WWTree_TTbar_madgraph										tt_bar          210                	831.76          10215131\n");
    in_file.write("WWTree_TTbar													tt_bar          210                	831.76          187626200\n");
    in_file.write("#WWTree_Higgs750												ggHx750	      	1			 		0.6398			96200\n");
    in_file.write("#WWTree_VBFHiggs750											qqHx750	      	12		 			0.1915			100000\n");
    in_file.write("#WWTree_Higgs1000											ggHx1000	    1			 		0.1233			400000\n");
    in_file.write("#WWTree_VBFHiggs1000											qqHx1000	    12		 			0.08732			399634\n");
    in_file.write("#WWTree_RSGraviton800										RSGrav800GeV    1                  	1.16691        	31906\n");
    in_file.write("#WWTree_BulkGraviton1400										RSGrav1000      1                 	0.000045798     50000\n");
    in_file.write("#WWTree_WprimeToWZ1400										Wprime1000      1                  	0.031283        50000\n");
    in_file.write("####################################################################################\n");
    in_file.write("######### SIGNAL VALUES\n");
    in_file.write("####################################################################################\n");
    
    for k in range(2):
        
        signalName=sample_value[k][sample_number][2];
        reducedName=sample_value[k][sample_number][3];
        color=sample_value[k][sample_number][4];
        xsec=sample_value[k][sample_number][5];
        NumberEntriesBefore=sample_value[k][sample_number][6];
        in_file.write("%s\t\t\t\t\t\t\t\t%s\t\t\t%.0f\t\t%f\t\t%f\n"%(signalName,reducedName,color,xsec,NumberEntriesBefore));
    
    in_file.close();
    return filename;           
            
            
def make_latex_table(cuts_number,E_data,cuttype,nsample,ofile):
 
    
    sample=E_data[cuttype][0][nsample][0];
    mass_str=str(E_data[cuttype][0][nsample][1]);
    
    ofile.write("\n\n\n\n");
    ofile.write("\changefontsizes{9pt}\n");
    ofile.write("\\begin{frame}[allowframebreaks]\n");
    ofile.write("\changefontsizes{8pt}\n");
    
    if cuttype:
           ofile.write("\\frametitle{%s %s - Control Plots - Single Cut }\n"%(sample,mass_str));   
           ofile.write("\\framesubtitle{$\mu$-channel  Ntuple: \\texttt{WWTree\_22sep\_jecV7\_lowmass/WWTree\_mu/ } }\n");
    else:
           ofile.write("\\frametitle{%s %s - Control Plots - Consecutive Cuts }\n"%(sample,mass_str));   
           ofile.write("\\framesubtitle{$\mu$-channel  Ntuple: \\texttt{WWTree\_22sep\_jecV7\_lowmass/WWTree\_mu/ } }\n");
    
    
    latex_table=[[0 for i in range(number_Events_type+5)] for j in range(cuts_number+1)];
    latex_table[0][0]="Cut";
    latex_table[0][1]="W+Jet Py";
    latex_table[0][2]="W+JetHe";
    latex_table[0][3]="TTbarPo";
    latex_table[0][4]="TTbarMC";
    latex_table[0][5]="VVqcd";
    latex_table[0][6]="WWewk";
    latex_table[0][7]="Stop";
    latex_table[0][8]="AllBkgPy";
    latex_table[0][9]="AllBkgHe";
    latex_table[0][10]="Signalgg";
    latex_table[0][11]="SignalVBF";
    latex_table[0][12]="SigPy";
    latex_table[0][13]="SigHe";
    latex_table[0][14]="$\\frac{S_j}{S_{j-1}}$ Pythia";
    latex_table[0][15]="$\\frac{S_j}{S_{j-1}}$ Herwig";
    
        
    
    c=0;
    for c in range(cuts_number):
        latex_table[c+1][0]=c+1;
        ev=0;
        for ev in range(number_Events_type):
            latex_table[c+1][ev+1]=E_data[cuttype][c][nsample][ev+2];
            #print "cuttype: %.0f\t CutNumber:%.0f\t SampleNumber:%.0f\t Line:%.0f\t VALUE:%f"%(cuttype,c,nsample,ev+2,latex_table[c+1][ev+1])

        sig_n_py=E_data[cuttype][c][nsample][number_Events_type+2]
        sig_n_he=E_data[cuttype][c][nsample][number_Events_type+3]
        latex_table[c+1][12]=sig_n_py;
        latex_table[c+1][13]=sig_n_he;
        
        #latex_table[c+1][12]=latex_table[c+1][11]/(1+TMath.Sqrt(latex_table[c+1][8]));
        #latex_table[c+1][13]=latex_table[c+1][11]/(1+TMath.Sqrt(latex_table[c+1][9]));
    
        
        if c:
           sig_n_1_py=E_data[cuttype][c-1][nsample][number_Events_type+2];
           sig_n_1_he=E_data[cuttype][c-1][nsample][number_Events_type+3];
           
           if (sig_n_1_py and sig_n_1_he):
              sig_rel_py=sig_n_py/sig_n_1_py;
              sig_rel_he=sig_n_he/sig_n_1_he;
           else:
              sig_rel_py=0;
              sig_rel_he=0;
        else:
           sig_rel_py=sig_n_py;
           sig_rel_he=sig_n_he;
           
        latex_table[c+1][14]=sig_rel_py;
        latex_table[c+1][15]=sig_rel_he;
        
       
    if cuttype:
       table_string_name="Single Cut %s%s"%(sample,mass_str);
    else:
       table_string_name="Consecutive Cuts %s%s"%(sample,mass_str);
       

    ofile.write("\n\n\n");
    ofile.write("\\begin{table}[H]\n");
    ofile.write("\\begin{center}\n");
    ofile.write("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
    ofile.write("\hline \multicolumn{7}{|c|}{%s} \\\ \n"%(table_string_name));
    
    '''
    for r in range(cuts_number+1):
        for c in range(number_Events_type+3):
            if not r:
                   print "Riga:%.0f\tColonna:%.0f\tVALUE:%s\n"%(r,c,latex_table[r][c])
            else:
                   print "Riga:%.0f\tColonna:%.0f\tVALUE:%f\n"%(r,c,latex_table[r][c])
    '''
    for r in range(cuts_number+1):
        if r:
           ofile.write("\hline %.0f & %.4f & %.5f & %.5f & %.5f & %.4f & %.5f \\\ \n"%(latex_table[r][0],latex_table[r][1],latex_table[r][3],latex_table[r][4],latex_table[r][5],latex_table[r][6],latex_table[r][7]));
        else:
           ofile.write("\hline %s & %s & %s & %s & %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][1],latex_table[r][3],latex_table[r][4],latex_table[r][5],latex_table[r][6],latex_table[r][7]));
        
            
    
    ofile.write("\hline\n");
    ofile.write("\end{tabular}\n");
    #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
    ofile.write("\end{center}\n");
    ofile.write("\end{table}\n");
    
    ofile.write("\\framebreak\n");
    
    
                
    if cuttype:
       ofile.write("\n\n\n");
       ofile.write("\\begin{table}[H]\n");
       ofile.write("\\begin{center}\n");
       ofile.write("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
       ofile.write("\hline \multicolumn{7}{|c|}{%s} \\\ \n"%(table_string_name));
    
    
    
       for r in range(cuts_number+1):
           if r:
              ofile.write("\hline %.0f & %f & %f & %f & %f & %f & %f \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
           else:
              ofile.write("\hline %s & %s & %s & %s & %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
        
            
    
       ofile.write("\hline\n");
       ofile.write("\end{tabular}\n");
       #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
       ofile.write("\end{center}\n");
       ofile.write("\end{table}\n");
       
    else:
       ofile.write("\n\n\n");
       ofile.write("\\begin{table}[H]\n");
       ofile.write("\\begin{center}\n");
       ofile.write("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
       ofile.write("\hline \multicolumn{7}{|c|}{%s} \\\ \n"%(table_string_name));
    
    
    
       for r in range(cuts_number+1):
           if r:
              ofile.write("\hline %.0f & %f & %f & %f & %f & %f & %f \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
           else:
              ofile.write("\hline %s & %s & %s & %s & %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
        
            
    
       ofile.write("\hline\n");
       ofile.write("\end{tabular}\n");
       #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
       ofile.write("\end{center}\n");
       ofile.write("\end{table}\n");
       
       ofile.write("\\framebreak\n");
       
       ofile.write("\n\n\n");
       ofile.write("\\begin{table}[H]\n");
       ofile.write("\\begin{center}\n");
       ofile.write("\\begin{tabular}{|c|c|c|}\n");
       ofile.write("\hline \multicolumn{3}{|c|}{%s} \\\ \n"%(table_string_name));
    
    
    
       for r in range(cuts_number+1):
           if r:
              ofile.write("\hline %.0f & %f & %f  \\\ \n"%(latex_table[r][0],latex_table[r][14],latex_table[r][15]));
           else:
              ofile.write("\hline %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][14],latex_table[r][15]));
        
            
    
       ofile.write("\hline\n");
       ofile.write("\end{tabular}\n");
       #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
       ofile.write("\end{center}\n");
       ofile.write("\end{table}\n");       
       
    ofile.write("\end{frame}\n");

def replace_latex(in_string):
    out1=in_string.replace("&&", "\&\&");
    out2=out1.replace("_", "\_");
    out3=out2.replace(">", "\\texttt{>}");
    out4=out3.replace("<", "\\texttt{<}");
    return out4


    
    
    
    
   
    
    
    
###########################################################################################
######## MAIN FUNCTION 
###########################################################################################



if __name__ == '__main__':

    lumi=options.lumi;
    lumi_str=str("%.0f"%lumi);
    lumi_str_all=str("%f"%lumi);
    print "\n\n\nWelcome! \nControl Plots Maker!\n"
    print "\nLuminosity:\t%f"%lumi
    
    
    #########################################################
    ######### MAKING DIRECTORY
    #########################################################
    print "\n\n\n----------- Check or making directory ---------------------\n"

    ntuple_dir="output/Ntuple_%s"%(options.ntuple);
    if not os.path.isdir(ntuple_dir):
           pd1 = subprocess.Popen(['mkdir',ntuple_dir]);
           pd1.wait();


    lumi_dir=ntuple_dir+"/Lumi_%s"%(lumi_str);
    if not os.path.isdir(lumi_dir):
           pd2 = subprocess.Popen(['mkdir',lumi_dir]);
           pd2.wait();
       


    #cuts_file_dir="cfg/DataMCComparison_InputCfgFile/MATTEO_Cuts_dir";
    cuts_file_dir="cfg/DataMCComparison_InputCfgFile";
    if not os.path.isdir(cuts_file_dir):
           pd3 = subprocess.Popen(['mkdir',cuts_file_dir]);
           pd3.wait();
       
       
    control_plots_dir=lumi_dir+"/ControlPlots";
    if os.path.isdir(control_plots_dir):
       pd10 = subprocess.Popen(['rm','-r',control_plots_dir]);
       pd10.wait();   
       
    
    #if not os.path.isdir(control_plots_dir):
    pd4 = subprocess.Popen(['mkdir',control_plots_dir]);
    pd4.wait();
    
    
    summaryF = control_plots_dir+"/Summary_ControlPlots.txt";
    output_summaryFile=open(summaryF,'w+');
    output_summaryFile.write("\n\nSUMMARY CONTROL PLOTS\n\n");
    
    summary_latex = control_plots_dir+"/Summary_latex_ControlPlots.txt";
    output_summary_latex_File=open(summary_latex,'w+');
    #output_summaryFile.write("\n\nSUMMARY CONTROL PLOTS\n\n");


    ######################################################
    ######## CUTS ITEMIZE
    ######################################################

    #cuts_itemize=["issignal==1","v_pt>200","pfMET>40","l_pt>40","ungroomed_jet_pt>200","nBTagJet_medium<1","jet_mass pr > 65 && jet mass pr < 105","abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001","jet_tau2tau1 < 0.45"];

    

    conjunction=" && ";
    
    counter=0;
    for line in cuts_itemize:
        counter+=1;
    

    cuts_number=counter;   
    
    
    
    ### Make the SampleListFile
    value = make_sample_table(cuts_number);

    cut_string1="";
    cut_string2="";
    counter=int(0);
    table=value[0];
    n_BulkGraviton=value[1];
    n_Higgs=value[2];
    n_attributes=value[3];
    n_total_sample=n_BulkGraviton+n_Higgs;
    ScaleFactor_table=value[4];
    Events_table_main=value[5];
    
    cuts_table_main=[[0 for i in range(cuts_number)] for j in range(2)];
    for counter in range(cuts_number):
        
        
        ### Make the cutFiles
        ###### index 1: consecutive cuts
        ###### index 2: single cut   
        
        
           
           
        if counter==0:
           cut_string1=cuts_itemize[counter];
           cut_string2=cuts_itemize[counter];
    
        else:
           cut_string1=cut_string1+conjunction+cuts_itemize[counter];
           cut_string2=cuts_itemize[counter];

        
        
        counter_cut_number=counter+1;
        
        cuts_table_main[0][counter]=cut_string1;
        cuts_table_main[1][counter]=cut_string2;
        
        cuts_file1=cuts_file_dir+"/MATTEO_cuts_file1_%s.txt"%(str(counter_cut_number));
        output_cuts_file1=open(cuts_file1,'w+');
        output_cuts_file1.write(cut_string1);
        output_cuts_file1.close();

        cuts_file2=cuts_file_dir+"/MATTEO_cuts_file2_%s.txt"%(str(counter_cut_number));
        output_cuts_file2=open(cuts_file2,'w+');
        output_cuts_file2.write(cut_string2);
        output_cuts_file2.close();
    
        print "\n\n\n\n------------------Making CutsFile---------------------\n"
        print "Total number of cuts:\t%.0f"% cuts_number
        print "Current cut:\t%.0f"%(counter+1)
        print "\nCut file 1:\t%s"%cuts_file1
        print "Cut String 1:\t%s"%cut_string1
        print "\nCut file 2:\t%s"%cuts_file2
        print "Cut String 2:\t%s"%cut_string2
        print "\n------------------------------------------------------"
        
        
        output_summaryFile.write("\n\n\n\n------------------Making CutsFile---------------------\n");
        output_summaryFile.write("Total number of cuts:\t%.0f"% cuts_number);
        output_summaryFile.write("Current cut:\t%.0f"%(counter+1));
        output_summaryFile.write("\nCut file 1:\t%s"%cuts_file1);
        output_summaryFile.write("Cut String 1:\t%s"%cut_string1);
        output_summaryFile.write("\nCut file 2:\t%s"%cuts_file2);
        output_summaryFile.write("Cut String 2:\t%s"%cut_string2);
        output_summaryFile.write("\n------------------------------------------------------");
    

    
    
    
    
        control_cuts1_dir=control_plots_dir+"/Consecutive_Cuts_%s"%(str(counter_cut_number));
        if not os.path.isdir(control_cuts1_dir):
               pd5 = subprocess.Popen(['mkdir',control_cuts1_dir]);
               pd5.wait();
    
        control_cuts2_dir=control_plots_dir+"/Single_Cuts_%s"%(str(counter_cut_number));
        if not os.path.isdir(control_cuts2_dir):
               pd6 = subprocess.Popen(['mkdir',control_cuts2_dir]);
               pd6.wait();
           
        #cuts_file1=cuts_file_dir+"/MATTEO_cuts_file1_%s.txt"%(str(counter+1));
        #cuts_file2=cuts_file_dir+"/MATTEO_cuts_file2_%s.txt"%(str(counter+1));
    
        print "\n\n\n----------------------Making Control Plots-----------------------\n"
        print "\nProcessing Cuts %0.f of %0.f\n"%(counter_cut_number,cuts_number)
        print "Using Cuts File 1:\t%s"%cuts_file1
        print "Using Cuts File 2:\t%s"%cuts_file2
        print "\n-----------------------------------------------------------------\n\n\n"
        
        
        output_summaryFile.write("\n\n\n----------------------Making Control Plots-----------------------\n");
        output_summaryFile.write("\nProcessing Cuts %0.f of %0.f\n"%(counter_cut_number,cuts_number));
        output_summaryFile.write("Using Cuts File 1:\t%s"%cuts_file1);
        output_summaryFile.write("Using Cuts File 2:\t%s"%cuts_file2);
        output_summaryFile.write("\n-----------------------------------------------------------------\n\n\n");
        
        
        
        
        
        
        
        
        for k in range(n_total_sample):
            SampleList_filename=make_SampleList_file(options.ntuple,k,table);
            print "\n\nMade SampleList File:\t%s\n"%SampleList_filename
            output_summaryFile.write("\n\nMade SampleList File:\t%s\n"%SampleList_filename);
            
            sample=table[0][k][0];
            mass=table[0][k][1];        
            mass_str=str("%.0f"%mass)
            output_summaryFile.write("\n-------------------------------------------------------------\n");
            output_summaryFile.write("PROCESSING:\t%s%s\n"%(sample,mass_str));
            #output_summaryFile.write("\n-----------------------------------------------------------------\n\n\n");
            
            final_dir1=control_cuts1_dir+"/"+sample+mass_str;
            if not os.path.isdir(final_dir1):
               pd7 = subprocess.Popen(['mkdir',final_dir1]);
               pd7.wait();
                 
            final_dir2=control_cuts2_dir+"/"+sample+mass_str;
            if not os.path.isdir(final_dir2):
               pd8 = subprocess.Popen(['mkdir',final_dir2]);
               pd8.wait();
                            
            output_summaryFile.write("\n-----------------------------------------------------------------\n");
            output_summaryFile.write("CONSECUTIVE CUTS");
            output_summaryFile.write("\n-----------------------------------------------------------------\n");    
            output_summary_latex_File.write("\nCONSECUTIVE CUT\n");
            
            #make_Input_file(ntuple,channel,lumi,cuts_file,SampleList_filename,table,k,ScaleFactor_table,cut_number,cuttype):
            info_input_file1= make_Input_file(options.ntuple,options.channel,options.lumi,cuts_file1,SampleList_filename,table,k,ScaleFactor_table,counter_cut_number,1);
            #data_mc_comparison_input_file1=info_input_file1[0];
            Events_output=run_log(info_input_file1,final_dir1,table,k,options.lumi,output_summaryFile,Events_table_main,1,counter,output_summary_latex_File,cuts_number)
            
            output_summaryFile.write("\n-----------------------------------------------------------------\n");
            output_summaryFile.write("SINGLE CUT");
            output_summaryFile.write("\n-----------------------------------------------------------------\n");  
                   
            output_summary_latex_File.write("\nSINGLE CUT\n");
                
            info_input_file2 = make_Input_file(options.ntuple,options.channel,options.lumi,cuts_file2,SampleList_filename,table,k,ScaleFactor_table,counter_cut_number,2);
            #data_mc_comparison_input_file2=info_input_file2[0];
            Events_output=run_log(info_input_file2,final_dir2,table,k,options.lumi,output_summaryFile,Events_table_main,2,counter,output_summary_latex_File,cuts_number)
            
            
    
    
    
    
    latex_file = control_plots_dir+"/Latex_ControlPlots.txt";
    latex_out_file=open(latex_file,'w+');
    
    latex_out_file.write("\documentclass{beamer}\n");
    latex_out_file.write("\usetheme{Boadilla}\n");
    latex_out_file.write("\usecolortheme{seahorse}\n");
    latex_out_file.write("\\title{ControlPlots}\n");
    latex_out_file.write("\\author{Matteo Rappo}\n");
    latex_out_file.write("\setbeamertemplate{navigation symbols}{}\n");
    latex_out_file.write("\usepackage[latin1]{inputenc}\n");
    latex_out_file.write("\usepackage[english,italian]{babel}\n");
    latex_out_file.write("\usepackage{amsmath}\n");
    latex_out_file.write("\usepackage{enumerate}\n");
    latex_out_file.write("\usepackage{amsfonts}\n");
    latex_out_file.write("\usepackage{amssymb}\n");
    latex_out_file.write("\usepackage{float}\n");
    latex_out_file.write("\usepackage{placeins}\n");
    latex_out_file.write("\usepackage{subfig}\n");
    latex_out_file.write("\usepackage{multirow,makecell}\n");
    latex_out_file.write("\usepackage{array,booktabs}\n");
    latex_out_file.write("\usepackage{comment}\n");
    latex_out_file.write("\usepackage{scrextend}\n");
    latex_out_file.write("\usepackage{verbatim,longtable}\n");
    latex_out_file.write("\setbeamertemplate{caption}[numbered]\n");
    latex_out_file.write("\\newcolumntype{P}[1]{>{\centering\\arraybackslash}p{#1}}\n");
    latex_out_file.write("\\newcolumntype{M}[1]{>{\centering\\arraybackslash}m{#1}}\n");
    latex_out_file.write("\\newcolumntype{D}[1]{>{\\arraybackslash}m{#1}}\n");
    latex_out_file.write("\\newcolumntype{C}[1]{>{\centering\let\\newline\\\\arraybackslash\hspace{0pt}}m{#1}}\n");
    latex_out_file.write("\n");
    latex_out_file.write("\n");
    latex_out_file.write("\n");
    latex_out_file.write("\changefontsizes{9pt}\n");
    latex_out_file.write("\\begin{document}\n");
    latex_out_file.write("\n");
    latex_out_file.write("\n");
     






    latex_out_file.write("\n\n\n");
    latex_out_file.write("\changefontsizes{9pt}\n");
    latex_out_file.write("\\begin{frame}[t,allowframebreaks]\n");
    latex_out_file.write("\\frametitle{Control Plots - Consecutive Cuts }\n");   
    latex_out_file.write("\\framesubtitle{$\mu$-channel  Ntuple: \\texttt{WWTree\_22sep\_jecV7\_lowmass/WWTree\_mu/ } }\n");
    latex_out_file.write("\changefontsizes{7pt}\n");
    #latex_out_file.write("\\begin{table}[H]\n");
    #latex_out_file.write("\\begin{center}\n");
    latex_out_file.write("\\begin{longtable}{|M{10pt}|D{310pt}|}\n");
    latex_out_file.write("\hline \multicolumn{2}{|c|}{Consecutive Cuts} \\\ \n");
    
    
    j=0;
    for j in range(cuts_number):
        cn=j+1;
        tmp=replace_latex(cuts_table_main[0][j]);
        #tmp2="\\texttt{"+tmp1+"}";
        latex_out_file.write("\hline %.0f & %s \\\ \n"%(cn,tmp));
         
            
    latex_out_file.write("\hline\n");
    #latex_out_file.write("\end{tabular}\n");
    #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
    #latex_out_file.write("\end{center}\n");
    latex_out_file.write("\end{longtable}\n");
    latex_out_file.write("\end{frame}\n");
    
    
    latex_out_file.write("\n\n\n");
    latex_out_file.write("\changefontsizes{9pt}\n");
    latex_out_file.write("\\begin{frame}\n");
    latex_out_file.write("\\frametitle{Control Plots - Single Cut }\n");   
    latex_out_file.write("\\framesubtitle{$\mu$-channel  Ntuple: \\texttt{WWTree\_22sep\_jecV7\_lowmass/WWTree\_mu/ } }\n");
    latex_out_file.write("\changefontsizes{7pt}\n");
    latex_out_file.write("\\begin{table}[H]\n");
    latex_out_file.write("\\begin{center}\n");
    latex_out_file.write("\\begin{tabular}{|M{10pt}|D{310pt}|}\n");
    latex_out_file.write("\hline \multicolumn{2}{|c|}{Single Cut} \\\ \n");
    
    j=0;
    for j in range(cuts_number):
        cn=j+1;
        tmp=replace_latex(cuts_table_main[1][j]);
        #tmp2="\\texttt{"+tmp1+"}";
        latex_out_file.write("\hline %.0f & %s \\\ \n"%(cn,tmp));
          
            
    
    latex_out_file.write("\hline\n");
    latex_out_file.write("\end{tabular}\n");
    #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
    latex_out_file.write("\end{center}\n");
    latex_out_file.write("\end{table}\n");
    latex_out_file.write("\end{frame}\n");




            


    
    
    
    nsample=0;
    for nsample in range(n_total_sample):
        
        make_latex_table(cuts_number,Events_output,0,nsample,latex_out_file);
        
        make_latex_table(cuts_number,Events_output,1,nsample,latex_out_file);
        
    latex_out_file.write("\end{document}\n");
    output_summary_latex_File.close()
    latex_out_file.close();  
    output_summaryFile.close();
        
        
        
        
        
        
         
        
        
      
            
            
            
            
            
    
