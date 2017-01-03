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

Events_type_global=["Wjets_Pythia_Events_g",      # 0
                    "Wjets_Herwig_Events_g",      # 1
                    "TTbar_Powegh_Events_g",      # 2
                    "TTbar_MC_Events_g",          # 3
                    "VV_QCD_Events_g",            # 4
                    "WW_EWK_Events_g",            # 5
                    "STop_Events_g",              # 6
                    "All_bkg_Pythia_g",           # 7
                    "All_bkg_Herwig_g",           # 8
                    "Signal_gg_g",                # 9
                    "Signal_VBF_g"];              # 10

number_Events_type=11;



       
BulkGraviton_xsec=[0.177400,0.0331548,0.008993];
VBF_BulkGraviton_xsec=[0.01089,0.00217,0.000655];
Higgs_xsec=[0.33639,0.06765];
VBF_Higgs_xsec=[0.03354,0.02375];

xsec_BulkGraviton=[0.177400,0.0331548,0.008993];
xsec_VBF_BulkGraviton=[0.01089,0.00217,0.000655];
xsec_Higgs=[0.33639,0.06765];
xsec_VBF_Higgs=[0.03354,0.02375];
		    
NumEntriesBefore_BulkGraviton=[49600,50000,50000];
NumEntriesBefore_VBF_BulkGraviton=[50000,50000,50000];
NumEntriesBefore_Higgs=[399600,400000];
NumEntriesBefore_VBF_Higgs=[398400,400000];
		    
ScaleFactor_BulkGraviton=[900,2000,6000];
ScaleFactor_Higgs=[25,120]

sampleValue =[["BulkGraviton","BG_600",600,xsec_BulkGraviton[0],NumEntriesBefore_BulkGraviton[0],ScaleFactor_BulkGraviton[0]],
              ["BulkGraviton","BG_800",800,xsec_BulkGraviton[1],NumEntriesBefore_BulkGraviton[1],ScaleFactor_BulkGraviton[1]],
              ["BulkGraviton","BG_1000",1000,xsec_BulkGraviton[2],NumEntriesBefore_BulkGraviton[2],ScaleFactor_BulkGraviton[2]],
              ["Higgs","H_650",650,xsec_Higgs[0],NumEntriesBefore_Higgs[0],ScaleFactor_Higgs[0]],
              ["Higgs","H_1000",1000,xsec_Higgs[1],NumEntriesBefore_Higgs[1],ScaleFactor_Higgs[1]]];
              
              
sampleValue_VBF =[["BulkGraviton","VBF_BG_600",600,xsec_VBF_BulkGraviton[0],NumEntriesBefore_VBF_BulkGraviton[0],ScaleFactor_BulkGraviton[0]],
                  ["BulkGraviton","VBF_BG_800",800,xsec_VBF_BulkGraviton[1],NumEntriesBefore_VBF_BulkGraviton[1],ScaleFactor_BulkGraviton[1]],
                  ["BulkGraviton","VBF_BG_1000",1000,xsec_VBF_BulkGraviton[2],NumEntriesBefore_VBF_BulkGraviton[2],ScaleFactor_BulkGraviton[2]],
                  ["Higgs","VBF_H_650",650,xsec_VBF_Higgs[0],NumEntriesBefore_VBF_Higgs[0],ScaleFactor_Higgs[0]],
                  ["Higgs","VBF_H_1000",1000,xsec_VBF_Higgs[1],NumEntriesBefore_VBF_Higgs[1],ScaleFactor_Higgs[1]]];


total_sample_value=[sampleValue,sampleValue_VBF];

lumi=options.lumi;
lumi_str=str("%.0f"%lumi);
cuts_itemize=["issignal==1 && abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001 && v_pt>200 && pfMET>40 && l_pt>40 && ungroomed_jet_pt>200 && nBTagJet_medium <1 && jet_mass_pr > 65 && jet_mass_pr < 105","jet_tau2tau1 < 0.45"];
#cuts_itemize=["deltaR_lak8jet>(TMath::Pi()/2.0) && abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001","TMath::Abs(deltaphi_METak8jet)>2.0","TMath::Abs(deltaphi_Vak8jet)>2.0","v_pt>200","pfMET>40","l_pt>40","ungroomed_jet_pt>200","nBTagJet_medium <1","jet_tau2tau1 < 0.45","jet_mass_pr > 65 && jet_mass_pr < 105 "];
#cuts_itemize=["issignal==1 && v_pt>200","pfMET>40 && l_pt>40 && ungroomed_jet_pt>200 && nBTagJet_medium <1 && jet_mass_pr > 65 && jet_mass_pr < 105 && abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)>0.001","jet_tau2tau1 < 0.45"];


###########################################################################################
######## FUNCTION DEFINITION
###########################################################################################



    
    
    




def run_log(Sample_Number_ll,Cut_Number_ll,C_Type_ll,Cut_Total_Number_ll,Significance_Table_ll,Input_Cfg_File_ll,Latex_File_ll,Summary_Output_File_ll,Directory_Path_ll):
    
    Cut_Type_ll=C_Type_ll-1;
    Sample_ll=sampleValue[Sample_Number_ll][0];
    Mass_ll=sampleValue[Sample_Number_ll][2];        
    Mass_ll_str=str("%.0f"%Mass_ll);
    process_name=Sample_ll+Mass_ll_str;
   
  
    
    Latex_File_ll.write("%s\n"%Sample_ll);
    Latex_File_ll.write("%s\n"%Mass_ll_str);
    Latex_File_ll.write("%s\n"%Cut_Type_ll);
    Latex_File_ll.write("%.0f\n"%Cut_Number_ll);
            
            
    Log_Dir_ll=Directory_Path_ll+"/log";
    if not os.path.isdir(Log_Dir_ll):
           pd1 = subprocess.Popen(['mkdir',Log_Dir_ll]);
           pd1.wait();
            
    Data_Dir_ll=Directory_Path_ll+"/data";
    if not os.path.isdir(Data_Dir_ll):
           pd2 = subprocess.Popen(['mkdir',Data_Dir_ll]);
           pd2.wait();
    
    cn_print=Cut_Number_ll+1; 
    if Cut_Type_ll:
       print_string=("\tSingle CUT: %.0f of %.0f")%(cn_print,Cut_Total_Number_ll);     
    else:
       print_string=("\tRecursive CUTS: %.0f of %.0f")%(cn_print,Cut_Total_Number_ll);
    
    log_file1=Log_Dir_ll+"/Log_ControlPlots_%s_%s%s.txt"%(options.channel,Sample_ll,Mass_ll_str)
    log_file2=Log_Dir_ll+"/Error_Log_ControlPlots_%s_%s%s.txt"%(options.channel,Sample_ll,Mass_ll_str)
     
          
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
          
    sys.stdout.write("\n\n\n\n-------------------------------------------------------------------------------------------------------\n\n")
    sys.stdout.write("STARTING\t\t")
    sys.stdout.write(process_name)
    sys.stdout.write(print_string)
    sys.stdout.write("\n\n")
    
    #Summary_Output_File_ll.write("\n\n--------------------------------\n\n")
    Summary_Output_File_ll.write("STARTING\t\t")
    Summary_Output_File_ll.write(process_name)
    Summary_Output_File_ll.write(print_string)
    Summary_Output_File_ll.write("\n\n--------------------------------\n\n")                
    
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
    
    
    
    p1 = subprocess.Popen(['./bin/DataMCComparisonPlot.exe',Input_Cfg_File_ll[0]],stdout=subprocess.PIPE,stderr=output_log2)
    p1.wait()
    start=0;
    for line in p1.stdout:
        #sys.stdout.write(line)
        output_log1.write(line)
        if line.find('Event Scaled To Lumi') !=-1:
           start=1;
        if start:   
           Summary_Output_File_ll.write(line);
        counter_events_type=0;
        for ev in Events_type:
            if line.find(ev) !=-1:
               cut_string = line.split(ev);
               new_string = cut_string[1]
               print line
               val=float(new_string);
               Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][counter_events_type]=val;
               #print "Cut_Type_ll: %.0f\t Cut_Number_ll:%.0f\t Sample_Number_ll:%.0f\t Line:%.0f\t VALUE:%f"%(Cut_Type_ll,cnumber,k,cn+2,val)
               Latex_File_ll.write("%f\n"%val);
               
            counter_events_type=counter_events_type+1;
    p1.wait();    

    sig_Pytia=Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-1]/(1+TMath.Sqrt(Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-4]));
    sig_Herwig=Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-1]/(1+TMath.Sqrt(Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-3]));
    Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type]=sig_Pytia;
    Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type+1]=sig_Herwig;
    print "\n\n------------Calculate Significance----------------\n"
    print "\nqq signal: %f\n"% Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-1]
    print "\nPythia Bkg: %f\n"% Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-4]
    print "\nHewig Bkg: %f\n"% Significance_Table_ll[Cut_Type_ll][Cut_Number_ll][Sample_Number_ll][number_Events_type-3]
    print "\nSig Pythia: %f\n"% sig_Pytia
    print "\nSig Herwig: %f\n\n"% sig_Herwig    
    #print "\n\n--------------------------------------------------\n\n"
    Latex_File_ll.write("%f\n"%sig_Pytia);
    Latex_File_ll.write("%f\n"%sig_Herwig);
    Summary_Output_File_ll.write("Significance Pythia:    %f\n"%sig_Pytia);
    Summary_Output_File_ll.write("Significance Herwig:    %f\n"%sig_Herwig);
    
    
    Latex_File_ll.write("\n\n");
    
    
    print "\nInputFile: %s"%Input_Cfg_File_ll[0]   
    print "\nOutputDirectory: %s"%Input_Cfg_File_ll[1]
    
    
    print "\nVBF_Sample: %s \tReducedName: %s \tMass: %.0f \txSec: %f \tNumbEntBefore: %.0f \tScaleFactor: %.0f\n"%(total_sample_value[1][Sample_Number_ll][0],total_sample_value[1][Sample_Number_ll][1],total_sample_value[1][Sample_Number_ll][2],total_sample_value[1][Sample_Number_ll][3],total_sample_value[1][Sample_Number_ll][4],total_sample_value[1][Sample_Number_ll][5])
    
    print "\nNormal_Sample: %s \tReducedName: %s \tMass: %.0f \txSec: %f \tNumbEntBefore: %.0f \tScaleFactor: %.0f\n"%(total_sample_value[0][Sample_Number_ll][0],total_sample_value[0][Sample_Number_ll][1],total_sample_value[0][Sample_Number_ll][2],total_sample_value[0][Sample_Number_ll][3],total_sample_value[0][Sample_Number_ll][4],total_sample_value[0][Sample_Number_ll][5])
    
    path_dir_in_tmp=Input_Cfg_File_ll[1];
    path_dir_in=path_dir_in_tmp+"Run2_MCDataComparisonRSGraviton2000_mu_plot"; 
    #path_dir_in="output/run2/MCDATAComparisonPlot_mu_22sep_%s/Run2_MCDataComparisonRSGraviton2000_mu_plot"%process_name;
    path_dir_out=Directory_Path_ll;
    p2 = subprocess.Popen(['cp','-r',path_dir_in,path_dir_out])
    p2.wait()



    data_in=Directory_Path_ll+"/Run2_MCDataComparisonRSGraviton2000_mu_plot";
    data_out=Directory_Path_ll+"/data";
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
            
                        
    #sys.stdout.write("\n\n------------------------------------------------------------------\n\n")
    sys.stdout.write("ENDED\t\t")
    sys.stdout.write(process_name)
    sys.stdout.write(print_string)
    sys.stdout.write("\n\n-------------------------------------------------------------------------------------------------------\n\n\n\n")
    
    Summary_Output_File_ll.write("\n\n--------------------------------\n\n")
    Summary_Output_File_ll.write("ENDED\t\t")
    Summary_Output_File_ll.write(process_name)
    Summary_Output_File_ll.write(print_string)
    Summary_Output_File_ll.write("\n\n--------------------------------\n\n")
    
    
    return Significance_Table_ll;
    
     
      

        
      
    
            


         
def make_latex_table(Cuts_Total_Number_tk,Sig_Data_Table_tk,ctype_tk,Sample_Number_tk,Output_File_tk,NtupleTTName_tk):
 
    Cut_Type_tk=ctype_tk-1;
    Sample_tk=sampleValue[Sample_Number_tk][0];
    Mass_Str_tk=str(sampleValue[Sample_Number_tk][2]);
    
    Output_File_tk.write("\n\n\n\n");
    Output_File_tk.write("\changefontsizes{9pt}\n");
    Output_File_tk.write("\\begin{frame}[allowframebreaks]\n");
    Output_File_tk.write("\changefontsizes{8pt}\n");
    
    if Cut_Type_tk:
           Output_File_tk.write("\\frametitle{%s %s - Control Plots - Single Cut }\n"%(Sample_tk,Mass_Str_tk));   
           Output_File_tk.write("\\framesubtitle{$\mu$-channel  Ntuple: "+NtupleTTName_tk+" }\n");
    else:
           Output_File_tk.write("\\frametitle{%s %s - Control Plots - Consecutive Cuts }\n"%(Sample_tk,Mass_Str_tk));   
           Output_File_tk.write("\\framesubtitle{$\mu$-channel  Ntuple: "+NtupleTTName_tk+" }\n");
    
    
    latex_table=[[0 for i in range(number_Events_type+5)] for j in range(Cuts_Total_Number_tk+1)];
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
    
        
    
    n_cut=0;
    for n_cut in range(Cuts_Total_Number_tk):
        latex_table[n_cut+1][0]=n_cut+1;
        ev=0;
        for ev in range(number_Events_type):
            latex_table[n_cut+1][ev+1]=Sig_Data_Table_tk[Cut_Type_tk][n_cut][Sample_Number_tk][ev];
            #print "Cut_Type_tk: %.0f\t CutNumber:%.0f\t SampleNumber:%.0f\t Line:%.0f\t VALUE:%f"%(Cut_Type_tk,c,Sample_Number_tk,ev+2,latex_table[c+1][ev+1])

        sig_n_py=Sig_Data_Table_tk[Cut_Type_tk][n_cut][Sample_Number_tk][number_Events_type]
        sig_n_he=Sig_Data_Table_tk[Cut_Type_tk][n_cut][Sample_Number_tk][number_Events_type+1]
        latex_table[n_cut+1][12]=sig_n_py;
        latex_table[n_cut+1][13]=sig_n_he;
        
        #latex_table[c+1][12]=latex_table[c+1][11]/(1+TMath.Sqrt(latex_table[c+1][8]));
        #latex_table[c+1][13]=latex_table[c+1][11]/(1+TMath.Sqrt(latex_table[c+1][9]));
    
        
        if n_cut:
           sig_n_1_py=Sig_Data_Table_tk[Cut_Type_tk][n_cut-1][Sample_Number_tk][number_Events_type];
           sig_n_1_he=Sig_Data_Table_tk[Cut_Type_tk][n_cut-1][Sample_Number_tk][number_Events_type+1];
           
           if (sig_n_1_py and sig_n_1_he):
              sig_rel_py=sig_n_py/sig_n_1_py;
              sig_rel_he=sig_n_he/sig_n_1_he;
           else:
              sig_rel_py=0;
              sig_rel_he=0;
        else:
           sig_rel_py=sig_n_py;
           sig_rel_he=sig_n_he;
           
        latex_table[n_cut+1][14]=sig_rel_py;
        latex_table[n_cut+1][15]=sig_rel_he;
        
       
    if Cut_Type_tk:
       table_string_name="Single Cut %s%s"%(Sample_tk,Mass_Str_tk);
    else:
       table_string_name="Consecutive Cuts %s%s"%(Sample_tk,Mass_Str_tk);
       

    Output_File_tk.write("\n\n\n");
    Output_File_tk.write("\\begin{table}[H]\n");
    Output_File_tk.write("\\begin{center}\n");
    Output_File_tk.write("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
    Output_File_tk.write("\hline \multicolumn{7}{|c|}{%s} \\\ \n"%(table_string_name));
    
    '''
    for r in range(Cuts_Total_Number_tk+1):
        for c in range(number_Events_type+3):
            if not r:
                   print "Riga:%.0f\tColonna:%.0f\tVALUE:%s\n"%(r,c,latex_table[r][c])
            else:
                   print "Riga:%.0f\tColonna:%.0f\tVALUE:%f\n"%(r,c,latex_table[r][c])
    '''
    for r in range(Cuts_Total_Number_tk+1):
        if r:
           Output_File_tk.write("\hline %.0f & %.4f & %.5f & %.5f & %.5f & %.4f & %.5f \\\ \n"%(latex_table[r][0],latex_table[r][1],latex_table[r][3],latex_table[r][4],latex_table[r][5],latex_table[r][6],latex_table[r][7]));
        else:
           Output_File_tk.write("\hline %s & %s & %s & %s & %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][1],latex_table[r][3],latex_table[r][4],latex_table[r][5],latex_table[r][6],latex_table[r][7]));
        
            
    
    Output_File_tk.write("\hline\n");
    Output_File_tk.write("\end{tabular}\n");
    #Output_File_tk.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
    Output_File_tk.write("\end{center}\n");
    Output_File_tk.write("\end{table}\n");
    
    Output_File_tk.write("\\framebreak\n");
    
    
                
    if Cut_Type_tk:
       Output_File_tk.write("\n\n\n");
       Output_File_tk.write("\\begin{table}[H]\n");
       Output_File_tk.write("\\begin{center}\n");
       Output_File_tk.write("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
       Output_File_tk.write("\hline \multicolumn{7}{|c|}{%s} \\\ \n"%(table_string_name));
    
    
       r=0;  
       for r in range(Cuts_Total_Number_tk+1):
           if r:
              Output_File_tk.write("\hline %.0f & %f & %f & %f & %f & %f & %f \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
           else:
              Output_File_tk.write("\hline %s & %s & %s & %s & %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
        
            
    
       Output_File_tk.write("\hline\n");
       Output_File_tk.write("\end{tabular}\n");
       #Output_File_tk.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
       Output_File_tk.write("\end{center}\n");
       Output_File_tk.write("\end{table}\n");
       
    else:
       Output_File_tk.write("\n\n\n");
       Output_File_tk.write("\\begin{table}[H]\n");
       Output_File_tk.write("\\begin{center}\n");
       Output_File_tk.write("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
       Output_File_tk.write("\hline \multicolumn{7}{|c|}{%s} \\\ \n"%(table_string_name));
    
    
       r=0;
       for r in range(Cuts_Total_Number_tk+1):
           if r:
              Output_File_tk.write("\hline %.0f & %f & %f & %f & %f & %f & %f \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
           else:
              Output_File_tk.write("\hline %s & %s & %s & %s & %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][8],latex_table[r][9],latex_table[r][10],latex_table[r][11],latex_table[r][12],latex_table[r][13]));
        
            
    
       Output_File_tk.write("\hline\n");
       Output_File_tk.write("\end{tabular}\n");
       #Output_File_tk.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
       Output_File_tk.write("\end{center}\n");
       Output_File_tk.write("\end{table}\n");
       
       Output_File_tk.write("\\framebreak\n");
       
       Output_File_tk.write("\n\n\n");
       Output_File_tk.write("\\begin{table}[H]\n");
       Output_File_tk.write("\\begin{center}\n");
       Output_File_tk.write("\\begin{tabular}{|c|c|c|}\n");
       Output_File_tk.write("\hline \multicolumn{3}{|c|}{%s} \\\ \n"%(table_string_name));
    
    
       r=0;
       for r in range(Cuts_Total_Number_tk+1):
           if r:
              Output_File_tk.write("\hline %.0f & %f & %f  \\\ \n"%(latex_table[r][0],latex_table[r][14],latex_table[r][15]));
           else:
              Output_File_tk.write("\hline %s & %s & %s \\\ \n"%(latex_table[r][0],latex_table[r][14],latex_table[r][15]));
        
            
    
       Output_File_tk.write("\hline\n");
       Output_File_tk.write("\end{tabular}\n");
       #Output_File_tk.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
       Output_File_tk.write("\end{center}\n");
       Output_File_tk.write("\end{table}\n");       
       
    Output_File_tk.write("\end{frame}\n");


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

    Lumi_mm=options.lumi;
    Lumi_mm_str=str("%.0f"%Lumi_mm);
    Lumi_mm_str_all=str("%f"%Lumi_mm);
    print "\n\n\nWelcome! \nControl Plots Maker!\n"
    print "\nNtuple:\t%s\n"%options.ntuple
    print "Luminosity:\t%f"%Lumi_mm

    
    
    
    
    #########################################################
    ######### CHECK INITIAL VALUE
    #########################################################    
    
    
    for i in range(2):
        for k in range(5):
        
            if i:
               print "\nVBF_Sample: %s \tReducedName: %s \tMass: %.0f \txSec: %f \tNumbEntBefore: %.0f \tScaleFactor: %.0f\n"%(total_sample_value[i][k][0],total_sample_value[i][k][1],total_sample_value[i][k][2],total_sample_value[i][k][3],total_sample_value[i][k][4],total_sample_value[i][k][5])
            else:
               print "\nNormal_Sample: %s \tReducedName: %s \tMass: %.0f \txSec: %f \tNumbEntBefore: %.0f \tScaleFactor: %.0f\n"%(total_sample_value[i][k][0],total_sample_value[i][k][1],total_sample_value[i][k][2],total_sample_value[i][k][3],total_sample_value[i][k][4],total_sample_value[i][k][5])
    
    
    
    
    
    
    
    
    
    #########################################################
    ######### MAKING DIRECTORY
    #########################################################
    print "\n\n\n----------- Check or making directory ---------------------\n"

    Ntuple_Dir_mm="output/Ntuple_%s"%(options.ntuple);
    if not os.path.isdir(Ntuple_Dir_mm):
           pd1 = subprocess.Popen(['mkdir',Ntuple_Dir_mm]);
           pd1.wait();


    Lumi_Dir_mm=Ntuple_Dir_mm+"/Lumi_%s"%(Lumi_mm_str);
    if not os.path.isdir(Lumi_Dir_mm):
           pd2 = subprocess.Popen(['mkdir',Lumi_Dir_mm]);
           pd2.wait();
       


    Cuts_File_Dir_mm="cfg/DataMCComparison_InputCfgFile";
    if not os.path.isdir(Cuts_File_Dir_mm):
           pd3 = subprocess.Popen(['mkdir',Cuts_File_Dir_mm]);
           pd3.wait();
       
       
    Control_Plots_Dir_mm=Lumi_Dir_mm+"/ControlPlots";
    if os.path.isdir(Control_Plots_Dir_mm):
       pd10 = subprocess.Popen(['rm','-r',Control_Plots_Dir_mm]);
       pd10.wait();   
       
    
    pd4 = subprocess.Popen(['mkdir',Control_Plots_Dir_mm]);
    pd4.wait();
    
    cfg_file_removal=Cuts_File_Dir_mm+"/MATTEO_*";
    pd5b = subprocess.Popen(['rm','-r',cfg_file_removal]);
    pd5b.wait();
    
    
    #########################################################
    ######### MAKE OUTPUT FILE
    #########################################################    
    summaryF_mm = Control_Plots_Dir_mm+"/Summary_ControlPlots.txt";
    Output_Summary_File_mm=open(summaryF_mm,'w+');
    Output_Summary_File_mm.write("\n\nSUMMARY CONTROL PLOTS\n\n");
    
    summary_latex_mm = Control_Plots_Dir_mm+"/Summary_latex_ControlPlots.tex";
    Output_Summary_Latex_File_mm=open(summary_latex_mm,'w+');
    
    latex_file = Control_Plots_Dir_mm+"/Latex_ControlPlots.tex";
    Output_Beamer_Latex_File_mm=open(latex_file,'w+');




    

    
    #########################################################
    ######### MAKE CONFIGURATION FILE
    ######################################################### 
    
    # Make Cuts File
    counter=0;
    for line in cuts_itemize:
        counter=counter+1;
    

    Cuts_Total_Number_mm=counter;  
    
    cut_string1="";
    cut_string2="";
    conjunction=" && ";
    
    i=j=0;
    Cuts_filename_table=[[0 for j in range(Cuts_Total_Number_mm)]for i in range(2)]; 
    
    i=j=0;
    cuts_table_main=[[0 for j in range(Cuts_Total_Number_mm)]for i in range(2)];
    cut_counter=0;
    for cut_counter in range(Cuts_Total_Number_mm):
        
        
        ### Make the cutFiles
        ###### index 1: consecutive cuts
        ###### index 2: single cut   
        
        
           
           
        if cut_counter==0:
           cut_string1=cuts_itemize[cut_counter];
           cut_string2=cuts_itemize[cut_counter];
    
        else:
           cut_string1=cut_string1+conjunction+cuts_itemize[cut_counter];
           cut_string2=cuts_itemize[cut_counter];

        
        
        Plus_Cut_Counter=cut_counter+1;
        
        cuts_table_main[0][cut_counter]=cut_string1;
        cuts_table_main[1][cut_counter]=cut_string2;
        
        cuts_file1=Cuts_File_Dir_mm+"/MATTEO_cuts_file1_%s.txt"%(str(Plus_Cut_Counter));
        output_cuts_file1=open(cuts_file1,'w+');
        output_cuts_file1.write(cut_string1);
        output_cuts_file1.close();

        cuts_file2=Cuts_File_Dir_mm+"/MATTEO_cuts_file2_%s.txt"%(str(Plus_Cut_Counter));
        output_cuts_file2=open(cuts_file2,'w+');
        output_cuts_file2.write(cut_string2);
        output_cuts_file2.close();
    
        print "\n\n\n\n------------------Making CutsFile---------------------\n"
        print "Total number of cuts:\t%.0f"% Cuts_Total_Number_mm
        print "Current cut:\t%.0f"%(Plus_Cut_Counter)
        print "\nCut file 1:\t%s"%cuts_file1
        print "Cut String 1:\t%s"%cut_string1
        print "\nCut file 2:\t%s"%cuts_file2
        print "Cut String 2:\t%s"%cut_string2
        print "\n------------------------------------------------------"
        
        
        Output_Summary_File_mm.write("\n\n\n\n------------------Making CutsFile---------------------\n");
        Output_Summary_File_mm.write("Total number of cuts:\t%.0f"% Cuts_Total_Number_mm);
        Output_Summary_File_mm.write("Current cut:\t%.0f"%(Plus_Cut_Counter));
        Output_Summary_File_mm.write("\nCut file 1:\t%s"%cuts_file1);
        Output_Summary_File_mm.write("Cut String 1:\t%s"%cut_string1);
        Output_Summary_File_mm.write("\nCut file 2:\t%s"%cuts_file2);
        Output_Summary_File_mm.write("Cut String 2:\t%s"%cut_string2);
        Output_Summary_File_mm.write("\n------------------------------------------------------");
        
        
        Cuts_filename_table[0][cut_counter]=cuts_file1;
        Cuts_filename_table[1][cut_counter]=cuts_file2;
    
    
    
    # Make VariableList
    FileName_VariableList_mm="cfg/DataMCComparison_InputCfgFile/MATTEO_VariableList_76x.txt";
    Output_VariableList_mm=open(FileName_VariableList_mm,'w+');
    Output_VariableList_mm.write("############################################################################\n");
    Output_VariableList_mm.write("##  Variable						Nbin		Min		Max			Label\n");
    Output_VariableList_mm.write("############################################################################\n");
    Output_VariableList_mm.write("# l_pt							50			0		1000		pT_{l}_(GeV)\n");
    Output_VariableList_mm.write("# nPV								25			0		50			nPV\n");
    Output_VariableList_mm.write("# l_pt							25			0		500			pT_{l}_(GeV)\n");
    Output_VariableList_mm.write("# l_eta							25			-2.5	2.5			#eta_{l}\n");
    Output_VariableList_mm.write("# l_eta							20			-2.5	2.5			#eta_{l}\n");
    Output_VariableList_mm.write("# l_phi							30			-3.14	3.14		#phi_{l}\n");
    Output_VariableList_mm.write("# v_pt							25			200		700			pT^{W}_{l}_(GeV)\n");
    Output_VariableList_mm.write("v_mt								20			0		400			mT^{W}_{l}_(GeV)\n");
    Output_VariableList_mm.write("# pfMET							28			0		560			MET[GeV]\n");
    Output_VariableList_mm.write("# pfMETpuppi						28			0		560			MET[GeV]\n");
    Output_VariableList_mm.write("# pfMET							50      0       1000      MET[GeV]\n");
    Output_VariableList_mm.write("# pfMETpuppi_Phi					20   -3.15      3.15	  #phi_{Puppi MET}\n");
    Output_VariableList_mm.write("# pfMET_Phi						30   -3.14      3.14	  #phi_{MET}\n");
    Output_VariableList_mm.write("# ungroomed_jet_pt				32    100       740      pT^{AK8}_(GeV)\n");
    Output_VariableList_mm.write("# ungroomed_jet_eta				25    -2.5      2.5        #eta^{AK8}\n");
    Output_VariableList_mm.write("# ungroomed_jet_phi				30    -3.14     3.14     #phi^{AK8}\n");
    Output_VariableList_mm.write("# ungroomed_PuppiAK8_jet_pt		32    100       740      pT^{puppi AK8}_(GeV)\n");
    Output_VariableList_mm.write("# ungroomed_PuppiAK8_jet_eta		25    -2.5      2.5        #eta^{puppi AK8}\n");
    Output_VariableList_mm.write("# ungroomed_PuppiAK8_jet_phi		30    -3.14     3.14     #phi^{puppi AK8}\n");
    Output_VariableList_mm.write("# jet_mass_pr						15      65       105     Jet_Pruned_Mass_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# jet_mass_pr						22      40       150    Jet_Pruned_Mass_(GeV)\n");
    Output_VariableList_mm.write("# jet_mass_so						22      40       150    Jet_Softdrop_Mass_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# PuppiAK8_jet_mass_pr			22      40       150    puppiAK8_Jet_Pruned_Mass_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# PuppiAK8_jet_mass_so			22      40       150    puppiAK8_Jet_Softdrop_Mass_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# mass_lvj_type0					40    0       3000    M_{WW}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# mass_lvj_type2					40    0       3000    M_{WW}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# mass_lvj_type0_PuppiAK8			40    0       3000    M_{WW}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# mass_lvj_type2_PuppiAK8			40    0       3000    M_{WW}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# mass_lvj_type0_met				56    200       3000    M_{WW}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# mass_lvj_type2_met				56    200       3000    M_{WW}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("# nu_pz_type0						30   -500      500        pZ^{#nu}[GeV]\n");
    Output_VariableList_mm.write("# nu_pz_type2                   30   -500      500        pZ^{#nu}[GeV]\n");
    Output_VariableList_mm.write("#nu_pz_type0_met                   60   -500      500        pZ^{#nu}[GeV]\n");
    Output_VariableList_mm.write("#nu_pz_type2_met                   60   -500      500        pZ^{#nu}[GeV]\n");
    Output_VariableList_mm.write("#nbjets_csvl_veto                  5      0       5        N_{bjet}^{csvl}\n");
    Output_VariableList_mm.write("#nbjets_csvm_veto                  5      0       5        N_{bjet}^{csvm}\n");
    Output_VariableList_mm.write("#nbjets_csvt_veto                  5      0       5        N_{bjet}^{csvt}\n");
    Output_VariableList_mm.write("#numberJetBin                      5      0       5        N_{jets}\n");
    Output_VariableList_mm.write("# jet_tau2tau1                     25     0.      1.       #tau_{2}/#tau_{1}\n");
    Output_VariableList_mm.write("# PuppiAK8_jet_tau2tau1                     25     0.      1.       puppiAK8_#tau_{2}/#tau_{1}\n");
    Output_VariableList_mm.write("# jet2_pt				 25	0	500	 pT^{AK4}_{1}_(GeV)\n");
    Output_VariableList_mm.write("# jet2_btag				 25	0	1	 btag^{AK4}_{1}_(GeV)\n");
    Output_VariableList_mm.write("# jet3_pt				 25	0	500	 pT^{AK4}_{1}_(GeV)\n");
    Output_VariableList_mm.write("# jet3_btag				 25	0	1	 btag^{AK4}_{1}_(GeV)\n");
    Output_VariableList_mm.write("#jet_tau2tau1_exkT                30     0.1      1.      #tau_{2}/#tau_{1}_extkT\n");
    Output_VariableList_mm.write("#jet_tau2tau1_pr                  30     0.1      1.      #tau_{2}/#tau_{1}_pruned\n");
    Output_VariableList_mm.write("#jet_massdrop_pr                  35     0.1      1.      Pruned_Mass_Drop_(GeV/c^{2})\n");
    Output_VariableList_mm.write("#jet_qjetvol                      35     0        1.         QjetVolatility\n");
    Output_VariableList_mm.write("#jet_charge                       50     -2.5     2.5       Jet_Charge\n");
    Output_VariableList_mm.write("#jet_charge_k05                       50     -2.0     2.0       Jet_Charge_k05\n");
    Output_VariableList_mm.write("#jet_charge_k07                       50     -1.5     1.5       Jet_Charge_k07\n");
    Output_VariableList_mm.write("#jet_charge_k10                       45     -0.8     0.8       Jet_Charge_k10\n");
    Output_VariableList_mm.write("#jet_GeneralizedECF                   35     0.     0.6       Jet_Generalized_ECF\n");
    Output_VariableList_mm.write("#ttb_ca8_mass_pr                   25      40      130      Jet_Pruned_Mass_(GeV/c^{2})\n");
    Output_VariableList_mm.write("#ttb_ht                            35     0      700       \n");
    Output_VariableList_mm.write("#ttb_ca8_ungroomed_pt              30   200      600       pT^{AK8}[GeV]\n");
    Output_VariableList_mm.write("#ttb_ca8_tau2tau1_exkT             30    0.1      1.0      #tau_{2}/#tau_{1}_exkT\n");
    Output_VariableList_mm.write("#ttb_ca8_tau2tau1_pr               30    0.1      1.0       #tau_{2}/#tau_{1}_pruned\n");
    Output_VariableList_mm.write("#ttb_ca8_tau2tau1                  30    0.1      1.0       #tau_{2}/#tau_{1}\n");
    Output_VariableList_mm.write("#ttb_ca8_charge                    50     -2.5     2.5       Jet_Charge\n");
    Output_VariableList_mm.write("#ttb_ca8_charge_k05                50     -2.0     2.0       Jet_Charge_k05\n");
    Output_VariableList_mm.write("#ttb_ca8_charge_k07                50     -1.5     1.5       Jet_Charge_k07\n");
    Output_VariableList_mm.write("#ttb_ca8_charge_k10                45     -0.8     0.8       Jet_Charge_k10\n");
    Output_VariableList_mm.write("#ttb_ca8_GeneralizedECF            30     0.     0.5       Jet_Generalized_ECF\n");
    Output_VariableList_mm.write("#ttb_ca8_mu                        20    0.1      0.7       Pruned_Mass_Drop_(GeV/c^{2})\n");
    Output_VariableList_mm.write("#ttb_mlvj                          40   400     1400       M_{WW}(GeV/c^{2})\n");
    Output_VariableList_mm.write("#vbf_maxpt_j1_bDiscriminatorCSV  50      0       1          j1_bDiscriminator\n");
    Output_VariableList_mm.write("vbf_maxpt_j1_eta                50      -5      5          #eta_{j1}\n");
    Output_VariableList_mm.write("vbf_maxpt_j1_pt                 50      0       300        pT_{j1}_(GeV)\n");
    Output_VariableList_mm.write("#vbf_maxpt_j1_QGLikelihood       50      0       1          j1_QGLikelihood\n");
    Output_VariableList_mm.write("#vbf_maxpt_j2_bDiscriminatorCSV  50      0       1          j2_bDiscriminator\n");
    Output_VariableList_mm.write("vbf_maxpt_j2_eta                50      -5      5          #eta_{j2}\n");
    Output_VariableList_mm.write("vbf_maxpt_j2_pt                 50      0       300        pT_{j2}_(GeV)\n");
    Output_VariableList_mm.write("#vbf_maxpt_j2_QGLikelihood       50      0       1          j2_QGLikelihood\n");
    Output_VariableList_mm.write("abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta) 35    0       9     #Delta#eta_{jj}\n");
    Output_VariableList_mm.write("vbf_maxpt_jj_m                  40      0       1500        M_{jj}_(GeV/c^{2})\n");
    Output_VariableList_mm.write("vbf_maxpt_jj_phi                50      -3.14   3.14       #phi_{jj}\n");
    Output_VariableList_mm.write("vbf_maxpt_jj_eta                30      -4.7       4.7     #eta_{jj}\n");
    Output_VariableList_mm.write("#mass_ungroomedjet_closerjet      30      80     400         M_{top}^{had}\n");
    Output_VariableList_mm.write("#mass_leptonic_closerjet          30      100     400   	    M_{top}^{lep}\n");
    Output_VariableList_mm.close();
    # Make InputFile and SampleListFile
    Ntuple_mm=options.ntuple;       
    Sample_Counter=0;
    Channel_mm=options.channel;
            
    for line in sampleValue:
        Sample_Counter=Sample_Counter+1;

    
    
    # Make SampleListFile     
    Sample_Total_Number_mm=Sample_Counter;
    Sample_Counter_mm=i=j=k=0;
    Cfg_FileName_Table_mm=[[[0 for k in range(2*Sample_Total_Number_mm)]for j in range(Cuts_Total_Number_mm)]for i in range(2)]
    for Sample_Counter_mm in range(Sample_Total_Number_mm):
    
        Sample_mm=sampleValue[Sample_Counter_mm][0];
        Mass_mm=sampleValue[Sample_Counter_mm][2];
        Mass_str_mm=str("%.0f"%Mass_mm);
            
        FileName_Sample_mm="cfg/DataMCComparison_InputCfgFile/MATTEO_SampleList_76x_%s_%s%s.txt"%(Ntuple_mm,Sample_mm,Mass_str_mm);
    
    
        Output_SampleFile_mm=open(FileName_Sample_mm,'w+');
        Output_SampleFile_mm.write("####################################################################################################\n");
        Output_SampleFile_mm.write("###   Sample Name				Reduced Name	Color		XSec (pb)		NumEntriesBefore\n");
        Output_SampleFile_mm.write("####################################################################################################\n");
        Output_SampleFile_mm.write("WWTree_data_golden_2p1			DATA			1			1				1\n");
        Output_SampleFile_mm.write("#WWTree_WJets					W+Jets			2			61526.7			24156124\n");
        Output_SampleFile_mm.write("WWTree_WJets100					W+Jets          2			1347			10205377\n");
        Output_SampleFile_mm.write("WWTree_WJets200					W+Jets          2			360.			4949568\n");
        Output_SampleFile_mm.write("WWTree_WJets400					W+Jets          2			48.9			1943664\n");
        Output_SampleFile_mm.write("#WWTree_WJets600					W+Jets          2			18.77			1041358\n");
        Output_SampleFile_mm.write("WWTree_WJets600bis				W+Jets          2			12.8			3767766\n");
        Output_SampleFile_mm.write("WWTree_WJets800					W+Jets          2			5.26			1568277\n");    
        Output_SampleFile_mm.write("WWTree_WJets1200				W+Jets          2			1.33			246239\n");    
        Output_SampleFile_mm.write("WWTree_WJets2500				W+Jets          2			0.03089			251982\n");    
        Output_SampleFile_mm.write("WWTree_WW						WW              4			118.7			988418\n");
        Output_SampleFile_mm.write("WWTree_WZ						WZ              4			47.13			1000000\n");
        Output_SampleFile_mm.write("WWTree_ZZ						ZZ              4			16.523			985600\n");
        Output_SampleFile_mm.write("WWTree_WW_excl					WW              4			49.997			1924400\n");
        Output_SampleFile_mm.write("WWTree_WZ_excl					WZ              4			10.71			25704656\n");
        Output_SampleFile_mm.write("WWTree_ZZ_excl					ZZ              4			3.22			15301695\n");
        Output_SampleFile_mm.write("WWTree_sch						STop            7			3.65792			998400\n");
        Output_SampleFile_mm.write("WWTree_tch_bar					STop            7			26.0659			1630900\n");
        Output_SampleFile_mm.write("WWTree_tch						STop            7			43.79844		3299200\n");
        Output_SampleFile_mm.write("WWTree_tWch						STop            7			35.6			1000000\n");
        Output_SampleFile_mm.write("WWTree_tWch_bar					STop            7			35.6			999400\n");
        Output_SampleFile_mm.write("#WWTree_TTbar_amcatnlo			tt_bar          210			831.76			38475776\n");
        Output_SampleFile_mm.write("#WWTree_TTbar_madgraph			tt_bar          210			831.76			10215131\n");
        Output_SampleFile_mm.write("WWTree_TTbar					tt_bar          210			831.76			187626200\n");
        Output_SampleFile_mm.write("#WWTree_Higgs750				ggHx750	      	1			0.6398			96200\n");
        Output_SampleFile_mm.write("#WWTree_VBFHiggs750				qqHx750	      	12			0.1915			100000\n");
        Output_SampleFile_mm.write("#WWTree_Higgs1000				ggHx1000	    1			0.1233			400000\n");
        Output_SampleFile_mm.write("#WWTree_VBFHiggs1000			qqHx1000	    12			0.08732			399634\n");
        Output_SampleFile_mm.write("#WWTree_RSGraviton800			RSGrav800GeV    1			1.16691			31906\n");
        Output_SampleFile_mm.write("#WWTree_BulkGraviton1400		RSGrav1000      1			0.000045798		50000\n");
        Output_SampleFile_mm.write("#WWTree_WprimeToWZ1400			Wprime1000      1			0.031283		50000\n");
        Output_SampleFile_mm.write("####################################################################################\n");
        Output_SampleFile_mm.write("######### SIGNAL VALUES\n");
        Output_SampleFile_mm.write("####################################################################################\n");
        
        
        
        if Sample_mm=="BulkGraviton":
           SignalName_mm="WWTree_BulkGraviton"+Mass_str_mm;
           SignalName_VBF_mm="WWTree_VBFBulkGraviton"+Mass_str_mm;
                
        else:
           SignalName_mm="WWTree_Higgs"+Mass_str_mm;
           SignalName_VBF_mm="WWTree_VBFHiggs"+Mass_str_mm;
        
        ReducedName_mm=sampleValue[Sample_Counter_mm][1]
        ReducedName_VBF_mm=sampleValue_VBF[Sample_Counter_mm][1]
        
        Color_VBF_mm=1;
        Color_mm=13;
        
        Xsec_mm=sampleValue[Sample_Counter_mm][3];
        Xsec_VBF_mm=sampleValue_VBF[Sample_Counter_mm][3];
        
        NumberEntriesBefore_mm=sampleValue[Sample_Counter_mm][4];
        NumberEntriesBefore_VBF_mm=sampleValue_VBF[Sample_Counter_mm][4];
        
        ScaleFactor_mm=sampleValue[Sample_Counter_mm][5];
        
        Output_SampleFile_mm.write("%s\t\t\t\t\t%s\t\t\t%.0f\t\t%f\t\t%f\n"%(SignalName_mm,ReducedName_mm,Color_mm,Xsec_mm,NumberEntriesBefore_mm));
        Output_SampleFile_mm.write("%s\t\t\t\t\t%s\t\t\t%.0f\t\t%f\t\t%f\n"%(SignalName_VBF_mm,ReducedName_VBF_mm,Color_VBF_mm,Xsec_VBF_mm,NumberEntriesBefore_VBF_mm));
    
        Output_SampleFile_mm.close();
        
        print "\n\nMade SampleList File:\t%s\n"%FileName_Sample_mm
        Output_Summary_File_mm.write("\n\nMade SampleList File:\t%s\n"%FileName_Sample_mm);
       
        if Channel_mm=="mu":
           leptonT="muon";
        
        elif Channel_mm=="el":
           leptonT="electron";   
        
        else:
           leptonT="emu";
    
        
        if options.nodata:
           withData="true";
        else:
           withData="false";

    
        Cut_Number_mm=0;
        CutType_mm=0;
        for CutType_mm in range(2):
          for Cut_Number_mm in range(Cuts_Total_Number_mm):
          
            
            Cuts_FileName_mm=Cuts_filename_table[CutType_mm][Cut_Number_mm];
            Cut_Type_String_mm=str("%.0f"%(CutType_mm+1));
            Cut_Number_String_mm=str("%.0f"%(Cut_Number_mm+1));
            Dir_Data_Saved_mm="output/run2/MCDATAComparisonPlot_%s_%s_%s%s_%s_%s/"%(Channel_mm,Ntuple_mm,Sample_mm,Mass_mm,Cut_Type_String_mm,Cut_Number_String_mm);
       
            Cfg_Input_FileName_mm="cfg/DataMCComparison_InputCfgFile/MATTEO_DataMCComparison_InputCfgFile_%s%s_%s_%s.cfg"%(Sample_mm,Mass_str_mm,Cut_Type_String_mm,Cut_Number_String_mm)
            Output_SampleFile_mm_sample=open(Cfg_Input_FileName_mm,'w+');
            Output_SampleFile_mm_sample.write("[Input]\n\n");
            Output_SampleFile_mm_sample.write(("InputDirectory = /afs/cern.ch/user/l/lbrianza/work/public/%s/WWTree_%s/\n")%(Ntuple_mm,Channel_mm));
            Output_SampleFile_mm_sample.write("TreeName = otree\n");
            Output_SampleFile_mm_sample.write(("LeptonType = %s\n")%leptonT);
            Output_SampleFile_mm_sample.write(("InputSampleList = %s\n")%(FileName_Sample_mm));
            Output_SampleFile_mm_sample.write("InputVariableList = %s\n"%FileName_VariableList_mm);
            Output_SampleFile_mm_sample.write(("InputCutList = %s\n")%Cuts_FileName_mm);
            Output_SampleFile_mm_sample.write(("SignalqqHName = %s\n")%ReducedName_VBF_mm);
            Output_SampleFile_mm_sample.write(("SignalggHName = %s\n")%ReducedName_mm);
            Output_SampleFile_mm_sample.write(("WithoutData = %s\n\n\n")%withData);    
            Output_SampleFile_mm_sample.write("[Option]\n\n");    
            Output_SampleFile_mm_sample.write("BackgroundWeight = genWeight\n");
            Output_SampleFile_mm_sample.write("BackgroundWeightMCatNLO = 1\n");
            Output_SampleFile_mm_sample.write("SignalggHWeight = 1\n");
            Output_SampleFile_mm_sample.write("SignalqqHWeight = 1\n");
            Output_SampleFile_mm_sample.write("SignalGravitonWeight = genWeight\n");
            Output_SampleFile_mm_sample.write(("Lumi = %f\n")%Lumi_mm);
            Output_SampleFile_mm_sample.write("ttbarControlplots = false\n");
            Output_SampleFile_mm_sample.write("SignalScaleFactor = %f\n"%(ScaleFactor_mm));
            Output_SampleFile_mm_sample.write("NormalizeSignalToData = false\n");
            Output_SampleFile_mm_sample.write("NormalizeBackgroundToData = false\n\n\n");
            Output_SampleFile_mm_sample.write("[Output]\n\n");
            Output_SampleFile_mm_sample.write(("OutputRootDirectory = %s\n"%(Dir_Data_Saved_mm)));
            Output_SampleFile_mm_sample.write("OutputRootFile = Run2_MCDataComparisonRSGraviton2000_mu.root\n");
            Output_SampleFile_mm_sample.write("\n");
            Output_SampleFile_mm_sample.close();
            
            print "\n\nMade CFG File:\t%s\n"%Cfg_Input_FileName_mm
            Output_Summary_File_mm.write("\n\nMade CFG File:\t%s\n"%Cfg_Input_FileName_mm);
            Cfg_FileName_Table_mm[CutType_mm][Cut_Number_mm][Sample_Counter_mm]=Cfg_Input_FileName_mm;
            Cfg_FileName_Table_mm[CutType_mm][Cut_Number_mm][Sample_Total_Number_mm+Sample_Counter_mm]=Dir_Data_Saved_mm;
    
    
    print "\n\n----------- Check CFG name ----------------\n"
    i=j=k=0;
    for i in range(2):
        for j in range(Cuts_Total_Number_mm):
            for k in range(Sample_Total_Number_mm):
                print "\nInputFilename: %s\t Directory: %s\n"%(Cfg_FileName_Table_mm[i][j][k],Cfg_FileName_Table_mm[i][j][Sample_Total_Number_mm+k])
    print "\n\n-------------------------------------------\n"
    i=j=k=0;
    Significance_Table_mm=[[[[0 for z in range(number_Events_type+4) ]for k in range(Sample_Total_Number_mm)]for j in range(Cuts_Total_Number_mm)]for i in range(2)];
    
    
    
    #########################################################
    ######### MAKE CONTROL PLOTS
    #########################################################
    Cut_Number_mm=0
    for Cut_Number_mm in range(Cuts_Total_Number_mm):
        
        # Make directory
        Plus_Cut_Counter=Cut_Number_mm+1;
        Plus_Cut_Counter_str=str("%.0f"%(Plus_Cut_Counter));
        control_cuts1_dir=Control_Plots_Dir_mm+"/Consecutive_Cuts_%s"%(str(Plus_Cut_Counter));
        if not os.path.isdir(control_cuts1_dir):
               pd5 = subprocess.Popen(['mkdir',control_cuts1_dir]);
               pd5.wait();
    
        control_cuts2_dir=Control_Plots_Dir_mm+"/Single_Cuts_%s"%(str(Plus_Cut_Counter));
        if not os.path.isdir(control_cuts2_dir):
               pd6 = subprocess.Popen(['mkdir',control_cuts2_dir]);
               pd6.wait();
           
        cuts_file1=Cuts_filename_table[0][Cut_Number_mm];
        cuts_file2=Cuts_filename_table[1][Cut_Number_mm];
    
        print "\n\n\n----------------------Making Control Plots-----------------------\n"
        print "\nProcessing Cuts %0.f of %0.f\n"%(Plus_Cut_Counter,Cuts_Total_Number_mm)
        print "Using Cuts File 1:\t%s"%cuts_file1
        print "Using Cuts File 2:\t%s"%cuts_file2
        print "\n-----------------------------------------------------------------\n\n\n"
        
        
        Output_Summary_File_mm.write("\n\n\n----------------------Making Control Plots-----------------------\n");
        Output_Summary_File_mm.write("\nProcessing Cuts %0.f of %0.f\n"%(Plus_Cut_Counter,Cuts_Total_Number_mm));
        Output_Summary_File_mm.write("Using Cuts File 1:\t%s"%cuts_file1);
        Output_Summary_File_mm.write("Using Cuts File 2:\t%s"%cuts_file2);
        Output_Summary_File_mm.write("\n-----------------------------------------------------------------\n\n\n");
        
        
        
        
        
        
        
        
        # Run ControPlots code
        n_sample=0;
        for n_sample in range(Sample_Total_Number_mm):
            

            
            sample=sampleValue[n_sample][0];
            mass=sampleValue[n_sample][2];
            mass_str=str("%.0f"%mass);       
            
            Output_Summary_File_mm.write("\n\n\n\n\n\n\n\n-----------------------------------------------------------------------------------------------------------------\n");
            Output_Summary_File_mm.write("PROCESSING:\t%s%s\n"%(sample,mass_str));
            #Output_Summary_File_mm.write("\n-----------------------------------------------------------------\n\n\n");
            print "\n\n\n\n\n-----------------------------------------------------------------------------------------------------------------\n"
            print "PROCESSING:\t%s%s\n"%(sample,mass_str)
                        
            
            final_dir1=control_cuts1_dir+"/"+sample+mass_str;
            if not os.path.isdir(final_dir1):
               pd7 = subprocess.Popen(['mkdir',final_dir1]);
               pd7.wait();
                 
            final_dir2=control_cuts2_dir+"/"+sample+mass_str;
            if not os.path.isdir(final_dir2):
               pd8 = subprocess.Popen(['mkdir',final_dir2]);
               pd8.wait();
            
            print "\nFinalDir1: %s"%final_dir1
            print "\nFinalDir2: %s"%final_dir2
            print "\n-----------------------------------------------------------------------------------------------------------------\n\n\n\n"
            
                            
            Output_Summary_File_mm.write("\n-----------------------------------------------------------------\n");
            Output_Summary_File_mm.write("CONSECUTIVE CUTS");
            Output_Summary_File_mm.write("\n-----------------------------------------------------------------\n");    
            Output_Summary_Latex_File_mm.write("\nCONSECUTIVE CUT\n");
            cfg_file_1=[Cfg_FileName_Table_mm[0][Cut_Number_mm][n_sample],Cfg_FileName_Table_mm[0][Cut_Number_mm][Sample_Total_Number_mm+n_sample]];
            #def run_log(sampleNumber,cutNumber,ctype,total_cutNumber,Significance_Table_mm_log,input_cfg_file,latex_file,Output_Summary_File_mm,in_cc_dir):
            Significance_Table_mm=run_log(n_sample,Cut_Number_mm,1,Cuts_Total_Number_mm,Significance_Table_mm,cfg_file_1,Output_Summary_Latex_File_mm,Output_Summary_File_mm,final_dir1)
            
            
            
            
            
            Output_Summary_File_mm.write("\n-----------------------------------------------------------------\n");
            Output_Summary_File_mm.write("SINGLE CUT");
            Output_Summary_File_mm.write("\n-----------------------------------------------------------------\n");  
            Output_Summary_Latex_File_mm.write("\nSINGLE CUT\n");
            cfg_file_2=[Cfg_FileName_Table_mm[1][Cut_Number_mm][n_sample],Cfg_FileName_Table_mm[1][Cut_Number_mm][Sample_Total_Number_mm+n_sample]];
            #def run_log(sampleNumber,cutNumber,ctype,total_cutNumber,Significance_Table_mm_log,input_cfg_file,latex_file,Output_Summary_File_mm,in_cc_dir):
            Significance_Table_mm=run_log(n_sample,Cut_Number_mm,2,Cuts_Total_Number_mm,Significance_Table_mm,cfg_file_2,Output_Summary_Latex_File_mm,Output_Summary_File_mm,final_dir2)
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #########################################################
    ######### MAKE BEAMER OUTPUT
    #########################################################
    
    # Beamer Settings
    Output_Beamer_Latex_File_mm.write("\documentclass{beamer}\n");
    Output_Beamer_Latex_File_mm.write("\usetheme{Boadilla}\n");
    Output_Beamer_Latex_File_mm.write("\usecolortheme{seahorse}\n");
    Output_Beamer_Latex_File_mm.write("\\title{ControlPlots}\n");
    Output_Beamer_Latex_File_mm.write("\\author{Matteo Rappo}\n");
    Output_Beamer_Latex_File_mm.write("\setbeamertemplate{navigation symbols}{}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage[latin1]{inputenc}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage[english,italian]{babel}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{amsmath}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{enumerate}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{amsfonts}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{amssymb}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{float}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{placeins}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{subfig}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{multirow,makecell}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{array,booktabs}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{comment}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{scrextend}\n");
    Output_Beamer_Latex_File_mm.write("\usepackage{verbatim,longtable}\n");
    Output_Beamer_Latex_File_mm.write("\setbeamertemplate{caption}[numbered]\n");
    Output_Beamer_Latex_File_mm.write("\\newcolumntype{P}[1]{>{\centering\\arraybackslash}p{#1}}\n");
    Output_Beamer_Latex_File_mm.write("\\newcolumntype{M}[1]{>{\centering\\arraybackslash}m{#1}}\n");
    Output_Beamer_Latex_File_mm.write("\\newcolumntype{D}[1]{>{\\arraybackslash}m{#1}}\n");
    Output_Beamer_Latex_File_mm.write("\\newcolumntype{C}[1]{>{\centering\let\\newline\\\\arraybackslash\hspace{0pt}}m{#1}}\n");
    Output_Beamer_Latex_File_mm.write("\n");
    Output_Beamer_Latex_File_mm.write("\n");
    Output_Beamer_Latex_File_mm.write("\n");
    Output_Beamer_Latex_File_mm.write("\changefontsizes{9pt}\n");
    Output_Beamer_Latex_File_mm.write("\\begin{document}\n");
    Output_Beamer_Latex_File_mm.write("\n");
    Output_Beamer_Latex_File_mm.write("\n");
     


    # Ntuple Name in \texttt{} environment
    tmp_ntuple_texttt=replace_latex(options.ntuple);
    Ntuple_Name_texttt="\\texttt{"+tmp_ntuple_texttt+"}";
    
   

    # Slides with Cuts
    
    ## Consecutive Cuts
    Output_Beamer_Latex_File_mm.write("\n\n\n");
    Output_Beamer_Latex_File_mm.write("\changefontsizes{9pt}\n");
    Output_Beamer_Latex_File_mm.write("\\begin{frame}[t,allowframebreaks]\n");
    Output_Beamer_Latex_File_mm.write("\\frametitle{Control Plots - Consecutive Cuts }\n");   
    Output_Beamer_Latex_File_mm.write("\\framesubtitle{$\mu$-channel  Ntuple: "+Ntuple_Name_texttt+" }\n");
    Output_Beamer_Latex_File_mm.write("\changefontsizes{7pt}\n");
    #Output_Beamer_Latex_File_mm.write("\\begin{table}[H]\n");
    #Output_Beamer_Latex_File_mm.write("\\begin{center}\n");
    Output_Beamer_Latex_File_mm.write("\\begin{longtable}{|M{10pt}|D{310pt}|}\n");
    Output_Beamer_Latex_File_mm.write("\hline \multicolumn{2}{|c|}{Consecutive Cuts} \\\ \n");
    
    
    j=0;
    for j in range(Cuts_Total_Number_mm):
        cn=j+1;
        tmp=replace_latex(cuts_table_main[0][j]);
        #tmp2="\\texttt{"+tmp1+"}";
        Output_Beamer_Latex_File_mm.write("\hline %.0f & %s \\\ \n"%(cn,tmp));
         
            
    Output_Beamer_Latex_File_mm.write("\hline\n");
    #Output_Beamer_Latex_File_mm.write("\end{tabular}\n");
    #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
    #Output_Beamer_Latex_File_mm.write("\end{center}\n");
    Output_Beamer_Latex_File_mm.write("\end{longtable}\n");
    Output_Beamer_Latex_File_mm.write("\end{frame}\n");
    
    ## Single Cut
    Output_Beamer_Latex_File_mm.write("\n\n\n");
    Output_Beamer_Latex_File_mm.write("\changefontsizes{9pt}\n");
    Output_Beamer_Latex_File_mm.write("\\begin{frame}\n");
    Output_Beamer_Latex_File_mm.write("\\frametitle{Control Plots - Single Cut }\n");   
    Output_Beamer_Latex_File_mm.write("\\framesubtitle{$\mu$-channel  Ntuple: "+Ntuple_Name_texttt+" }\n");
    Output_Beamer_Latex_File_mm.write("\changefontsizes{7pt}\n");
    Output_Beamer_Latex_File_mm.write("\\begin{table}[H]\n");
    Output_Beamer_Latex_File_mm.write("\\begin{center}\n");
    Output_Beamer_Latex_File_mm.write("\\begin{tabular}{|M{10pt}|D{310pt}|}\n");
    Output_Beamer_Latex_File_mm.write("\hline \multicolumn{2}{|c|}{Single Cut} \\\ \n");
    
    j=0;
    for j in range(Cuts_Total_Number_mm):
        cn=j+1;
        tmp=replace_latex(cuts_table_main[1][j]);
        #tmp2="\\texttt{"+tmp1+"}";
        Output_Beamer_Latex_File_mm.write("\hline %.0f & %s \\\ \n"%(cn,tmp));
          
            
    
    Output_Beamer_Latex_File_mm.write("\hline\n");
    Output_Beamer_Latex_File_mm.write("\end{tabular}\n");
    #ofile.write("\\caption{Significanza in funzione delle Input Variables for Cut Optimization per diversi valori di risonanza}\n");
    Output_Beamer_Latex_File_mm.write("\end{center}\n");
    Output_Beamer_Latex_File_mm.write("\end{table}\n");
    Output_Beamer_Latex_File_mm.write("\end{frame}\n");




            


    
    
    # Significance slides
    nsample=0;
    for nsample in range(Sample_Total_Number_mm):
        
        # Consecutive Cuts
        make_latex_table(Cuts_Total_Number_mm,Significance_Table_mm,1,nsample,Output_Beamer_Latex_File_mm,Ntuple_Name_texttt);
        
        # Single Cut
        make_latex_table(Cuts_Total_Number_mm,Significance_Table_mm,2,nsample,Output_Beamer_Latex_File_mm,Ntuple_Name_texttt);
        
    Output_Beamer_Latex_File_mm.write("\end{document}\n");
    Output_Summary_Latex_File_mm.close()
    Output_Beamer_Latex_File_mm.close();  
    Output_Summary_File_mm.close();
        
        
         
        
        
      
            
            
            
            
            
    
