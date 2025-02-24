import sys
from time import sleep as sl
import glob
import math,statistics

infolder = sys.argv[1]# "/Users/z3521839/Library/CloudStorage/OneDrive-UNSW/PycharmProjects/Lab_utils/treelab_pranita/inputs"
# infolder = "/Users/mjohnpayne/Library/CloudStorage/OneDrive-UNSW/PycharmProjects/Lab_utils/treelab_pranita/inputs"
# of = "/Users/mjohnpayne/Library/CloudStorage/OneDrive-UNSW/PycharmProjects/Lab_utils/treelab_pranita/testout.txt"
of = sys.argv[2]
outfile = open(of,"w")
outfile.write("Strain\tshig/eiec\tserotype\tEcOnc10_relative_cov\tgene_based_pathotypes\tgenes_present\n")


inshiga = glob.glob(infolder+"/*_shig_out")
instec = glob.glob(infolder+"/*_stec_out")
inkma = glob.glob(infolder+"/*_kma_out.res")
shigaids = set([x.split("/")[-1].replace("_shig_out","") for x in inshiga])

stecids = set([x.split("/")[-1].replace("_stec_out","") for x in instec])

kmaids = set([x.split("/")[-1].replace("_kma_out.res","") for x in inkma])

inall = shigaids.intersection(stecids,kmaids)

if inall==shigaids==stecids==kmaids:
    print("all present")
else:
    print("common",len(inall))
    print("shig",len(shigaids))
    print("stec",len(stecids))
    print("kma",len(kmaids))


housekeeping = ["recA","purA","mdh","icd","gyrB","fumC","adk"]
pathotypeing = ["eae","escV","stx1","stx2","bfpB","aggR","aaiC","daaD","afaD","fyuA","traT","hlyA","papC","elt","est","invA","ipaH"]

for strain in inall:
    stec = False
    shiga = False
    srna = False
    type = ""
    srna_cov = 0
    relsrnacov = 0
    shigout = f"{infolder}/{strain}_shig_out"
    stecout = f"{infolder}/{strain}_stec_out"
    kmaout = f"{infolder}/{strain}_kma_out.res"
    kmares = open(kmaout, "r").read().splitlines()
    kmatypels = []
    kmagenels = []
    for line in kmares:
        col=line.split("\t")
        depthls = []

        if col[0] in housekeeping:
            depth = col[8].replace(" ","")
            depthls.append(float(depth))
        if col[0] == "EcOnc10":
            srna = True
            srna_cov = float(col[8].replace(" ",""))
        if col[0] in pathotypeing:
            kmagenels.append(col[0])
    kmagenels = list(set(kmagenels))
    genomedepth = statistics.mean(depthls)
    if srna:
        relsrnacov = srna_cov/genomedepth
    shigres = open(shigout,"r").read().splitlines()
    print(strain)
    print(shigres[:10])
    res = shigres[1].split("\t")
    sero = res[4]
    stecres = open(stecout,"r").read().splitlines()
    stecres = stecres[1].split("\t")

    if ("eae" in kmagenels or "escV" in kmagenels) and ("stx1" in kmagenels or "stx2" in kmagenels) and "bfp" not in kmagenels:
        genetype="EHEC"
    elif ("eae" in kmagenels or "escV" in kmagenels) and ("stx1" not in kmagenels and "stx2" not in  kmagenels):
        genetype="EPEC"
    elif ("aggR" in kmagenels and "aaiC" in kmagenels) and ("stx1" not in kmagenels and "stx2" not in  kmagenels):
        genetype="EAEC"
    elif ("eae" in kmagenels or "escV" in kmagenels) and ("eae" not in kmagenels and "escV" not in kmagenels): #### not either or not both
        genetype="DAEC"
    elif len(set(["fyuA","traT","hlyA","papC"]).intersection(kmagenels)) > 1:
        genetype = "UPEC"
    elif "elt" in kmagenels or "est" in kmagenels:
        genetype = "ETEC"
    elif ("invA" in kmagenels and "ipaH" in kmagenels) and (("stx1" not in kmagenels and "stx2" not in  kmagenels) and "bfp" not in kmagenels):
        genetype="EIEC"
    else:
        genetype="other"




    if res[3] == "Shigella/EIEC Unclustered":
        # print(strain,"EIEC",sero)
        shiga = True
        tooltype = "Shigella/EIEC"
    elif res[4].startswith("S"):
        # print(strain,"shigella",sero)
        shiga = True
        tooltype = "Shigella"
        detail = res[4]
    elif res[3].startswith("C"):
        # print(strain, "EIEC",sero)
        shiga = True
        tooltype = "EIEC"
        detail = res[4]
    elif stecres[1] == "Unclustered STEC":
        stec=True
        # print(strain,stecres[1])
        tooltype = "STEC"
        detail = stecres[3]
    elif stecres[1] == "Other_Ecoli":
        stec = False
        tooltype = "Other"
        detail = stecres[3]
        # print(strain, type,genomedepth,relsrnacov)
    else:
        tooltype = "STEC"
        stec = True
        detail = stecres[3]
    outfile.write("\t".join([strain,tooltype,detail, str(relsrnacov),genetype,"|".join(kmagenels)])+"\n")
outfile.close()











