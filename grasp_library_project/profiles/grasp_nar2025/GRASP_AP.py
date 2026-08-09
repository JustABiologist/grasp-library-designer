import sys
import pandas as pd
import csv
from datetime import datetime
import os
import dnacauldron as dc
#import locale
import _locale
_locale._getdefaultlocale = (lambda *args: ['en_US', 'utf8'])
from Bio import BiopythonWarning
import warnings
warnings.simplefilter('ignore',BiopythonWarning)
#print(locale.getpreferredencoding(False))


def listToString(list_in):  # root function; convert lists to strings for easy printing
    code_out_str = ""
    return code_out_str.join(list_in)


def DNA_to_PPRcode(code):
    # translates DNA or RNA into PPR binding code.
    code = code.upper()
    PPR_code = []
    dnaA = ["A"]
    dnaT = ["T", "U"]
    dnaC = ["C"]
    dnaG = ["G"]
    for x in code:
        if x in dnaA:
            PPR_code.append("TN")
        elif x in dnaT:
            PPR_code.append("ND")
        elif x in dnaC:
            PPR_code.append("NN")
        elif x in dnaG:
            PPR_code.append("TD")
        else:
            sys.exit("Error: noncanonical base detected. Please input ACGTU only.")
    PPR_codenew = listToString(PPR_code)
    return PPR_codenew


def PPR_sizecheck(PPRstr, GAP_ef):
    sizeError = ''
    if GAP_ef == '1':
        return len(PPRstr) if len(PPRstr) == 9 or len(PPRstr) == 14 or len(PPRstr) == 19 \
            else sizeError + "Error: PPR target unsuitable size (must be 9, 14, or 19 bases)"
    elif GAP_ef == '2':
        return len(PPRstr) if len(PPRstr) == 16 or len(PPRstr) == 21 or len(PPRstr) == 26 \
            else sizeError + "Error: PPR target unsuitable size for editing factor (must be 16, 21, or 26 bases)"
    return sizeError


def codify(PPRcode):
    # turns PPR code string into part-searchable variables of a list
    if len(PPRcode) > 14:
        trimcode = PPRcode[1:-1]  # removes first and last index
    else:
        trimcode = PPRcode
    duplex = [trimcode[i:i + 2] for i in range(0, len(trimcode), 2)]  # makes list of pairs of AAs (5th and Last)
    if len(PPRcode) > 14:
        duplex.insert(0, PPRcode[0])
        duplex.append(PPRcode[-1:])
    return duplex


def partpicker(PPRcode):
    code_part = ["T", "N", "DN", "DT", "NN", "NT"]
    code_compute = ["_5T", "_5N", "_LD5N", "_LD5T", "_LN5N", "_LN5T"]
    part_prefix9S = ["pPR-1_1A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_1E",
                     "pPR-1_2A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_2E_L"]
    part_prefix14S = ["pPR-1_1A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_14E",
                      "pPR-1_14A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_1E",
                      "pPR-1_2A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_2E_L"]
    part_prefix19S = ["pPR-1_1A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_14E",
                      "pPR-1_14A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_19E",
                      "pPR-1_19A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_1E",
                      "pPR-1_2A", "pPR-1_B", "pPR-1_C", "pPR-1_D", "pPR-1_2E_L"]
    PPRlen = len(PPRcode)
    lastindex = PPRlen - 1

    if PPRlen == 10:
        for i in range(lastindex):  # all but the very last part are treated here
            PPRcode[i] = part_prefix9S[i] + code_compute[code_part.index(PPRcode[i])]
        if PPRcode[lastindex] == "D":
            PPRcode[lastindex] = part_prefix9S[lastindex] + "D"
        elif PPRcode[lastindex] == "N":
            PPRcode[lastindex] = part_prefix9S[lastindex] + "N"
    elif PPRlen == 15:
        for i in range(lastindex):
            PPRcode[i] = part_prefix14S[i] + code_compute[code_part.index(PPRcode[i])]
        if PPRcode[lastindex] == "D":
            PPRcode[lastindex] = part_prefix14S[lastindex] + "D"
        elif PPRcode[lastindex] == "N":
            PPRcode[lastindex] = part_prefix14S[lastindex] + "N"
    elif PPRlen == 20:
        for i in range(lastindex):
            PPRcode[i] = part_prefix19S[i] + code_compute[code_part.index(PPRcode[i])]
        if PPRcode[lastindex] == "D":
            PPRcode[lastindex] = part_prefix19S[lastindex] + "D"
        elif PPRcode[lastindex] == "N":
            PPRcode[lastindex] = part_prefix19S[lastindex] + "N"

    return PPRcode


def S2prompt():
    global S2choice
    S2choice = str(
        input("The nonspecific S2#12 motif is capable of binding and permitting editing with any nucleotide at -4."
              "Use this S2 motif for all assemblies? [Y/N]"))
    answersy = ['y', 'Y', 'yes', 'YES', 'ye', 'yEs', 'Yes', 'yeS', 'YEs', 'yES', 'YeS', '1']
    answersn = ['n', 'N', 'No', 'no', 'NO', 'nO', '0']

    if S2choice in answersy:
        S2choice = 'Y'
    elif S2choice in answersn:
        S2choice = 'N'
    else:
        print("Invalid response. Choose again.")
        S2prompt()
    return S2choice


def editingpartpicker(RemainderPPRcode):
    Remainderparts = []
    P2parts = {"TN": "P2(A)", "NN": "P2(C)", "TD": "P2(G)", "ND": "P2(U)"}
    if S2choice == 'N':
        S2parts = {"TN": "S2GN(A)", "NN": "S2GG(C)", "TD": "S2AD(G)", "ND": "S2GD(U)"}
        Remainderparts.append(P2parts[RemainderPPRcode[0]] + 'L2' + S2parts[RemainderPPRcode[2]])
        Remainderparts.append("E1E2DYW_PG")
    else:
        Remainderparts.append(P2parts[RemainderPPRcode[0]] + 'L2S2#12(N)')
        Remainderparts.append("E1E2DYW_PG")

    return Remainderparts


def ppr_assort(PPRcode, PPRlen):  # picker. receives just a STRING of residues (PPR binding code, trims to 9 motifs)
    if PPRlen - 7 == 9 or PPRlen - 7 == 14 or PPRlen - 7 == 19:
        BindingPPRcode = codify(PPRcode[0:(PPRlen * 2 - 14)])
        RemainderPPRcode = codify(PPRcode[(PPRlen * 2 - 14):])
        # print(PPRcode)
        # print(BindingPPRcode)
        # print(RemainderPPRcode)
    else:
        PPRcode = codify(PPRcode)
    if PPRlen == 9 or PPRlen == 14 or PPRlen == 19:
        PPRcode = partpicker(PPRcode)
    else:

        PPRcode = partpicker(BindingPPRcode)
        RemainderParts = editingpartpicker(RemainderPPRcode)
        PPRcode.append(RemainderParts)
        # print(PPRcode)
    return PPRcode


def GAP_file(filepath):
    namelist = []
    seqlist = []
    try:
        with open(filepath, 'r') as openfile:
            for i in openfile:
                i = i.rstrip()
                if '>' in i:
                    # print(i)
                    namelist.append(i[1:])
                elif i == '':
                    pass
                elif nt_check(i):
                    # print(i)
                    seqlist.append(i)
                else:
                    pass
    except:
        print('Filepath error, check path spelling.')
        GAP_intro()
        return
    return namelist, seqlist


def nt_check(nt_str):
    nt_str = nt_str.upper()
    # print(nt_str)
    nt = ["A", "C", "G", "T", "U"]
    for i in nt_str:
        if i in nt:
            pass
        else:
            sys.exit("Error: Invalid character present. Revise target sequences and retry.")

    return True


def assembly_reformat(gross_list):
    reformat_gross_list = []
    for assembly in gross_list:
        working_assembly = [assembly[i:i + 5] for i in range(0, len(assembly), 5)]
        reformat_gross_list.append(working_assembly)

    return reformat_gross_list


def textwrite_prompt():
    textwrite = str(input("Write output to a .txt file? [Y/N]: "))
    answersy = ['y', 'Y', 'yes', 'YES', 'ye', 'yEs', 'Yes', 'yeS', 'YEs', 'yES', 'YeS', '1']
    answersn = ['n', 'N', 'No', 'no', 'NO', 'nO', '0']

    if textwrite in answersy:
        textwrite = 'Y'
    elif textwrite in answersn:
        textwrite = 'N'
    else:
        print("Invalid response. Choose again.")
        textwrite_prompt()
    return textwrite


def protocol_printer(namelist, assembly_list):
    textwrite = textwrite_prompt()
    # print(namelist)
    for i in range(len(assembly_list)):
        edit_assembly = ''
        protein = assembly_list[i]
        if 'E1E2DYW_PG' in protein[-1:][0][0]:
            edit_assembly = protein[-1:]
            # print(edit_assembly[0][0])
            protein = protein[:-1]
        #             print("EHLLPO")
        #             print(protein)
        assemblynum = 0
        protlen = len(protein)
        # print('protein:')
        # print(protein)
        # print(protein)
        for assembly in protein:
            # print(assembly)
            assemblynum += 1
            if assemblynum == 1:
                partnum = 'pPR0_1'
            elif assemblynum == protlen:
                partnum = 'pPR0_2'
            elif assemblynum == 3 and protlen == 4:
                partnum = 'pPR0_19'
            elif assemblynum == 2 and protlen == 3 or protlen == 4:
                partnum = 'pPR0_14'

            if len(assembly_list) <= 100:
                print(pd.DataFrame({
                    (namelist[i] + " - " + str(assemblynum) + "/" + str(protlen) +
                     ' - ' + partnum): assembly,
                }))

            if textwrite == 'Y':
                filename = ('assembly_' + datenow + '.txt')
                with open(filename, 'a') as outfile:
                    outfile.write(namelist[i] + " - " +
                                  str(assemblynum) + "/" + str(protlen) + ' - ' + partnum + '\n')
                    for minus1 in assembly:
                        outfile.write((minus1) + '\n')
                    outfile.write('\n')
        if edit_assembly != '':
            if len(assembly_list) <= 100:
                print(pd.DataFrame({
                    (namelist[i] + ' - Parts for editing: '): edit_assembly[0][0],
                }))

        if textwrite == 'Y':
            # filename = ('assembly' + datetime.today().strftime("%Y-%m-%d_%H-%M-%S") + '.txt')
            with open(filename, 'a') as outfile:
                outfile.write(namelist[i] + ' - Parts for editing: \n')
                for minus1 in edit_assembly[0][0]:
                    outfile.write((minus1) + '\n')
                outfile.write('\n')
    if textwrite == 'Y':
        print("Done! Output located in: " + os.getcwd())
    else:
        print("Done!")
    # GAP_intro()
    return


def GAP_editfactor():
    editfactorprompt = str(input("Are you assembling binding PPRs [1] or editing PPRS [2]? "))
    return editfactorprompt


def GAP_intro():
    gross_part_list = []
    print("\n" + "1: To generate an assembly plan for a single target, type/paste a nucleotide sequence." + '\n' +
          "Alternatively, paste the path to a text file (.txt, .fasta, etc). On Windows, this can be done by "
          "holding SHIFT and right-clicking on a file, then copying as path." + '\n' +
          "Alternatively, place file in " + os.getcwd() + " and copy its name." + '\n' +
          "Files with multiple sequences must be loaded in FASTA format." + '\n' +
          "Note that editing and non-editing factors cannot be included in the same file at this time." + '\n' +
          "Ambiguity codes are not supported. Type 'Exit' to exit.")
    GAPstr = str(input("Enter here: "))
    if GAPstr == 'Exit':
        print("Exiting... \n")
        intro()
        return
    GAP_ef = GAP_editfactor()
    if GAP_ef == '1':
        pass
    elif GAP_ef == '2':
        S2prompt()
    else:
        print("Improper selection. '1' or '2' only.")
        GAP_ef = GAP_editfactor()

    if '.fasta' in GAPstr or '.txt' in GAPstr:
        if GAPstr[0] == '"':
            GAPstr = GAPstr[1:-1]
        # print(GAPstr)
        namelist, seqlist = GAP_file(GAPstr)
        tracker = len(seqlist)
        tracknum = 0
        print("Processing...")
        for i in seqlist:
            nt_check(i)
            PPRCODE = DNA_to_PPRcode(i)
            PPRlen = PPR_sizecheck(i, GAP_ef)
            if type(PPRlen) == str:
                print(PPRlen)  # this is an error code. this is also what superior coding looks like
                GAP_intro()
                return
            PPRpartsub = ppr_assort(PPRCODE, PPRlen)
            gross_part_list.append(PPRpartsub)
        assembly_list = assembly_reformat(gross_part_list)
        protocol_printer(namelist, assembly_list)

    else:
        nt_check(GAPstr)
        PPRCODE = DNA_to_PPRcode(GAPstr)
        # print(PPRCODE)
        PPRlen = PPR_sizecheck(GAPstr, GAP_ef)
        if type(PPRlen) == str:
            print(PPRlen)  # this is an error code. this is also what superior coding looks like
            GAP_intro()
            return
        PPRpartlist = ppr_assort(PPRCODE, PPRlen)
        #print(PPRpartlist)
        gross_part_list.append(PPRpartlist)
        assembly_list = assembly_reformat(gross_part_list)
        for protein in assembly_list:
            assemblynum = 0
            # print('protein:')
            # print(protein)
            for assembly in protein:
                if "E1E2DYW_PG" in assembly[0]:
                    assemblynum += 1
                    print(pd.DataFrame({
                        (str(assemblynum) + "/" + str(
                            len(protein))): assembly[0],
                    }))
                else:
                    assemblynum += 1
                    print(pd.DataFrame({
                        (str(assemblynum) + "/" + str(
                            len(protein))): assembly,
                    }))
        print("Done!")
    print("Task completed. Returning to menu.")
    intro()


def DNAcauldron_plan_prepper(sheetname):
    repository = dc.SequenceRepository()
    repository.import_records(folder="GRASP_seqs/", use_file_names_as_ids=True)

    assembly_plan = dc.AssemblyPlan.from_spreadsheet(
        assembly_class=dc.Type2sRestrictionAssembly,
        path=sheetname,
        name=""
    )
    plan_simulation = assembly_plan.simulate(sequence_repository=repository)

    plan_simulation.write_report("GRASP_APO_" + datenow)

    with open("GRASP_APO_" + datenow + "\\readme.txt", 'a') as readme:
        readme.write("Find a copy of all .gb files required for assembly in 'part_records' \n"
                     "Find all level 0 and level 1 assemblies in 'all_construct_records'")

    print("Assembly completed. \n"
          "Outputs: \n" +
          os.getcwd() + "\GRASP_part_list_" + datenow + ".csv -- DNAcauldron-compatible parts list. \n" +
          os.getcwd() + "\GRASP_APO_" + datenow + " -- Assembly plan and constructs.")


def DNA_cauldron_spreadsheet(namelist, proteinlist, makeplan):

    lvl1parts = ['pAGM9121 - T7promTrxA',
                 'pAGM9121 - Term_noTarget', 'modified_1-1R_pICH47802_lc_p15A_ori_',
                 "E1E2DYW_sub_5"]

    sheetname = "GRASP_part_list_" + datenow + ".csv"

    with open(sheetname, 'w', newline='', encoding="utf-8") as infile:

        inwriter = csv.writer(infile, lineterminator='\n')
        inwriter.writerow(["construct","parts",'','','','','','','',''])
        for i in range(len(proteinlist)):
            proteinlist[i][0] += "_AGGT"
            pllen = len(proteinlist[i])
            if pllen >= 10 and pllen < 12: ## for 9-motif binding PPRs
                inwriter.writerow([namelist[i] + "_pPR0_1"] + proteinlist[i][0:5] + ['pAGM9121'])
                inwriter.writerow([namelist[i] + "_pPR0_2"] + proteinlist[i][5:10] + ['pAGM9121'])
                if pllen > 10: ## for 9-motif editing PPRs
                    inwriter.writerow([namelist[i] + "_pPR1"] + [namelist[i] + "_pPR0_1"]+ [namelist[i] + "_pPR0_2"] +
                                      lvl1parts + ["pAGM9121 - " + proteinlist[i][10][0]])
            elif pllen >= 15 and pllen < 17: ## for 14-motif binding PPRs
                inwriter.writerow([namelist[i] + "_pPR0_1"] + proteinlist[i][0:5] + ['pAGM9121'])
                inwriter.writerow([namelist[i] + "_pPR0_14"] + proteinlist[i][5:10] + ['pAGM9121'])
                inwriter.writerow([namelist[i] + "_pPR0_2"] + proteinlist[i][10:15] + ['pAGM9121'])
                if pllen > 15: ## for 9-motif editing PPRs
                    inwriter.writerow([namelist[i] + "_pPR1"] + [namelist[i] + "_pPR0_1"] + [namelist[i] + "_pPR0_14"]
                                      + [namelist[i] + "_pPR0_2"] + lvl1parts + ["pAGM9121 - " + proteinlist[i][15][0]])
            elif pllen >= 20 and pllen < 22: ## for 19-motif binding PPRs
                inwriter.writerow([namelist[i] + "_pPR0_1"] + proteinlist[i][0:5] + ['pAGM9121'])
                inwriter.writerow([namelist[i] + "_pPR0_14"] + proteinlist[i][5:10] + ['pAGM9121'])
                inwriter.writerow([namelist[i] + "_pPR0_19"] + proteinlist[i][10:15] + ['pAGM9121'])
                inwriter.writerow([namelist[i] + "_pPR0_2"] + proteinlist[i][15:20] + ['pAGM9121'])
                if pllen > 20: ## for 19-motif editing PPRs
                    inwriter.writerow([namelist[i] + "_pPR1"] + [namelist[i] + "_pPR0_1"] + [namelist[i] + "_pPR0_14"]
                                      + [namelist[i] + "_pPR0_19"] + [namelist[i] + "_pPR0_2"]
                                      + lvl1parts + ["pAGM9121 - " + proteinlist[i][20][0]])

    if makeplan == 1:
        DNAcauldron_plan_prepper(sheetname)

def GAP_intro_DNAcauldron(makeplan):
    gross_part_list = []
    print("\n" + "1: To generate an assembly plan for a single target, type/paste a nucleotide sequence." + '\n' +
          "Alternatively, paste the path to a text file (.txt, .fasta, etc). On Windows, this can be done by "
          "holding SHIFT and right-clicking on a file, then copying as path." + '\n' +
          "Alternatively, place file in " + os.getcwd() + " and copy its name." + '\n' +
          "Files with multiple sequences must be loaded in FASTA format." + '\n' +
          "Note that editing and non-editing factors cannot be included in the same file at this time." + '\n' +
          "Ambiguity codes are not supported."
          "Type 'Exit' to exit.")
    GAPstr = str(input("Enter here: "))

    if GAPstr == 'Exit':
        print("Exiting... \n")
        intro()
        return
    GAP_ef = GAP_editfactor()
    if GAP_ef == '1':
        pass
    elif GAP_ef == '2':
        S2prompt()
    else:
        print("Improper selection. '1' or '2' only.")
        GAP_ef = GAP_editfactor()

    if '.fasta' in GAPstr or '.txt' in GAPstr:
        if GAPstr[0] == '"':
            GAPstr = GAPstr[1:-1]
        # print(GAPstr)
        namelist, seqlist = GAP_file(GAPstr)
        print("Processing...")
        for i in seqlist:
            nt_check(i)
            PPRCODE = DNA_to_PPRcode(i)
            PPRlen = PPR_sizecheck(i, GAP_ef)
            if type(PPRlen) == str:
                print(PPRlen)  # this is an error code. this is also what superior coding looks like
                GAP_intro()
                return
            PPRpartsub = ppr_assort(PPRCODE, PPRlen)
            gross_part_list.append(PPRpartsub)
        # print(namelist)
        # print(seqlist)
        # print(gross_part_list)
        #makeplan = 1
        DNA_cauldron_spreadsheet(namelist, gross_part_list, makeplan)

    else:
        nt_check(GAPstr)
        PPRCODE = DNA_to_PPRcode(GAPstr)
        # print(PPRCODE)
        PPRlen = PPR_sizecheck(GAPstr, GAP_ef)
        if type(PPRlen) == str:
            print(PPRlen)  # this is an error code. this is also what superior coding looks like
            GAP_intro()
            return
        PPRpartlist = ppr_assort(PPRCODE, PPRlen)
        gross_part_list.append(PPRpartlist)
        #makeplan = 1
        DNA_cauldron_spreadsheet(["your_assembly"], gross_part_list, makeplan)


def tutorial():
    gross_part_list = []
    tuto_bp = "AAUGACGGU"
    bp_split = [i for i in tuto_bp]
    bp_pprcode = DNA_to_PPRcode(tuto_bp)
    bp_pprcode_d = [bp_pprcode[i:i + 2] for i in range(0, len(bp_pprcode), 2)]
    bp_pprcode_codify = codify(bp_pprcode)
    print("This program takes the input of a 9, 14, or 19 nucleotide DNA/RNA string, converts it to the PPR binding "
          "code, and returns a list of level -1 'pPR-1' in blocks of 5, which are used in level 0 assembly of 'ppR0'"
          "plasmids. \nThis also works for targets intended for RNA editing; in these cases, 7 additional nucleotides"
          " (totalling 16, 21, or 26) must be included.")
    print("For example, " + tuto_bp + " is input. This is converted to PPR Binding Code like so:")
    print(' |'.join(bp_split))
    print('|'.join(bp_pprcode_d))
    print("Then, the binding code is rearranged to fit the GRASP Assembly Strategy.")
    print('|'.join(bp_pprcode_codify))
    print("From here, GRASP level -1 motif parts are assigned. For a 9-nucleotide target, 10 level -1 parts are needed."
          "These are divided into two level 0 assemblies each of 5 parts:")
    PPRpartlist = ppr_assort(bp_pprcode, 9)
    gross_part_list.append(PPRpartlist)
    assembly_list = assembly_reformat(gross_part_list)
    for protein in assembly_list:
        assemblynum = 0
        for assembly in protein:
            assemblynum += 1
            print(pd.DataFrame({
                ("Level 0 Assembly number " + str(assemblynum) + "/" + str(
                    len(protein))): assembly,
            }))

    intro()
    return

def spreadsheet_to_cauldron():
    print("Place your spreadsheet inside " + os.getcwd() + " or direct to its location. \n"
       "If you wish to use custom or MoClo parts, place the .gb files inside " + os.getcwd() +
          "\GRASP_Seqs, and use the same names within your spreadsheet. ")
    ssheet = str(input("Spreadsheet name or path: "))
    if ssheet[0] == '"':
        ssheet = ssheet[1:-1]

    with open(ssheet, 'r') as sheetin:
        sheetread = list(csv.reader(sheetin))

    lenrow1 = len(sheetread[0])
    maxrow = max([len(i) for i in sheetread])

    if lenrow1 < maxrow:
        print("ERROR: First row is too short ("+ str(lenrow1) +" vs max " + str(maxrow) + "), assembly may fail. \n"
      "DNAcauldron seems to need empty cells in row 1 equal to (or greater than) the maximum number of cells in the table")

    DNAcauldron_plan_prepper(ssheet)

def intro():
    global datenow
    datenow = datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    print("GRASP Assembly Planner (GAP) v0.5 \n")
    print("Options:" + '\n' + "1: Generate a spreadsheet and Assembly Plan from one or multiple PPR targets using DNAcauldron "
        "(https://edinburgh-genome-foundry.github.io/DnaCauldron/). Use if assembling with GRASP Kit parts only.\n" +
        "2: Generate an Assembly Plan Spreadsheet compatible with DNAcauldron without performing assembly. "
        "Use if you are using customised parts or MoClo parts not included in the GRASP Kit \n" +
        "3: Run DNAcauldron on a pre-existing spreadsheet (e.g. one made in Function 2). \n"
        "4: Print a basic plan from one or multiple PPR targets without assembly. Use to give an idea of how GRASP"
        "Assembly looks. \n" +
        "5: Show an explanation of this program: how a 9-motif binding PPR is created for 9-base target \n" +
        "6/other: Exit the program")
    selection = str(input("Select an option: "))
    if selection == '1':
        makeplan = 1
        GAP_intro_DNAcauldron(makeplan)
    elif selection == '2':
        makeplan = 0
        GAP_intro_DNAcauldron(makeplan)
    elif selection == '3':
        spreadsheet_to_cauldron()
    elif selection == '4':
        GAP_intro()
    elif selection == '5':
        tutorial()
    else:
        exit_opt = str(input("Exit? (Y/N) "))
        if exit_opt in ["N", "n", "no", "No", "NO", "nO"]:
            intro()
        else:
            return

intro()