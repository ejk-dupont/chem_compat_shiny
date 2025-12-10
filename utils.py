import re
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcOxidationNumbers
import pandas as pd
import requests

def rg_list():
   reactivegroups = ["Acetals, Ketals, Hemiacetals, and Hemiketals (70)",
            "Acids, Carboxylic (3)",
            "Acids, Strong Non-oxidizing (1)",
            "Acids, Strong Oxidizing (2)",
            "Acids, Weak (60)",
            "Acrylates and Acrylic Acids (71)",
            "Acyl Halides, Sulfonyl Halides, and Chloroformates (40)",
            "Alcohols and Polyols (4)",
            "Aldehyde (5)",
            "Alkynes, with Acetylenic Hydrogen (63)",
            "Alkynes, with No Acetylenic Hydrogen (64)",
            "Amides and Imides (6)",
            "Amines, Aromatic (68)",
            "Amines, Phosphines, and Pyridines (7)",
            "Anhydrides (37)",
            "Aryl Halides (66)",
            "Azo, Diazo, Azido, Hydrazine, and Azide Compounds (8)",
            "Bases, Strong (10)",
            "Bases, Weak (61)",
            "Carbamates (9)",
            "Carbonate Salts (62)",
            "Chlorosilanes (55)",
            "Conjugated Dienes (65)",
            "Cyanides, Inorganic (11)",
            "Diazonium Salts (25)",
            "Epoxides (34)",
            "Esters, Sulfate Esters, Phosphate Esters, Thiophosphate Esters, and Borate Esters (13)",
            "Ethers (14)",
            "Fluoride Salts, Soluble (48)",
            "Fluorinated Organic Compounds (47)",
            "Halogenated Organic Compounds (17)",
            "Halogenating Agents (59)",
            "Hydrocarbons, Aliphatic Saturated (29)",
            "Hydrocarbons, Aliphatic Unsaturated (28)",
            "Hydrocarbons, Aromatic (16)",
            "Isocyanates and Isothiocyanates (18)",
            "Ketones (19)",
            "Metal Hydrides, Metal Alkyls, Metal Aryls, and Silanes (35)",
            "Metals, Alkali, Very Active (21)",
            "Metals, Elemental and Powder, Active (22)",
            "Metals, Less Reactive (23)",
            "Nitrate and Nitrite Compounds, Inorganic (69)",
            "Nitrides, Phosphides, Carbides, and Silicides (51)",
            "Nitriles (26)",
            "Nitro, Nitroso, Nitrate, and Nitrite Compounds, Organic (27)",
            "Non-Redox-Active Inorganic Compounds (46)",
            "Not Chemically Reactive (98)",
            "Organometallics (42)",
            "Oxidizing Agents, Strong (44)",
            "Oxidizing Agents, Weak (49)",
            "Oximes (75)",
            "Peroxides, Organic (30)",
            "Phenolic Salts (72)",
            "Phenols and Cresols (31)",
            "Polymerizable Compounds (76)",
            "Quaternary Ammonium and Phosphonium Salts (73)",
            "Reducing Agents, Strong (45)",
            "Reducing Agents, Weak (50)",
            "Salts, Acidic (38)",
            "Salts, Basic (39)",
            "Siloxanes (58)",
            "Sulfides, Inorganic (33)",
            "Sulfides, Organic (20)",
            "Sulfite and Thiosulfate Salts (74)",
            "Sulfonates, Phosphonates, and Thiophosphonates, Organic (32)",
            "Thiocarbamate Esters and Salts/Dithiocarbamate Esters and Salts (12)",
            "Water (100)"]
   return reactivegroups

def rg_list_numsort():
    reactivegroups_numsort = ['Acids, Strong Non-oxidizing (1)',
        'Acids, Strong Oxidizing (2)',
        'Acids, Carboxylic (3)',
        'Alcohols and Polyols (4)',
        'Aldehyde (5)',
        'Amides and Imides (6)',
        'Amines, Phosphines, and Pyridines (7)',
        'Azo, Diazo, Azido, Hydrazine, and Azide Compounds (8)',
        'Carbamates (9)',
        'Bases, Strong (10)',
        'Cyanides, Inorganic (11)',
        'Thiocarbamate Esters and Salts/Dithiocarbamate Esters and Salts (12)',
        'Esters, Sulfate Esters, Phosphate Esters, Thiophosphate Esters, and Borate Esters (13)',
        'Ethers (14)',
        'Hydrocarbons, Aromatic (16)',
        'Halogenated Organic Compounds (17)',
        'Isocyanates and Isothiocyanates (18)',
        'Ketones (19)',
        'Sulfides, Organic (20)',
        'Metals, Alkali, Very Active (21)',
        'Metals, Elemental and Powder, Active (22)',
        'Metals, Less Reactive (23)',
        'Diazonium Salts (25)',
        'Nitriles (26)',
        'Nitro, Nitroso, Nitrate, and Nitrite Compounds, Organic (27)',
        'Hydrocarbons, Aliphatic Unsaturated (28)',
        'Hydrocarbons, Aliphatic Saturated (29)',
        'Peroxides, Organic (30)',
        'Phenols and Cresols (31)',
        'Sulfonates, Phosphonates, and Thiophosphonates, Organic (32)',
        'Sulfides, Inorganic (33)',
        'Epoxides (34)',
        'Metal Hydrides, Metal Alkyls, Metal Aryls, and Silanes (35)',
        'Anhydrides (37)',
        'Salts, Acidic (38)',
        'Salts, Basic (39)',
        'Acyl Halides, Sulfonyl Halides, and Chloroformates (40)',
        'Organometallics (42)',
        'Oxidizing Agents, Strong (44)',
        'Reducing Agents, Strong (45)',
        'Non-Redox-Active Inorganic Compounds (46)',
        'Fluorinated Organic Compounds (47)',
        'Fluoride Salts, Soluble (48)',
        'Oxidizing Agents, Weak (49)',
        'Reducing Agents, Weak (50)',
        'Nitrides, Phosphides, Carbides, and Silicides (51)',
        'Chlorosilanes (55)',
        'Siloxanes (58)',
        'Halogenating Agents (59)',
        'Acids, Weak (60)',
        'Bases, Weak (61)',
        'Carbonate Salts (62)',
        'Alkynes, with Acetylenic Hydrogen (63)',
        'Alkynes, with No Acetylenic Hydrogen (64)',
        'Conjugated Dienes (65)',
        'Aryl Halides (66)',
        'Amines, Aromatic (68)',
        'Nitrate and Nitrite Compounds, Inorganic (69)',
        'Acetals, Ketals, Hemiacetals, and Hemiketals (70)',
        'Acrylates and Acrylic Acids (71)',
        'Phenolic Salts (72)',
        'Quaternary Ammonium and Phosphonium Salts (73)',
        'Sulfite and Thiosulfate Salts (74)',
        'Oximes (75)',
        'Polymerizable Compounds (76)',
        'Not Chemically Reactive (98)']
    return reactivegroups_numsort

reactivenumbers = ["70","3","1","2","60","71","40","4","5","63","64","6","68","7","37","66","8","10","61","9","62","55","65","11","25","34","13","14","48","47","17","59","29","28","16","18","19","35","21","22","23","69","51","26","27","46","98","42","44","49","75","30","72","31","76","73","45","50","38","39","58","33","20","74","32","12","100"]

def element_list():
    elementgroups = ['All - alphabetical',
                     'All - numerical',
                     'Acids',
                     'Bases',
                     'Boron Compounds',
                     'Carbon Compounds',
                     'Halogen Compounds',
                     'Metals / Neutral Salts',
                     'Nitrogen Compounds',
                     'Oxidizing Agents',
                     'Oxygen Compounds',
                     'Phosphorus Compounds',
                     'Reducing Agents',
                     'Silicon Compounds',
                     'Sulfur Compounds']
    return elementgroups

def adsorbent_list():
    adsorbentlist=['Absorbent GP (wood pulp) (Cellulosic)',
                'Kimwipes (Cellulosic)',
                'K-Sorb (Cellulosic)',
                'Paper Towels (Cellulosic)',
                'Peat Moss (Cellulosic)',
                'Putzwolle (Cellulosic)',
                'Safety Sorb (Cellulosic)',
                'Sawdust (Cellulosic)',
                'SlikWik (Cellulosic)',
                'Terrasorb L (peat) (Cellulosic)',
                'Lite-Dri (Cellulosic)',
                'Budget Dry (Clay/Mineral)',
                'Kitty Litter (Clay/Mineral)',
                'Safe-T-Sorb (Clay/Mineral)',
                'Vermiculite (Clay/Mineral)',
                'Zorball (Clay/Mineral)',
                'Oil Dri (Clay/Mineral)',
                'Plastic Pigs (Polymeric)',
                'Sand (Sand)',
                'Soil (Dirt/Earth)',
                'Dirt (Dirt/Earth)']
    return adsorbentlist

def determinefunctionals (smilestr): #Automatically determines functional groups
    m = Chem.MolFromSmiles(smilestr)

    #oxidation states
    CalcOxidationNumbers(m)
    ox_data = pd.read_excel("oxidation_states.xlsx")
    ox_list = []
    for i, atom in enumerate(m.GetAtoms()):
        atom.SetProp("atomNote",f"{i}:{atom.GetProp('OxidationNumber')}")
        entry = [atom.GetSymbol(),atom.GetProp('OxidationNumber')]
        ox_list.append(entry)

    #pKa prediction setup
    data_pka = predict_pka(smilestr)
    pd_pka = pd.DataFrame.from_dict(data_pka)
    allacid_pka = (pd_pka['Acid'].dropna()).tolist()
    allbase_pka = (pd_pka['Base'].dropna()).tolist()
    str_pka = " "

    if not allacid_pka:
        acid_pka = 20
    elif len(allacid_pka) > 1:
        acid_pka = float(min(allacid_pka)) #find the lowest pKa (strongest acid)
    else:
        acid_pka = float(allacid_pka[0])
    
    if not allbase_pka:
        base_pka = -20
    elif len(allbase_pka) > 1:
        base_pka = float(max(allbase_pka)) #find the highest pKa (strongest base)
    else:
        base_pka = float(allbase_pka[0])
        
    all_pka = allacid_pka + allbase_pka
    all_pka = [x for x in all_pka if str(x) != 'nan'] #remove duplicates
    all_pka = sorted(set(all_pka)) #sort
    str_pka = ' '.join(all_pka) #string of predicted pKas printed on interface

    #halogen counts
    groups = []
    halo_str = ["F","Cl","Br","I"]
    nothalo_str = ["Fe", "Fm", "Fr", "In", "Ir"] #This is such terrible programming, but I don't care I'm tired and it works.
    halocount = 0
    for halo in halo_str:
        count = smilestr.count(halo)
        halocount += count
    for nothalo in nothalo_str:
        notcount = smilestr.count(nothalo)
        halocount = halocount - notcount
    
    #Let's assign functional groups...
    #Acetals, Ketals, Hemiacetals, and Hemiketals (70) 
    acetal = Chem.MolFromSmarts('[O][CX4][O]')
    if m.HasSubstructMatch(acetal):
        groups.append('Acetals, Ketals, Hemiacetals, and Hemiketals (70)')
        
    #Acids, Carboxylic (3) 
    carboxylic = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if m.HasSubstructMatch(carboxylic):
        groups.append('Acids, Carboxylic (3)')
        
    #Acids, Strong Non-oxidizing (1)
    stracidnonox1 = Chem.MolFromSmarts('[F,Cl,Br,I].[Al,Sb,Fe,Ga,P,Sn,Ti,V,Zn,Zr]')
    if m.HasSubstructMatch(stracidnonox1) or (smilestr in halo_str) or acid_pka < -2:
        groups.append('Acids, Strong Non-oxidizing (1)')

    #Acids, Strong Oxidizing (2)
    stracidox = Chem.MolFromSmarts('[$([*](=O)(=O)[O;H]),$([*](=O)(-O)[O;H]),$([*](=O)(=O)(=O)[O;H])]')
    stracidox3 = Chem.MolFromSmarts('[$(C(=O)[O][O;H]),$([O-][S](=O)(=O)[O-])]')
    if m.HasSubstructMatch(stracidox) or m.HasSubstructMatch(stracidox3):
        groups.append('Acids, Strong Oxidizing (2)')
        
    #Acids, Weak (60)
    phoswkacid = Chem.MolFromSmarts('[$([P](=O)[O]),$([P](=S)[O]),$(C#N)]')
    wkacid = Chem.MolFromSmarts('[$([Se](=O)[O]),$([S](=O)(=O)[O]),$([O][B][O]),$(C#N),$([S][C]=[S,N])]')
    if m.HasSubstructMatch(phoswkacid) or m.HasSubstructMatch(wkacid) or (10 > acid_pka > -2):
        groups.append('Acids, Weak (60)')
    
    #Acrylates and Acrylic Acids (71) (TERMINAL)
    acrylate = Chem.MolFromSmarts('[C;H2]=[C]-[C](=O)-[O]')
    if m.HasSubstructMatch(acrylate):
        groups.append('Acrylates and Acrylic Acids (71)')
        
    #Acyl Halides, Sulfonyl Halides, and Chloroformates (40) eh
    ahalide = Chem.MolFromSmarts('[CX3](=O)[F,Cl,Br,I]')
    sulfonyl = Chem.MolFromSmarts('[Sv6](=O)(=O)[F,Cl,Br,I]')
    if m.HasSubstructMatch(ahalide) or m.HasSubstructMatch(sulfonyl):
        groups.append('Acyl Halides, Sulfonyl Halides, and Chloroformates (40)')
        
    #Alcohols and Polyols (4)
    alcohol = Chem.MolFromSmarts('[CX4;C,H]-[O;H1]')
    if m.HasSubstructMatch(alcohol):
        groups.append('Alcohols and Polyols (4)')
        
    #Aldehyde (5)
    aldehyde = Chem.MolFromSmarts('[CX3;H1](=O)[#6]')
    if m.HasSubstructMatch(aldehyde):
        groups.append('Aldehyde (5)')
    
    #Alkynes, with Acetylenic Hydrogen (63)
    alkyneh = Chem.MolFromSmarts('C#[CH]')
    if m.HasSubstructMatch(alkyneh):
        groups.append('Alkynes, with Acetylenic Hydrogen (63)')
    
    #Alkynes, with No Acetylenic Hydrogen (64)
    alkyne = Chem.MolFromSmarts('[C!H]#[C!H]')
    if m.HasSubstructMatch(alkyne):
        groups.append('Alkynes, with No Acetylenic Hydrogen (64)')
    
    #Amides and Imides (6)
    amide1 = Chem.MolFromSmarts('[n,N][c,C](=O)[*;!O;!S]')
    amide2 = Chem.MolFromSmarts('[n,N][CH](=O)')
    if m.HasSubstructMatch(amide1) or m.HasSubstructMatch(amide2):
        groups.append('Amides and Imides (6)')
    
    #Amines, Aromatic (68)
    amine_aro = Chem.MolFromSmarts('[NX3][c]')
    if m.HasSubstructMatch(amine_aro):
        groups.append('Amines, Aromatic (68)')
    
    #Amines, Phosphines, and Pyridines (7)
    amine = Chem.MolFromSmarts('[Nv3;H2,H1,H0;!$(NC=O);!$(N=C=O);!$(NC=S);!$(N-N);!$(N-c),!$(N=O)]')
    #amine = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O);!$(NC=S);!$(N-N)]')
    phosphine = Chem.MolFromSmarts('[Pv3]')
    pyridine = Chem.MolFromSmarts('cnc')
    if m.HasSubstructMatch(amine) or m.HasSubstructMatch(phosphine) or m.HasSubstructMatch(pyridine):
        groups.append('Amines, Phosphines, and Pyridines (7)')
        
    #Anhydrides (37)
    anhydride = Chem.MolFromSmarts('[*](=[OX1])[OX2][*](=[OX1])')
    if m.HasSubstructMatch(anhydride):
        groups.append('Anhydrides (37)')
        
    #Aryl Halides (66)
    arylhalide = Chem.MolFromSmarts('c-[F,Cl,Br,I]')
    if m.HasSubstructMatch(arylhalide):
        groups.append('Aryl Halides (66)')
        
    #Azo, Diazo, Azido, Hydrazine, and Azide Compounds (8)
    azo = Chem.MolFromSmarts('[#7]=[#7]')
    hydrozine = Chem.MolFromSmarts('[#7][#7]')
    if m.HasSubstructMatch(azo) or m.HasSubstructMatch(hydrozine):
        groups.append('Azo, Diazo, Azido, Hydrazine, and Azide Compounds (8)')
        
    #Bases, Strong (10)
    alkali = Chem.MolFromSmarts('[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba]')
    oxid = Chem.MolFromSmarts('[O-O,O-CH3,Ti-O,Zr-O]')
    nitride = Chem.MolFromSmarts('[N-3,NH2-]')
    aziridine = Chem.MolFromSmarts('C1CN1')
    hydr_zinelamine = Chem.MolFromSmarts('[NH2][NH2,OH]')
    if (((m.HasSubstructMatch(alkali) and (m.HasSubstructMatch(Chem.MolFromSmarts('[O]')) or m.HasSubstructMatch(oxid) or \
        m.HasSubstructMatch(nitride) or m.HasSubstructMatch(amine) or m.HasSubstructMatch(amide1) or m.HasSubstructMatch(amide2))) \
        or m.HasSubstructMatch(aziridine) or m.HasSubstructMatch(hydr_zinelamine)) or base_pka > 12):
        groups.append('Bases, Strong (10)')
        
    #Bases, Weak (61)
    if 2 < base_pka < 12:
        groups.append('Bases, Weak (61)')
        
    #Carbamates (9)
    carbamate = Chem.MolFromSmarts('[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]')
    if m.HasSubstructMatch(carbamate):
        groups.append('Carbamates (9)') 
        
    #Carbonate Salts (62)
    carbsalts = Chem.MolFromSmarts('[O-][#6](=O)[O-]')
    if m.HasSubstructMatch(carbsalts):
        groups.append('Carbonate Salts (62)')
        
    #Chlorosilanes (55)
    chlorosilane = Chem.MolFromSmarts('[Cl]-[Si]')
    if m.HasSubstructMatch(chlorosilane):
        groups.append('Chlorosilanes (55)')
    
    #Conjugated Dienes (65)
    conjdiene = Chem.MolFromSmarts('[#6]=[#6]-[#6]=[#6]')
    if m.HasSubstructMatch(conjdiene):
        groups.append('Conjugated Dienes (65)')
        
    #Cyanides, Inorganic (11)
    cyaninorg = Chem.MolFromSmarts('[*&!#6][#6]#N')
    if m.HasSubstructMatch(cyaninorg):
        groups.append('Cyanides, Inorganic (11)')
        
    #Diazonium Salts (25)
    diasalt = Chem.MolFromSmarts('[N+]#N')
    if m.HasSubstructMatch(diasalt):
        groups.append('Diazonium Salts (25)')
        
    #Epoxides (34)
    epoxide = Chem.MolFromSmarts('C1CO1')
    if m.HasSubstructMatch(epoxide):
        groups.append('Epoxides (34)')   
        
    #Esters, Sulfate Esters, Phosphate Esters, Thiophosphate Esters, and Borate Esters (13)
    ester = Chem.MolFromSmarts('[O&H0][#6](=O)[*&!N&!O&!Cl]')
    sester = Chem.MolFromSmarts('[O][Sv6](=O)(=O)[O]')
    pester = Chem.MolFromSmarts('[O][Pv5](=[S,O])(O)[O]')
    bester = Chem.MolFromSmarts('[Bv3](O)(O)[O]')
    if m.HasSubstructMatch(ester) or m.HasSubstructMatch(sester) or m.HasSubstructMatch(pester) or m.HasSubstructMatch(bester): 
        groups.append('Esters, Sulfate Esters, Phosphate Esters, Thiophosphate Esters, and Borate Esters (13)')
    
    #Ethers (14)
    ether = Chem.MolFromSmarts('[OD2,oD2]([#6;!$(C=O)])[#6;!$(C=O)]')
    if m.HasSubstructMatch(ether) and not m.HasSubstructMatch(epoxide):
        groups.append('Ethers (14)')
        
    #Fluoride Salts, Soluble (48)
    fsalts = Chem.MolFromSmarts('[F-]')
    if m.HasSubstructMatch(fsalts):
        groups.append('Fluoride Salts, Soluble (48)')
    
    #Fluorinated Organic Compounds (47)
    forgo = Chem.MolFromSmarts('[#6][F]')
    if m.HasSubstructMatch(forgo) and not m.HasSubstructMatch(fsalts) :
        groups.append('Fluorinated Organic Compounds (47)')
        
    #Halogenated Organic Compounds (17)
    halide = Chem.MolFromSmarts('[#6][Cl,Br,I]')
    nonorg = Chem.MolFromSmarts('[#7,#14,#15,C-]')
    if m.HasSubstructMatch(halide) and not m.HasSubstructMatch(nonorg):
        groups.append('Halogenated Organic Compounds (17)')
        
    #Halogenating Agents (59)
    agents = Chem.MolFromSmarts('[Al,Sb,As,NH4,B,Br,Ge,P,Mo,Ti,V,Zn,Zr,O]')
    halo = Chem.MolFromSmarts('[F,Cl,Br,I]')
    if (m.HasSubstructMatch(halo) and m.HasSubstructMatch(agents)) or halocount > 1:
        groups.append('Halogenating Agents (59)')
        
    #Hydrocarbons, Aliphatic Saturated (29)
    hcsat = Chem.MolFromSmarts('[CH2]-[CH2]-[CH2]')
    methane = Chem.MolFromSmarts('[CH4]')
    ethane = Chem.MolFromSmarts('[CH3]-[CH3]')
    butane = Chem.MolFromSmarts('[CH3]-[CH2]-[CH3]')
    if m.HasSubstructMatch(hcsat) or m.HasSubstructMatch(methane) or m.HasSubstructMatch(ethane) or m.HasSubstructMatch(butane):
        groups.append('Hydrocarbons, Aliphatic Saturated (29)')
        
    #Hydrocarbons, Aliphatic Unsaturated (28)
    hcaunsat = Chem.MolFromSmarts('[C,H]-[C]=[C]-[*&!O]')
    ethene = Chem.MolFromSmarts('[CH2]=[CH2]')
    hcaunsat3 = Chem.MolFromSmarts('[C]-[C]#[C]')
    if m.HasSubstructMatch(hcaunsat) or m.HasSubstructMatch(ethene) or m.HasSubstructMatch(hcaunsat3):
        groups.append('Hydrocarbons, Aliphatic Unsaturated (28)')
        
    #Hydrocarbons, Aromatic (16)
    haromatic = Chem.MolFromSmarts('c1ccccc1')
    phenol = Chem.MolFromSmarts('c-[OH]')
    if m.HasSubstructMatch(haromatic) and not(m.HasSubstructMatch(amine_aro) or m.HasSubstructMatch(phenol) or m.HasSubstructMatch(arylhalide)):
        groups.append('Hydrocarbons, Aromatic (16)')
    
    
    #Isocyanates and Isothiocyanates (18)
    isocynate = Chem.MolFromSmarts('[*]-N=C=O')
    if m.HasSubstructMatch(isocynate):
        groups.append('Isocyanates and Isothiocyanates (18)')
    
    #Ketones (19)
    #ketone = Chem.MolFromSmarts('[#6][#6](=O)[#6][*&!N&!O]')
    ketone = Chem.MolFromSmarts('[#6][#6](=O)[#6]')
    if m.HasSubstructMatch(ketone):
        groups.append('Ketones (19)')
        
    #!! Metal Hydrides, Metal Alkyls, Metal Aryls, and Silanes (35)
    metalhydride = Chem.MolFromSmarts('[Al,Be,Ca,Cs,Cr,Ca,Cu,Fe,Li,Mg,Ni,Pd,K,Rb,Na,Ti,Tl,Zn][H,C,c]')
    silanes = Chem.MolFromSmarts('[Si][C,c]')
    if m.HasSubstructMatch(metalhydride) or m.HasSubstructMatch(silanes):
        groups.append('Metal Hydrides, Metal Alkyls, Metal Aryls, and Silanes (35)')
        
    #!! Metals, Alkali, Very Active (21)
    #metalkali = Chem.MolFromSmarts('[Li,Na,K,Rb,Cs]')
    metal_alkali = ['Li','Na','K','Rb','Cs']
    metal_alkali_flag = False
    for entry in ox_list:
        if entry[0] in metal_alkali:
            if entry[1] not in list(ox_data[ox_data["Atomic Letter"] == entry[0]]["Charges"].values[0].split(",")):
                metal_alkali_flag = False
    if metal_alkali_flag == True:
        groups.append('Metals, Alkali, Very Active (21)')
        
    #!! Metals, Elemental and Powder, Active (22)
    #metalelemental = Chem.MolFromSmarts('[Al,Sb,Ba,Be,Cd,Ca,Ce,Cr,Co,Ga,Hf,In,Fe,Mg,Ni,Se,Sr,Ti,Ta,Th,U,V,Y,Zn,Zr]')
    #if m.HasSubstructMatch(metalelemental):
    metal_active = ["Al","Sb","Ba","Be","Cd","Ca","Ce","Cr","Co","Ga","Hf","In","Fe","Mg","Ni","Se","Sr","Ti","Zn","Zr"]
    metal_active_flag = False
    for entry in ox_list:
        if entry[0] in metal_active:
            if entry[1] not in list(ox_data[ox_data["Atomic Letter"] == entry[0]]["Charges"].values[0].split(",")):
                metal_active_flag = True
    if metal_active_flag == True:
        groups.append('Metals, Elemental and Powder, Active (22)')
        
    #!! Metals, Less Reactive (23)
    #metalless = Chem.MolFromSmarts('[Ag,As,Au,Cu,Pb,Hg,Pd,Pt,Mn,Mo,Rh,Te,Tl,Sn,W]')
    metal_less = ["Ag","As","Au","Cu","Pb","Hg","Pd","Pt","Mn","Mo","Rh","Te","Tl","Sn","W"]
    metal_less_flag = False
    for entry in ox_list:
        if entry[0] in metal_less:
            if entry[1] not in list(ox_data[ox_data["Atomic Letter"] == entry[0]]["Charges"].values[0].split(",")):
                metal_less_flag = True
    if metal_less_flag == True:
        groups.append('Metals, Less Reactive (23)')
        
    #Nitrate and Nitrite Compounds, Inorganic (69)
    nitrate = Chem.MolFromSmarts('[O-][N+](=O)[O-]')
    nitrite = Chem.MolFromSmarts('[O-][N](=O)')
    metal = Chem.MolFromSmarts('[*+,*++;!N]')
    if (m.HasSubstructMatch(nitrate) or m.HasSubstructMatch(nitrite)) and m.HasSubstructMatch(metal):
        groups.append('Nitrate and Nitrite Compounds, Inorganic (69)')
        
    #Nitrides, Phosphides, Carbides, and Silicides (51)
    nitride = Chem.MolFromSmarts('[N-3,NH2-]')
    phosphide = Chem.MolFromSmarts('[P-3]')
    phosphide2 = Chem.MolFromSmarts('[PH2-]')
    carbides = Chem.MolFromSmarts('[C-]')
    silicide = Chem.MolFromSmarts('[Si-4]')
    if m.HasSubstructMatch(nitride) or m.HasSubstructMatch(phosphide) or m.HasSubstructMatch(phosphide2) or m.HasSubstructMatch(carbides) or m.HasSubstructMatch(silicide):
        groups.append('Nitrides, Phosphides, Carbides, and Silicides (51)')
    
    #Nitriles (26)
    nitrile = Chem.MolFromSmarts('[#6]C#N')
    if m.HasSubstructMatch(nitrile) and not m.HasSubstructMatch(nitrate) and not m.HasSubstructMatch(nitrite):
        groups.append('Nitriles (26)')
    
    #Nitro, Nitroso, Nitrate, and Nitrite Compounds, Organic (27)
    nitro = Chem.MolFromSmarts('[N](=O)')
    if m.HasSubstructMatch(nitro) and not m.HasSubstructMatch(metal):
        groups.append('Nitro, Nitroso, Nitrate, and Nitrite Compounds, Organic (27)')
        
    #Non-Redox-Active Inorganic Compounds (46)
    nonredox = Chem.MolFromSmarts('[F,Cl,Br,I].[Cd,Na,Ba,K,Ca,Mg]')
    if m.HasSubstructMatch(nonredox):
        groups.append('Non-Redox-Active Inorganic Compounds (46)')
        
    #Organometallics (42)
    orgmetal2 = Chem.MolFromSmarts('[Al,Be,Ca,Cs,Cr,Ca,Cu,Fe,Li,Mg,Ni,Pd,K,Rb,Na,Ti,Tl,Zn].[C,c]')
    if m.HasSubstructMatch(orgmetal2):
        groups.append('Organometallics (42)')
        
    #Oxidizing Agents, Strong (44)
    stroxidze = Chem.MolFromSmarts('[Pb4+]')
    oxdagent = Chem.MolFromSmarts('[$([O-][Cr](=O)),$([O-][Cl](=O)),$([O-][Br](=O)),$([O-][Mn](=O))]')
    peroxide = Chem.MolFromSmarts('[C](=O)[O][O,H]')
    if m.HasSubstructMatch(oxdagent) or m.HasSubstructMatch(peroxide) or m.HasSubstructMatch(stroxidze) or halocount > 1:
        groups.append('Oxidizing Agents, Strong (44)')
        
    #Oxidizing Agents, Weak (49)
    wkoxdagent1 = Chem.MolFromSmarts('[N,n][O;H]')
    wkoxdagent2 = Chem.MolFromSmarts('N(=O)[N,O]')
    wkoxdagent3 = Chem.MolFromSmarts('[$([N]=[N]=O),$(O=O),$(O=C=O)]')
    wkoxdagent4 = Chem.MolFromSmarts('[$([O-][As]([O-])([O-])=O),$([O-][S]([O-])(=O)=O),$([O-][Se]([O-])(=O)=O)]')
    if m.HasSubstructMatch(wkoxdagent1) or m.HasSubstructMatch(wkoxdagent2) or m.HasSubstructMatch(wkoxdagent3) or m.HasSubstructMatch(wkoxdagent4):
        groups.append('Oxidizing Agents, Weak (49)')
        
    #Oximes (75)
    oximes = Chem.MolFromSmarts('[o,O]-[n,N]=[c,C]')
    if m.HasSubstructMatch(oximes):
        groups.append('Oximes (75)') 
    
    #Peroxides, Organic (30)
    peroxide = Chem.MolFromSmarts('C-O-O-C')
    hperoxide = Chem.MolFromSmarts('[O;H]-[O;H]')
    if m.HasSubstructMatch(peroxide) or m.HasSubstructMatch(hperoxide):
        groups.append('Peroxides, Organic (30)')
        
    #Phenolic Salts (72)
    phsalt = Chem.MolFromSmarts('[*+].[O-][c,C]')
    if m.HasSubstructMatch(phsalt):
        if 'Hydrocarbons, Aromatic' in groups: groups.remove('Hydrocarbons, Aromatic')
        groups.append('Phenolic Salts (72)')
        
    #Phenols and Cresols (31)
    phenol = Chem.MolFromSmarts('c-[OH]')
    if m.HasSubstructMatch(phenol):
        if 'Hydrocarbons, Aromatic' in groups: groups.remove('Hydrocarbons, Aromatic')
        groups.append('Phenols and Cresols (31)')
        
    #!!!Polymerizable Compounds (76)
    butadiene = Chem.MolFromSmarts('C=[CX3;H1]-[CX3;H1]=C')
    styrene = Chem.MolFromSmarts('[C,R6]-C=C')
    others = Chem.MolFromSmarts('[$(C=CCN),$(O=C-C=O)]')
    #propylene_oxide = Chem.MolFromSmarts('[C,R3]-[O,R3]-[C,R3]-C') #?????
    if m.HasSubstructMatch(butadiene) or m.HasSubstructMatch(styrene) or m.HasSubstructMatch(others):
        groups.append('Polymerizable Compounds (76)')
        
    #Quaternary Ammonium and Phosphonium Salts (73)
    quadamm = Chem.MolFromSmarts('[NX4+,nX4+][*&!O]')
    quadphos = Chem.MolFromSmarts('[PX4+,pX4+][*&!O]')
    if m.HasSubstructMatch(quadamm) or m.HasSubstructMatch(quadphos):
        groups.append('Quaternary Ammonium and Phosphonium Salts (73)')
        
    #Reducing Agents, Strong (45)
    streduce = Chem.MolFromSmarts('[Li+,K+,Ba2+,Ca2+,Na+,Mg2+,Al3+]')
    if m.HasSubstructMatch(streduce):
        groups.append('Reducing Agents, Strong (45)')
        
    #Reducing Agents, Weak (50)
    wkreduce = Chem.MolFromSmarts('[Zn2+,Fe2+,Cr3+,Ni2+,Sn2+]')
    if m.HasSubstructMatch(wkreduce):
        groups.append('Reducing Agents, Weak (50)')
        
    #Salts, Acidic (38)
    acidsalt = Chem.MolFromSmarts('[H+].[F-,Cl-,Br-,I-]')
    acidsalt2 = Chem.MolFromSmarts('[$([O-,O][S](=O)(=O)[O-]),$([O-,O][P](=O)(O)[O-])]')
    if m.HasSubstructMatch(acidsalt) or m.HasSubstructMatch(acidsalt2):
        groups.append('Salts, Acidic (38)')
        
    #Salts, Basic (39)
    basicsalt = Chem.MolFromSmarts('[C](=O)[O-]')
    basicsalt2 = Chem.MolFromSmarts('[$([O-,O][C](=O)[O-]),$([C]#[N-]),$([O-;H1])]')
    if m.HasSubstructMatch(basicsalt) or m.HasSubstructMatch(basicsalt2):
        groups.append('Salts, Basic (39)')
        
    #Siloxanes (58)
    silox = Chem.MolFromSmarts('[Si][O]')
    silox2 = Chem.MolFromSmarts('[Si]=[O]')
    if m.HasSubstructMatch(silox) or m.HasSubstructMatch(silox2):
        groups.append('Siloxanes (58)')
    
    #Sulfides, Inorganic (33)
    sulfideinorg = Chem.MolFromSmarts('[S;!$(S=O)]')
    if m.HasSubstructMatch(sulfideinorg):
        groups.append('Sulfides, Inorganic (33)')
        
    #Sulfides, Organic (20)
    sulfideorg = Chem.MolFromSmarts('[S,s;!$(S=O)][c,C,H]')
    sulfideorg2 = Chem.MolFromSmarts('C=S')
    if m.HasSubstructMatch(sulfideorg) or m.HasSubstructMatch(sulfideorg2):
        if '33' in groups: groups.remove('33')
        groups.append('Sulfides, Organic (20)')
    
    #Sulfite and Thiosulfate Salts (74)
    sulfite = Chem.MolFromSmarts('[O][SX3](=O)O')
    thiosulfate = Chem.MolFromSmarts('[O][S](=S)(=O)[O]')
    if m.HasSubstructMatch(sulfite) or m.HasSubstructMatch(thiosulfate):
        groups.append('Sulfite and Thiosulfate Salts (74)')
        
    #Sulfonates, Phosphonates, and Thiophosphonates, Organic (32)
    sulfone = Chem.MolFromSmarts('[$([!O][S](=O)(=O)[!O]),$([O;!H][S](=O)(=O)[!O]),$([S]([O-])([O-])[!O])]')
    phosphonate = Chem.MolFromSmarts('[PX4](=[OX1])(-[OX2])(-[OX2])')
    thiophosphonate = Chem.MolFromSmarts('[PX4](=[S])(-[OX2])(-[OX2])')
    if (m.HasSubstructMatch(sulfone) or m.HasSubstructMatch(phosphonate) or m.HasSubstructMatch(thiophosphonate)):
        groups.append('Sulfonates, Phosphonates, and Thiophosphonates, Organic (32)')
                                 
    #Thiocarbamate Esters and Salts/Dithiocarbamate Esters and Salts (12)
    thiocarb = Chem.MolFromSmarts('[S,O][C](=[S,O])[N]')   
    thio = Chem.MolFromSmarts('[S,s]')   
    if m.HasSubstructMatch(thiocarb) and m.HasSubstructMatch(thio) :
        groups.append('Thiocarbamate Esters and Salts/Dithiocarbamate Esters and Salts (12)')
                                 
    #Water (100)
    water = Chem.MolFromSmarts('[OX2;H2]')
    if m.HasSubstructMatch(water):
        groups.append('Water (100)')
    
    return groups, str_pka

def predict_pka(smi): #using https://xundrug.cn/molgpka pKa prediction by graph-convolutional neural network
    #https://pubs.acs.org/doi/10.1021/acs.jcim.1c00075
    upload_url=r'http://xundrug.cn:5001/modules/upload0/'
    param={"Smiles" : ("tmg", smi)}
    headers={'token':'O05DriqqQLlry9kmpCwms2IJLC0MuLQ7'}
    response=requests.post(url=upload_url, files=param, headers=headers)
    jsonbool=int(response.headers['ifjson'])
    if jsonbool==1:
        res_json=response.json()
        if res_json['status'] == 200:
            pka_datas = res_json['gen_datas']
            return pka_datas
        else:
            raise RuntimeError("Error for pKa prediction")
    else:
        raise RuntimeError("Error for pKa prediction")

def otherdata_search(ask_string):
    df_data = pd.read_excel('other_database.xlsx',index_col=0) 
    
    if not re.match('[\d/-]+$', ask_string):
        results_rg = (df_data[df_data['OfficialChemicalName'] == ask_string.upper()]['CRW Reactivity Group(s)'])
        results_cas =  (df_data[df_data['OfficialChemicalName'] == ask_string.upper()]['CAS'])
        results_name = ask_string
    else:
        results_rg = (df_data[df_data['CAS'] == ask_string]['CRW Reactivity Group(s)'])
        results_name = (df_data[df_data['CAS'] == ask_string]['OfficialChemicalName'])
        results_cas = ask_string

    if len(results_rg) == 0:
        return ("None","None","None")
    else:
        if not re.match('[\d/-]+$', ask_string):
            results_cas = results_cas.iloc[0]
        else:
            results_name = results_name.iloc[0]

        results_rgs = list(results_rg.iloc[0].replace(" ","").split(","))
        rgs = []
        reactivegroups = rg_list()
        for num in results_rgs:
            index = reactivenumbers.index(num)
            rgs.append(reactivegroups[index])

        return results_name, results_cas, rgs

def mixlist_add(mixlist,name,cas,fg,smiles): #add to mixture list
    fg_sort = tuple(sorted(fg))
    mixlist_newrow = pd.DataFrame.from_records({'Name':[name],
                        'CAS':[cas],
                        'Reactive Groups':[fg_sort],
                        'SMILES':[smiles],
                        })
    mixlist = pd.concat([mixlist,mixlist_newrow],axis=0,ignore_index=True)
    return mixlist

def mixlist_bulkadd(mixlist,df_bulk): #bulk add to mixture list
    rg_all = []
    reactivegroups = rg_list()

    for i,row in df_bulk.iterrows():
        rg_str = []
        rg_nums = [x.strip() for x in str(row['Reactive Group Numbers']).split(',')]
        for num in rg_nums:
            index = reactivenumbers.index(num)
            rg_str.append(reactivegroups[index])
        rg_all.append(tuple(rg_str))

    df_bulk2 = df_bulk.drop(['Reactive Group Numbers'],axis=1)
    df_bulk2.insert(2,"Reactive Groups", rg_all)
    mixlist = pd.concat([mixlist,df_bulk2],axis=0,ignore_index=True)
    return mixlist

def mixlist_remove(mixlist,index_i): #remove from mixture list
    mixlist2 = mixlist.drop(mixlist.index[index_i])
    return mixlist2.reset_index(drop=True)

def mixlist_moveup(mixlist,index_i):
    df_holder = mixlist.iloc[index_i].copy()
    mixlist.iloc[index_i] = mixlist.iloc[index_i -1].copy()
    mixlist.iloc[index_i -1] = df_holder
    return mixlist.reset_index(drop=True)

def mixlist_movedown(mixlist,index_i):
    df_holder = mixlist.iloc[index_i].copy()
    mixlist.iloc[index_i] = mixlist.iloc[index_i +1].copy()
    mixlist.iloc[index_i +1] = df_holder
    return mixlist.reset_index(drop=True)

def backgroud(cell_value): #Display compat chart with color
        N_color = 'background-color: #FF6347'
        C_color = 'background-color: #FFD700'
        Y_color = 'background-color: #228B22'
        X_color = 'background-color: lightgray'
        no_color = 'background-color: white'
        
        if cell_value == "N":
            return N_color
        elif cell_value =="C":
            return C_color
        elif cell_value =="SR":
            return C_color
        elif cell_value =="Y":
            return Y_color
        elif cell_value =="X":
            return X_color
        else:
            return no_color

def makechart(dfinputs): #make compat chart
    rg_list = list(dfinputs["Reactive Groups"])
    flat_rgs = [group for groups in rg_list for group in groups]
    flat_rgs = [*set(flat_rgs)] #remove duplicates
    flat_rgs = sorted(flat_rgs)
    dfcrw4 = pd.read_excel('CRW4_Master.xlsx', sheet_name="CRW4Chart",index_col=0) 
    dfcrw4_comments = pd.read_excel('CRW4_Master.xlsx', sheet_name="CRW4Comments",index_col=0) 

    dffunctchart = pd.DataFrame(data="", index=flat_rgs, columns=flat_rgs) #make empty dataframe with columns and rows with functional groups
    for rowIndex, row in dffunctchart.iterrows(): #Fill smaller chemical compatability matrix
        for columnIndex, value in row.items():
            dffunctchart.loc[rowIndex, columnIndex] = dfcrw4.loc[rowIndex,columnIndex]

    dfchart = pd.DataFrame(data ="", columns = dfinputs['Name'].tolist(), index = dfinputs['Name'].tolist())
    dfchart_comments = pd.DataFrame(data ="", columns = dfinputs['Name'].tolist(), index = dfinputs['Name'].tolist())

    for index, row in dfinputs.iterrows(): #loop through chemicals (1)
        for index2, row2 in dfinputs.iterrows(): #loop through chemicals (2) again - to get all combos of chemical 1 & chemical 2
            if index2 < index: #skip other side of diagonal
                continue
            else:
                #print(dfinputs['Name'][index])
                #print(dfinputs['Name'][index2])
                list_comments = [] #keeps track of comments or message
                list_compatability = [] #keeps track of all compatabilities for a single chemical - chemical combination
                for functional in dfinputs['Reactive Groups'][index]: #loop through functional groups for chemical (1)
                    for functional2 in dfinputs['Reactive Groups'][index2]: #loop through functional groups for chemical (2)
                        compatability = []
                        if pd.isnull(dffunctchart.loc[functional][functional2]): #if empty in master, try other combination
                            compatability = (dffunctchart.loc[functional2][functional])
                            list_compatability.append(compatability)
                        else:
                            compatability = (dffunctchart.loc[functional][functional2])
                            list_compatability.append(compatability)

                #        if compatability == "C" or compatability == "SR" or compatability == "N": #only copy message for reactive combinations
                        if pd.isnull(dfcrw4_comments.loc[functional,functional2]): #if empty in master, try other combination
                            comments = (dfcrw4_comments.loc[functional2][functional])
                            list_comments.append(comments)
                            list_comments.append("\n --- \n")
                        else:
                            comments = (dfcrw4_comments.loc[functional][functional2])
                            list_comments.append(comments)
                            list_comments.append("\n --- \n")
                #            message = dfcrw4_comments.loc[functional,functional2]
                #                list_message.append(message)
                        
            if index == index2: #special case for the diagonals
                if "SR" in list_compatability:
                    dfchart.loc[dfinputs['Name'][index2],dfinputs['Name'][index]]=("SR")
                else:
                    dfchart.loc[dfinputs['Name'][index2],dfinputs['Name'][index]]=("X")
            elif "N" in list_compatability: #if there's any N reactivity between functional groups of chem(1) and chem(2), this takes priority
                dfchart.loc[dfinputs['Name'][index2],dfinputs['Name'][index]]=("N")
            elif "C" in list_compatability or "SR" in list_compatability: #if there's any C or SR reactivity between functional groups of chem(1) and chem(2), this takes 2nd priority
                dfchart.loc[dfinputs['Name'][index2],dfinputs['Name'][index]]=("C")
            else:
                dfchart.loc[dfinputs['Name'][index2],dfinputs['Name'][index]]=("Y") #otherwise, all chemistry should be compatible
            #comments added
            dfchart_comments.loc[dfinputs['Name'][index2],dfinputs['Name'][index]]=(list_comments)

    return dfchart_comments, dfchart, dfchart.style.applymap(backgroud)\
        .set_table_styles(
            [{"selector":"td, th", "props":[("border","1px solid grey !important")]},
             {"selector": "th.col_heading","props":[("writing-mode", "vertical-rl"),
                          ('transform', 'rotateZ(180deg)'), 
                          ('verticle-align', 'top'),
                          ('text-align','left')]
                },
            ]
        )\
        .set_properties(**{'text-align':'center'})

def makeadsorbchart(dfinputs, ads): #make adsorbent chart
    dfads = pd.read_excel('CRW4_Master.xlsx', sheet_name="Adsorbs",index_col=0)
    adslist =  list(ads)
    df_adschart = pd.DataFrame(data="",columns= adslist, index = dfinputs['Name'].tolist())
    
    index = 0
    for rowIndex, row in df_adschart.iterrows():
        rgs = list(dfinputs['Reactive Groups'][index])
        for adsorbent in adslist:
            result_list = []
            for rg in rgs:
                entry = dfads[adsorbent][rg]
                result_list.append(str(entry))
            
            if "N" in result_list:
                result = "N"
            elif "C" in result_list:
                result = "C"
            else:
                result = "Y"
            df_adschart [adsorbent][dfinputs['Name'][index]] = result
        index = index + 1

    df_adschart2 = df_adschart.T
    return df_adschart2.style.applymap(backgroud)\
        .set_table_styles(
            [{"selector":"td, th", "props":[("border","1px solid grey !important")]},
             {"selector": "th.col_heading","props":[("writing-mode", "vertical-rl"),
                          ('transform', 'rotateZ(180deg)'), 
                          ('verticle-align', 'top'),
                          ('text-align','left')]
                },
            ]
        )\
        .set_properties(**{'text-align':'center'})

def rg_details(rg):
    dfcrw4 = pd.read_excel('CRW4_Master.xlsx', sheet_name="ReactiveGroups",index_col=0)
    rgnum =  int(rg[rg.find("(")+1:rg.find(")")])
    return dfcrw4.loc[rgnum, :].values.flatten().tolist()
    
def newrg_list(catg):
    dfcrw4 = pd.read_excel('CRW4_Master.xlsx', sheet_name="ReactiveGroups", index_col=0)
    newrg_list=[]
    for columnIndex, value in dfcrw4["Category"].items():
        value_list = str(value).split(", ")
        if str(catg) in value_list:
            newrg_list.append(str(dfcrw4["Reactive Group"][columnIndex])+" (" + str(columnIndex) + ")")
    return newrg_list

def download_chart(df):
    df2 = df.iloc[:,:-1]
    rg_numbers = []
    for index, row in df2.iterrows():
        st_rg = "".join(map(str,row['Reactive Groups']))
        substrings = []
        split_str = st_rg.split("(")
        for s in split_str[1:]:
            split_s = s.split(")")
            if len(split_s) > 1:
                substrings.append(split_s[0])
            substrings_str = ", ".join(substrings)
        rg_numbers.append(substrings_str)
    df3 = df2.drop(['Reactive Groups'],axis=1)
    df3.insert(2,'Reactive Group Numbers', rg_numbers)

    df_chart = makechart(df)
    with pd.ExcelWriter('chart.xlsx') as excel_file:
        df3.to_excel(excel_file, sheet_name="Mixture", index=False)
        df_chart[2].to_excel(excel_file, sheet_name= "Compat Chart")

def listtostr(lst):
    if type(lst) is list: # apply conversion to list columns
        return";".join(map(str, lst))
    else:
        return lst